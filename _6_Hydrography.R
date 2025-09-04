###########################################################
## MATCH BATHYMETRY and HYDOGRAPHY
###########################################################

# Load CMEMS hydrography -- 3D
oce3d <- tidync("../env/cmems_mod_glo_phy_my_0.083deg_P1M.nc") %>%
  hyper_tibble(select_var = c("so","thetao","uo","vo")) %>%
  mutate(across(c(longitude, latitude, depth), as.numeric),
         lon = longitude, lat = latitude, date = time) %>%
  dplyr::select(lon, lat, depth, date, so, thetao, uo, vo) 

# Load CMEMS hydrography -- 2D
oce2d <- tidync("../env/cmems_mod_glo_phy_my_0.083deg_P1M.nc") %>%
  activate("D2,D1,D0") %>%  
  hyper_tibble(select_var = c("bottomT","mlotst","siconc","sithick")) %>%
  rename(lon = longitude, lat = latitude, date = time) %>%
  mutate(across(c(lon, lat), as.numeric)) 

# Load CMEMS chlorophyll
chl <- tidync("../env/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc") %>%
  activate("D2,D1,D0") %>%  
  hyper_tibble(select_var = c("CHL")) %>%
  rename(lon = longitude, lat = latitude, date = time) %>%
  mutate(across(c(lon, lat), as.numeric)) 

#############################

# Join OCE
oce.full <- merge(
  oce3d, oce2d, by=c("lon","lat","date"), all.x=T) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

# OCE averages
oce.avg <- oce.full %>%
  group_by(lat, lon, depth) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm=T)), .groups="drop") %>%
  st_as_sf(coords = c("lon","lat"), crs=4326, remove=F) %>%
  st_transform(32722) %>%
  rename(lon_oce = lon, lat_oce = lat,
         depth_oce = depth) %>%       
  filter(!is.na(lon_oce), !is.na(lat_oce), !is.na(depth_oce)) %>%
  mutate(
    x_oce = st_coordinates(.)[,1],          
    y_oce = st_coordinates(.)[,2])

# CHL averages
chl.avg <- chl %>%
  group_by(lat, lon) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm=T)), .groups="drop") %>%
  st_as_sf(coords = c("lon","lat"), crs=4326, remove=F) %>%
  st_transform(32722) %>%
  rename(lon_chl = lon, lat_chl = lat) %>%       
  filter(!is.na(lon_chl), !is.na(lat_chl)) %>%
  mutate(
    x_chl = st_coordinates(.)[,1],          
    y_chl = st_coordinates(.)[,2])

# Build matrix of coords
oce.mat <- cbind(oce.avg$x_oce, oce.avg$y_oce)
chl.mat <- cbind(chl.avg$x_chl, chl.avg$y_chl)

# kNN search: find nearest chl point for each oce cell
knn.chl <- get.knnx(chl.mat, oce.mat, k = 1)

# add CHL value + distance
oce.avg$CHL <- chl.avg$CHL[knn.chl$nn.index[,1]]
oce.avg$dist_to_chl <- knn.chl$nn.dist[,1]

# Prepare PS data, ensure right coord
ps118_sf <- ps118_combined %>%
  filter(!is.na(x), !is.na(y), !is.na(depth)) %>%
  st_as_sf(coords = c("x","y"), crs=32722) %>%
  mutate(
    depth_ps = abs(depth),             # <-- rename to depth_ps
    x_ps = st_coordinates(.)[,1],
    y_ps = st_coordinates(.)[,2])

# Prepare for matching
oce.mat <- cbind(st_coordinates(oce.avg), depth = oce.avg$depth_oce)
ps.mat  <- cbind(st_coordinates(ps118_sf), depth = ps118_sf$depth_ps)

# Nearest neighbor search (k=3) 
knn.oce <- get.knnx(oce.mat, ps.mat, k=3)

# Outcome
oce.match <- data.table(
  id_ps   = rep(seq_len(nrow(ps.mat)), each=3),
  id_oce  = as.vector(t(knn.oce$nn.index)),
  dist_3D = as.vector(t(knn.oce$nn.dist)),
  x_ps    = ps.mat[rep(seq_len(nrow(ps.mat)), each=3), 1],
  y_ps    = ps.mat[rep(seq_len(nrow(ps.mat)), each=3), 2],
  depth_ps = ps.mat[rep(seq_len(nrow(ps.mat)), each=3), 3],
  x_oce   = oce.mat[as.vector(t(knn.oce$nn.index)), 1],
  y_oce   = oce.mat[as.vector(t(knn.oce$nn.index)), 2],
  depth_oce = oce.mat[as.vector(t(knn.oce$nn.index)), 3])

# Compute distances
oce.match[, dist_horiz := sqrt((x_oce - x_ps)^2 + (y_oce - y_ps)^2)]
oce.match[, dist_vert  := abs(depth_oce - depth_ps)]

# Check outcome
oce.match[, .(
  horiz_min = min(dist_horiz),
  horiz_max = max(dist_horiz),
  horiz_median = median(dist_horiz),
  vert_min = min(dist_vert),
  vert_max = max(dist_vert),
  vert_median = median(dist_vert),
  d3_min = min(dist_3D),
  d3_max = max(dist_3D),
  d3_median = median(dist_3D))]
# horiz_median 3214m -- vert_median 323m

# Filter 
oce.filtered <- oce.match[
  dist_vert <= 100 & dist_horiz <= 10000][
  order(dist_3D), .SD[1], by = id_ps][, `:=`(
    retained_count = .N,
    total_count = nrow(oce.match),
    retained_pct = round(100 * .N / nrow(oce.match), 1))]

# Summary 
oce.filtered[, .(
  retained_count = unique(retained_count),
  total_count = unique(total_count),
  retained_pct = unique(retained_pct))]
# 14% retained

# For final merge: convert to data.table 
setDT(ps118_predict)  
setDT(oce.avg)      

# Append OCE data
oce.filtered <- oce.avg[
  oce.filtered, on =.(x_oce, y_oce), mult="first"]

# Append PS data
ps118_hydro <- ps118_predict[
  oce.filtered, on =.(x = x_ps, y = y_ps), mult="first"]


###########################################################
## LINKS BETWEEN CORAL, HYDRO + BATHY
###########################################################

# Subsample for plotting
set.seed(123)
pca.hydro <- ps118_hydro[sample(.N, min(.N, 1e4))]
#pca.hydro <- ps118_final[seq(1, .N, length.out = min(.N, 1e4))]

# Cap coral densities for plotting only
pca.hydro[, corals_capped := pmin(corals_predicted, 20)]

# Compute PCA on hydro variables only
hydroVars <- c("thetao","uo","vo","bottomT","siconc","mlotst","CHL")
pcs <- prcomp(as.matrix(pca.hydro[, ..hydroVars]), scale.=T)
pca.hydro[, c("PC1","PC2") := .(pcs$x[,1], pcs$x[,2])]

# K-means clustering on PCA scores
pca.hydro[, cluster := kmeans(pca.hydro[, .(PC1, PC2)], centers=4)$cluster]

# Prepare hydro variable arrows for biplot
rotate <- pcs$rotation[,1:2]  # loadings
arrow <- data.table(
  x=0, y=0,
  xend=rotate[,1]*3, yend=rotate[,2]*3,
  HydroVar=rownames(rotate))
label <- data.table(
  x=rotate[,1]*3.3, y=rotate[,2]*3.3,
  label=rownames(rotate))

# Include slope as an arrow (correlation with PCs)
slopeCor <- cor(pca.hydro$PS118_25m_slope, pca.hydro[, .(PC1, PC2)])

# Add to labels
arrow <- rbind(arrow, data.table(
  x=0, y=0,
  xend=slopeCor[1]*3, yend=slopeCor[2]*3, 
  HydroVar="slope"))
label <- rbind(label, data.table(
  x=slopeCor[1]*3.3, y=slopeCor[2]*3.3, 
  label="slope"))

# Compute median coral per cluster
median <- pca.hydro[, .(median = median(corals_capped)), by=cluster]

################################################

## PLOT
# PCA overlaid with densities & slope-corr
merge(pca.hydro, median, by="cluster") %>% ggplot() + 
  geom_point(aes(
    x=PC1, y=PC2, fill=median, shape=factor(cluster)), 
    alpha=0.6, size=3, stroke=0.05) +
  geom_segment(
    data=arrow, aes(x=x, y=y, xend=xend, yend=yend),
    arrow=arrow(length=unit(0.3,"cm")), color="gray55",
    linewidth=0.7) +
  geom_text(
    data=label,
    aes(x=x, y=y, label=label),
    size=3) +
  scale_fill_gradientn(
    #colors=c("aliceblue","springgreen2","yellow","orange2","mediumorchid3","gray22"),
    colors=c("aliceblue","palegreen2","gold","mediumorchid3","gray28"),
    values=scales::rescale(c(0,1.6,1.8,2.1,3)),
    breaks=c(0,1,2,3),
    limits=c(0.08,3),
    name="Median coral density") +
  scale_shape_manual(
    values=c(21,22,23,24), 
    name="Cluster") +
  xlab(paste0("PC1 (", round(100*summary(pcs)$importance[2,1],1),"%)")) +
  ylab(paste0("PC2 (", round(100*summary(pcs)$importance[2,2],1),"%)")) +  
  theme_classic() +
  theme(
    legend.position="right",
    axis.ticks = element_blank(),
    panel.grid.major = element_line(linewidth=0.2)) 

## FIGURE 9b
# Overarching links between hydroVars, slope and density
merge(pca.hydro, median, by="cluster") %>% melt(
  id.vars = c("cluster","PS118_25m_slope","median"),
  measure.vars = hydroVars,
  variable.name = "HydroVar",
  value.name = "Value") %>%
  mutate(slope_bin = cut(PS118_25m_slope, breaks=444)) %>%
  group_by(cluster, HydroVar, slope_bin, median) %>%
  summarise(
    Value = mean(Value), PS118_25m_slope = mean(PS118_25m_slope), .groups="drop") %>%
  filter(HydroVar %in% c("bottomT","mlotst","siconc")) %>%
  ggplot() +
  geom_point(aes(x=PS118_25m_slope, y=Value, color=median), size=1, alpha=0.8) +
  facet_wrap(~HydroVar, scales="free_y", nrow=1) +
  scale_color_gradientn(
    #colors=c("aliceblue","darkseagreen3","yellow","purple3","gray34"),
    # colors=c("aliceblue","springgreen2","yellow","orange2","mediumorchid3","gray22"),
    # values=scales::rescale(c(0,1.6,1.8,2,3)),
    colors=c("aliceblue","palegreen2","gold","mediumorchid3","gray28"),
    values=scales::rescale(c(0,1.6,1.8,2.1,3)),
    breaks=c(0,1,2,3),
    limits=c(0.08,3),
    name="Median coral") +
  geom_smooth(aes(
    x=PS118_25m_slope, y=Value), 
    method="loess", se=F, color="black") +
  labs(x="Slope (°)", y="Hydro variable value") +
  theme_classic() +
  theme(
    legend.position="right",
    panel.grid.major = element_line(linewidth=0.2))

################################################

## STATS
# Densities vs PCs
cor(pca.hydro$corals_capped, pca.hydro[, .(PC1, PC2)])

# Loadings
pcs$rotation[,1:2]

# Different densities per cluster?
pca.hydro %>%
  dunn_test(corals_predicted ~ cluster, p.adjust.method="BH") %>%
  add_significance() 

# Effect vs sample size; ICQ
pca.hydro %>%
  kruskal_effsize(corals_predicted ~ cluster)
pca.hydro %>% group_by(cluster) %>% summarise(
  median = median(corals_predicted), IQR = IQR(corals_predicted))

# Different slopes between clusters?
pca.hydro %>%
  dunn_test(PS118_25m_slope ~ cluster, p.adjust.method="BH") %>%
  add_significance() 

# Effect vs sample size; ICQ
pca.hydro %>%
  kruskal_effsize(PS118_25m_slope ~ cluster)
pca.hydro %>% group_by(cluster) %>% summarise(
  median = median(PS118_25m_slope),
  IQR = IQR(PS118_25m_slope))


###########################################################
 ## BIODIVERSITY HOTSPOTS
###########################################################

# Select organism densities 
predVars <- ps118_predict %>% dplyr::select(
  corals_predicted, demosponges_predicted, glassSponges_predicted, crinoids_predicted) %>%
  scale()

# Run PCA
pca.bio <- prcomp(predVars, scale.=T)

# Assign PC1 as "hotspot index"
ps118_predict$hotspot_index <- pca.bio$x[,1]

# Grid; average index
hotspot.grid <- ps118_predict %>%
  mutate(xr = round(x/500)*500,
         yr = round(y/500)*500) %>%
  group_by(xr, yr) %>%
  summarise(hotspot_index = mean(hotspot_index), .groups="drop") %>%
  filter(hotspot_index > 1)

# Spatial clustering  
clust <- hotspot.grid %>% 
  select(xr, yr) %>% 
  as.matrix() %>% 
  dbscan(eps=500, minPts=5)

# Assign clusters
hotspot.grid$cluster <- clust$cluster

# KNN match of coordinates
knn.hot <- get.knnx(
  data = hotspot.grid %>% select(xr, yr),
  query = ps118_predict %>% select(x, y),  
  k = 1)$nn.index

# Assign clusters to original data
ps118_predict$cluster <- hotspot.grid$cluster[knn.hot]

# Identify true hotspots
hottest <- ps118_predict %>%
  group_by(cluster) %>%
  summarise(
    n_points = n(),
    mean_hotspot = mean(hotspot_index, na.rm=T)) %>%
  filter(cluster != 0, n_points >= 500, mean_hotspot > 1.5)

# Assign hotspot flags
ps118_predict$hotspot_cat <- ifelse(
  ps118_predict$cluster %in% hottest$cluster,
  "Hotspot","Non-Hotspot")

# Create final hydro table
ps118_hydro <- merge(
  ps118_hydro,
  ps118_predict[, c("x","y","hotspot_index","cluster","hotspot_cat")],
  by = c("x","y"),
  all.x = TRUE)

#############################

## STATS
# Hotspot fraction
mean(ps118_predict$hotspot_cat =="Hotspot", na.rm=T) #6%

# Hotspot depth 
mean(ps118_predict$depth[ps118_predict$hotspot_cat == "Hotspot"], na.rm=T) #1672m

# Hotspot extent: points *25 (=resolution) /1000000 (=sqkm)
ps118_predict %>%
  group_by(hotspot_cat) %>%
  summarise(n_points = n()) %>%
  mutate(area = n_points*25/1000000)  # 40 sqkm

#############################

## PLOT
# Reload bathy data
ps118_bathy <- rast("../env/PS118_25m_bathy.img")
ps118_slope <- rast("../env/PS118_25m_slope.tif")

# Basemap
hotspot.rast <- rast(ps118_bathy)

# Write hotspot index into aligned raster
hotspot.rast[cellFromXY(hotspot.rast, ps118_predict[, c("x","y")])] <- 
  ifelse(ps118_predict$hotspot_cat == "Hotspot", 1, 0)

# Plot true hotspots
plot(ps118_bathy, col=gray.colors(100, 0.9, 0.1))
plot(hotspot.rast, add=T, col = c("transparent","yellow2"), legend=F)

# Export as GeoTiff
writeRaster(hotspot.rast, "hotspots.tif", overwrite=T)

#########################

# Modelling
cluster.summary <- ps118_hydro[, .(
  hotspot_index = mean(hotspot_index, na.rm=T),
  thetao = mean(thetao, na.rm=T),
  CHL = mean(CHL, na.rm=T),
  uo = mean(uo, na.rm=T),
  vo = mean(vo, na.rm=T),
  bottomT = mean(bottomT, na.rm=T),
  siconc = mean(siconc, na.rm=T),
  mlotst = mean(mlotst, na.rm=T)), by=cluster]

gam.hotspot <- gam(
  hotspot_index ~ s(thetao) + s(uo) + s(vo) +
    s(bottomT) + s(siconc) + s(mlotst) + s(CHL),
  data = cluster.summary,
  method = "REML")
summary(gam.hotspot)

