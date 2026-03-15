
#################################################################
 ## MATCH BATHYMETRY and HYDROGRAPHY
#################################################################

# Load CMEMS hydrography -- 3D
oce3d <- tidync("cmems_mod_glo_phy_my_0.083deg_P1M.nc") %>%
  hyper_tibble(select_var = c("so","thetao","uo","vo")) %>%
  mutate(across(c(longitude, latitude, depth), as.numeric),
         lon = longitude, lat = latitude, date = time) %>%
  dplyr::select(lon, lat, depth, date, so, thetao, uo, vo) 

# Load CMEMS hydrography -- 2D
oce2d <- tidync("cmems_mod_glo_phy_my_0.083deg_P1M.nc") %>%
  activate("D2,D1,D0") %>%  
  hyper_tibble(select_var = c("bottomT","mlotst","siconc")) %>%
  rename(lon=longitude, lat=latitude, date=time) %>%
  mutate(across(c(lon, lat), as.numeric)) 

# Load CMEMS chlorophyll
chl <- tidync("cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M.nc") %>%
  activate("D2,D1,D0") %>%  
  hyper_tibble(select_var = c("CHL")) %>%
  rename(lon=longitude, lat=latitude, date=time) %>%
  mutate(across(c(lon, lat), as.numeric)) 

#############################

# Join hydrography (OCE)
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
knn.chl <- get.knnx(chl.mat, oce.mat, k=1)

# add CHL value + distance
oce.avg$CHL <- chl.avg$CHL[knn.chl$nn.index[,1]]
oce.avg$dist_to_chl <- knn.chl$nn.dist[,1]

# Prepare PS data, ensure right coord
ps118_sf <- hotspots %>% #ps118_combined
  filter(!is.na(x), !is.na(y), !is.na(depthReal)) %>%
  st_as_sf(coords = c("x","y"), crs=32722) %>%
  mutate(
    #depth_ps = abs(depth),     
    depth_ps = depthReal,  
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
# horiz_median 3229m -- vert_median 336m

# Define candidates
oce.cand <- oce.match[dist_vert <= 100 & dist_horiz <= 10000]
#oce.cand[, retained_count := .N, by = id_ps]

# Set order; find nearest candidate per id_ps
setorder(oce.cand, id_ps, dist_3D)
oce.filtered <- oce.cand[!duplicated(id_ps)]

# Print stats 
cat(sprintf("Matched IDs: %d / %d (%.1f%%).",
  uniqueN(oce.filtered$id_ps),
  uniqueN(oce.match$id_ps),
  100 * uniqueN(oce.filtered$id_ps) / uniqueN(oce.match$id_ps)))

# For merging: convert to DT 
setDT(hotspots)  
setDT(oce.avg)      

# Append OCE averages
oce.filtered <- oce.avg[
  oce.filtered, on =.(x_oce, y_oce), mult="first"]

# Rename columns correctly 
oce.filtered[, depth_oce := NULL]
setnames(oce.filtered, "i.depth_oce","depth_oce")

# Append Hotspot data
ps118_hydro <- hotspots[
  oce.filtered, on =.(x = x_ps, y = y_ps), mult="first"]


#################################################################
 ## HYDROGRAPHIC REGIMES 
#################################################################

# Subsample 
set.seed(123)
pca.hydro <- ps118_hydro[sample(.N, min(.N, 1e4))]

# Capping to handle outliers (99th percentile for total densities)
pca.hydro[, total_capped := pmin(total_density, quantile(total_density, 0.99, na.rm=T))]

# PCA on hydro variables
hydroVars <- c("thetao","uo","vo","bottomT","siconc","mlotst","CHL","so")
pcs <- prcomp(as.matrix(pca.hydro[, ..hydroVars]), scale.=T)
pca.hydro[, c("PC1","PC2") := .(pcs$x[,1], pcs$x[,2])]
scores <- pcs$x[, 1:2]

# K-Means: ONCE, so kmeanClust and pca.hydro$cluster use same assignment
set.seed(123)
kmeans <- kmeans(scores, centers=4, nstart=25)

# Extract outcome
kmeanClust <- data.frame(
  Cluster = 1:nrow(kmeans$centers),
  Size    = as.vector(table(kmeans$cluster)),
  kmeans$centers)
print(kmeanClust)

# Silhouette test for documentation
kmeanSil = data.frame(
  k = 2:8, 
  avg_sil = sapply(2:8, function(k) {
    km <- kmeans(scores, centers = k, nstart=25)
    mean(silhouette(km$cluster, dist(scores))[, 3])}))

# Assign clusters to main data
pca.hydro[, regime := kmeans$cluster]

# Correlations (hydroVars, slope, total density) with clusters
slopeCor <- cor(pca.hydro$ps118_slope, pca.hydro[, .(PC1, PC2)])
depthCor <- cor(pca.hydro$depthReal, pca.hydro[, .(PC1, PC2)])
totalCor <- cor(pca.hydro$total_density, pca.hydro[, .(PC1, PC2)])

# Determine and select key variables
# based on summed PC1+PC2 loadings
# vo/uo omitted: low impact + directional (hard to interpret)
cor(pca.hydro$total_capped, pca.hydro[, .(PC1, PC2)])
rotate <- as.data.frame(pcs$rotation[,1:2]) %>% 
  rownames_to_column("var") %>% 
  filter((abs(PC1) + abs(PC2)) > 0.52) %>% 
  column_to_rownames("var")

# Format arrows for best plotting
arrow <- rbind(
  data.table(x=0, y=0, xend=rotate[,1]*4.4, yend=rotate[,2]*4.4, Var=rownames(rotate)),
  data.table(x=0, y=0, xend=slopeCor[1]*4.4, yend=slopeCor[2]*4.4, Var="Slope"),
  data.table(x=0, y=0, xend=totalCor[1]*4.4, yend=totalCor[2]*4.4, Var="Density"),
  data.table(x=0, y=0, xend=depthCor[1]*4.4, yend=depthCor[2]*4.4, Var="Depth"))

# Format labels
label <- data.table(x=arrow$xend*1.2, y=arrow$yend*1.2, label_org=arrow$Var, label = c(
  "T","bottomT","MLD","chl","slope","density","depth"))

# Summarise hydroVars per cluster
hydro.summary <- pca.hydro[, .(
  density = median(total_capped),
  density_iqr = IQR(total_capped),
  depth = median(depthReal),
  depth_iqr = IQR(depthReal),
  bottomT = median(bottomT),
  bottomT_iqr = IQR(bottomT),
  potT = median(thetao),
  potT_iqr = IQR(thetao),
  slope = median(ps118_slope),
  slope_iqr = IQR(ps118_slope),
  siconc  = median(siconc),
  siconc_iqr = IQR(siconc),
  chl = median(CHL),
  chl_iqr = IQR(CHL),
  MLD  = median(mlotst),
  MLD_iqr  = IQR(mlotst),
  n = .N), by = regime][order(-density)][, regime_new := .I]

## PLOT
# export size 4 x 5.5
merge(pca.hydro, hydro.summary[, .(regime, regime_new, density)], by="regime") %>%
  distinct(regime_new, PC1, PC2, .keep_all=T) %>%
  slice_sample(n = min(group_size(.))) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(
    aes(shape=factor(regime_new), colour=factor(regime_new), size=factor(regime_new)),
    position = position_jitter(width=0.03, height=0.03),  alpha=0.6) +
  scale_colour_manual(
    values = c("1"="purple4","2"="orange3","3"="honeydew3","4"="skyblue"),
    name = "Hydrographic regime") +
  scale_shape_manual(
    values = c(16, 17, 15, 18), 
    name = "Hydrographic regime") +
  geom_segment(
    data = arrow,
    aes(x = x, y = y, xend = xend*1.2, yend = yend*1.2),
    arrow = arrow(length = unit(0.2, "cm")), color="gray20",
    linewidth = 0.5) +
  geom_label(
    data = label,
    aes(x=x*1.35, y=y*1.35, label=label),
    size=2.4,  fontface="bold", fill="white", color="black", 
    label.size=0.1, label.padding=unit(0.10, "lines"),
    nudge_x=0.02 * sign(label$x),
    nudge_y=0.02 * sign(label$y)) +
  scale_size_manual(
    name="Median density", 
    values=c("1"=4.2,"2"=3.3,"3"=2.6,"4"=2.7), 
    labels=paste0(hydro.summary$regime_new," (",round(hydro.summary$density,2),")")) +
  xlab(paste0("PC1 (", round(100*summary(pcs)$importance[2,1],1),"%)")) +
  ylab(paste0("PC2 (", round(100*summary(pcs)$importance[2,2],1),"%)")) +
  theme_classic()

#####################

## STATS
# Different densities, depth, slope between regimes?
hydro.stats <- pca.hydro %>%
  select(regime, total_density, ps118_slope, depthReal, bottomT, thetao, siconc, mlotst, CHL) %>%  
  pivot_longer(cols= -regime, names_to="variable", values_to="value") %>%
  group_by(variable) %>%
  group_modify(~ {
    eff <- kruskal_effsize(value ~ regime, data= .x)
    med <- .x %>%
      group_by(regime) %>%
      summarise(
        median=median(value, na.rm=T),
        IQR=IQR(value, na.rm=T),.groups="drop")
    med %>%
      mutate(effsize=eff$effsize, magnitude=eff$magnitude)}) %>%
  ungroup() %>%
  arrange(variable, regime)

# Contribution to cluster separation
dunn.hydro <- hydroVars[hydroVars %in% arrow$Var] %>%
  lapply(function(v) {
    pca.hydro %>%
      select(regime, all_of(v)) %>%
      mutate(value= .data[[v]]) %>%  
      dunn_test(value ~ regime, p.adjust.method="BH") %>% 
      left_join(
        pca.hydro %>%
          select(regime, all_of(v)) %>%
          mutate(value= .data[[v]]) %>%
          kruskal_effsize(value ~ regime) %>%
          select(effsize, magnitude),
        by=character()) %>%
      mutate(variable=v)}) %>%
  bind_rows() %>%
  select(variable, group1, group2, statistic, p, p.adj, p.adj.signif, effsize, magnitude) %>%
  group_by(variable) %>%
  summarize(
    mean_effsize=mean(effsize),
    magnitude=first(magnitude)) %>%
  mutate(
    PC1_loading=pcs$rotation[variable, 1],
    PC2_loading=pcs$rotation[variable, 2]) %>%
  mutate(
    abs_PC1=abs(PC1_loading),
    abs_PC2=abs(PC2_loading),
    PC_rank=rank(-abs_PC1),             # Higher absolute loading=higher rank
    effsize_rank=rank(-mean_effsize),   # Higher effect size=higher rank
    combined_rank=round((effsize_rank + PC_rank) / 2)) %>%  # rounded mean
  arrange(combined_rank)

# Print for SI table
hydro.summary %>%
  transmute(
    Regime = paste0("Regime ", regime_new),
    `Density` = paste0(round(density, 2), " (", round(density_iqr, 2), ")"),
    `Depth (m)` = paste0(round(depth, 0), " (", round(depth_iqr, 0), ")"),
    `Bottom Temp` = paste0(round(bottomT, 2), " (", round(bottomT_iqr, 2), ")"),
    `Pot. Temp` = paste0(round(potT, 2), " (", round(potT_iqr, 2), ")"),
    `Slope (°)` = paste0(round(slope, 0), " (", round(slope_iqr, 0), ")"),
    `Sea-ice (%)` = paste0(round(siconc*100, 0), " (", round(siconc_iqr, 0), ")"),
    `Chl-a` = paste0(round(chl, 2), " (", round(chl_iqr, 2), ")"),
    `MLD (m)` = paste0(round(MLD, 1), " (", round(MLD_iqr, 1), ")")) %>%
  pivot_longer(-Regime, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = Regime, values_from = Value) %>% 
  write.table("./results/hydroRegimes.txt", sep="\t", row.names=F, quote=F)


#################################################################
 ## HOTSPOT STATS & HYDROGRAPHY LINKS 
#################################################################

## DENSITY VS DEPTH, TEMP, SLOPE
hydro.hot <- bind_rows(
  ps118_hydro %>%
    filter(cluster_all_hot_global > 0) %>%
    group_by(cluster_id=cluster_all_hot_global) %>%
    summarise(
      cluster_type="global",
      n = n(),
      slope = median(ps118_slope, na.rm=T),
      tempB = median(bottomT, na.rm=T),
      temp  = median(thetao, na.rm=T),
      depth = median(depthReal, na.rm=T),
      density = median(total_density, na.rm=T), .groups="drop"),
  ps118_hydro %>%
    filter(cluster_all_hot_depth > 0) %>%
    group_by(cluster_id = cluster_all_hot_depth) %>%
    summarise(
      cluster_type="depth",
      n       = n(),
      slope   = median(ps118_slope, na.rm=T),
      tempB   = median(bottomT, na.rm=T),
      temp    = median(thetao, na.rm=T),
      depth   = median(depthReal, na.rm=T),
      density = median(total_density, na.rm=T), .groups="drop"))

# FIGURE X - export size
ggplot(hydro.hot, aes(x=slope, y=tempB)) +
  geom_smooth(
    method="loess", span=1, se=T, 
    fill="grey84", colour="darkred", linewidth=0.8) +
  geom_point(
    aes(size = density, colour = depth), alpha=0.9) +
  scale_colour_viridis_c(
    option="mako", direction= -1, end=0.98, name="Median depth") +
  scale_size_continuous(
    range = c(2,7), name = "Median density (ind/m²)") +
  labs(
    x="Median slope (°)", y="Median bottom temperature (°C)") +
  theme_classic() +
  theme(
    axis.title = element_text(face="bold"),
    legend.position = "right")

# Enriched in hydrographic regime #1?
merge(pca.hydro, hydro.summary[, .(regime, regime_new)], by="regime") %>%
  mutate(
    hot_both  = hot_global & hot_depth,
    regime_new = factor(regime_new)) %>%
  group_by(regime_new) %>%
  summarise(
    obs   = sum(hot_both, na.rm=T),
    total = n(),.groups = "drop") %>%
  mutate(
    prop_obs = obs / total,
    expected_prop = sum(obs) / sum(total),
    enrichment = prop_obs / expected_prop,
    p_val = purrr::map_dbl(seq_len(n()), ~{
      a <- obs[.x]
      b <- total[.x] - obs[.x]
      c <- sum(obs) - a
      d <- (sum(total) - total[.x]) - c
      suppressWarnings(chisq.test(matrix(c(a,b,c,d), nrow=2), correct=F)$p.value)}),
    p_adj = p.adjust(p_val, method="BH")) %>%
  arrange(desc(enrichment))

############################

## COMPLETE STATS
# Dunn: ENV differences between zones?
set.seed(42)
dunn.hot <- hotspots %>%
  select(x, y, total_density, hot_global, hot_depth, ps118_slope, depthReal) %>%
  left_join(hottest %>% select(x, y, Hotspot), by = c("x","y")) %>%
  left_join(ps118_hydro %>% select(x, y, bottomT, siconc, thetao, mlotst, CHL), by = c("x","y")) %>%
  mutate(zone = case_when(
    !is.na(Hotspot) ~"Hotspot",
    hot_global & hot_depth ~"High-density cluster",
    TRUE ~ "Non-Hotspot")) %>%
  group_by(zone) %>%
  arrange(stats::runif(n()), .by_group=T) %>%
  slice_head(n=5000) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(depthReal, total_density, ps118_slope, bottomT, siconc, thetao, CHL, mlotst),
    names_to="variable", values_to="value") %>%
  filter(!is.na(value)) %>%
  group_by(variable) %>%
  filter(n_distinct(zone) >1) %>%
  group_modify(~{
    dunn <- dunn_test(.x, value~zone, p.adjust.method="BH")
    eff  <- wilcox_effsize(.x, value~zone) %>%
      select(group1, group2, effsize, magnitude)
    left_join(dunn, eff, by = c("group1","group2"))}) %>%
  ungroup() %>%
  select(variable, group1, group2, n1, n2, p.adj, effsize, magnitude)

# Feature-level medians (only Dunn-significant) 
hot.summary <- {
  vars_keep <- dunn.hot %>%
    filter(p.adj < 0.05) %>%
    distinct(variable) %>%
    pull(variable) %>%
    setdiff(grep("cor|demo|glass|ophi|sea", ., value=TRUE))
  df0 <- hotspots %>%
    dplyr::select(
      x, y, xr, yr,
      total_density, ps118_slope, depthReal,
      hot_global, hot_depth) %>%
    left_join(
      hottest %>% dplyr::select(x, y, Hotspot, hotspot_id),
      by = c("x","y")) %>%
    left_join(
      ps118_hydro %>% dplyr::select(x, y, bottomT, siconc, mlotst, thetao),
      by = c("x","y")) %>%
    mutate(zone = case_when(
      !is.na(Hotspot) ~ "Hotspot",
      hot_global & hot_depth ~"High-density cluster",
      TRUE ~"Non-Hotspot"))
  df_units <- df0 %>%
    mutate(
      unit_id = case_when(
        zone == "Hotspot" ~ paste0("HS_", hotspot_id),
        zone == "High-density cluster" ~ paste0("HD_", xr, "_", yr),
        TRUE ~ paste0("NH_", xr, "_", yr))) %>%
    group_by(zone, unit_id) %>%
    summarise(
      n_cells_unit = n(),
      across(
        c(total_density, ps118_slope, depthReal, bottomT, mlotst, siconc, thetao),
        ~ median(.x, na.rm=T),
        .names = "unit_median_{.col}"),.groups="drop")
  df_zone <- df_units %>%
    group_by(zone) %>%
    summarise(
      n_units = n(),
      n_cells = sum(n_cells_unit),
      across(
        starts_with("unit_median_"),
        ~ median(.x, na.rm=T),
        .names = "median_{.col}"), .groups="drop") %>%
    mutate(
      km2 = n_cells *625 / 1e6,
      perc_area = km2 / sum(km2) *100)
  df_zone %>%
    pivot_longer(
      cols = starts_with("median_unit_median_"),
      names_to = "variable",
      values_to = "median") %>%
    mutate(variable = sub("^median_unit_median_", "", variable)) %>%
    filter(variable %in% vars_keep) %>%
    mutate(across(c(median, km2, perc_area), ~ round(.x, 3))) %>%
    arrange(variable, zone)}

## PLOT: ENV heatmap
# Export size 3x3.5
hot.summary %>%
  mutate(
    variable = case_when(
      variable == "ps118_slope" ~ "Slope (°)",
      variable == "bottomT" ~ "Bottom temperature (°C)",
      variable == "thetao" ~ "Potential temperature (°C)",
      variable == "mlotst" ~ "Mixed-layer depth (m)",
      variable == "depthReal" ~ "Depth (m)",
      variable == "siconc" ~ "Sea-ice (%)",
      variable == "total_density" ~ "Total density"),
    zone = factor(zone, levels = c("Hotspot","High-density cluster","Non-Hotspot")),
    variable = factor(variable, levels=rev(c(
      "Total density","Slope (°)","Potential temperature (°C)",
      "Sea-ice (%)","Depth (m)","Mixed-layer depth (m)","Bottom temperature (°C)"))),
    label_text = case_when(
      variable == "Sea-ice (%)" ~ sprintf("%.0f", median *100),
      variable %in% c("Slope (°)","Total density","Depth (m)") ~ sprintf("%.0f", median),
      TRUE ~ sprintf("%.2f", median))) %>%
  group_by(variable) %>%
  mutate(
    median_scaled = (median - min(median, na.rm=T)) /
    (max(median, na.rm=T) - min(median, na.rm=T) + 1e-9)) %>%
  ungroup() %>%
  ggplot(aes(x = zone, y = variable, fill = median_scaled)) +
  geom_tile(color="white", linewidth=0.4) +
  geom_text(aes(label = label_text), size=3.1) +
  scale_fill_gradient(low="aliceblue", high="cadetblue3") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position="none")


#################################################################
 ## TEMPERATURE DECOUPLING AROUND HOTSPOTS
#################################################################

# Calculate delta
temp.delta <- ps118_hydro %>%
  mutate(delta_temp = thetao - bottomT)

# Plot distribution
set.seed(44)
temp.delta %>%
  group_by(hot_global) %>%
  sample_n(min(n(), 10000)) %>%
  ungroup() %>%
  ggplot(aes(x = delta_temp, fill = hot_global)) +
  geom_density(alpha = 0.5, color="white") + 
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_fill_manual(
    values = c("grey","purple3"), 
    labels = c("Non-Hotspots","Hotspots")) +
  labs(x = expression(theta[o] - bottomT), y="Density", fill="Group") +
  theme_minimal()

# GAM
gam.temp <- gam(
  delta_temp ~ s(depthReal) + s(ps118_slope) + hot_global, 
  data = temp.delta %>% sample_n(100000)) 
summary(gam.temp)

# Averages per zone
hotspots %>%
  left_join(hottest %>% select(x, y, Hotspot), by=c("x","y")) %>%
  left_join(ps118_hydro %>% select(x, y, thetao, bottomT), by = c("x","y")) %>%
  mutate(delta_temp = thetao - bottomT) %>%
  mutate(zone = case_when(
    !is.na(Hotspot) ~ "Hotspot",
    hot_global & hot_depth ~ "High-density cluster",
    TRUE ~ "Non-hotspot")) %>%
  group_by(zone) %>%
  summarise(
    mean_delta = mean(delta_temp, na.rm=T),
    sd_delta   = sd(delta_temp, na.rm=T),
    mean_slope = mean(ps118_slope, na.rm=T),
    n = n()) %>%
  arrange(desc(abs(mean_delta)))
