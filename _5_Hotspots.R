
#################################################################
 ## BENTHIC HOTSPOTS
#################################################################

# Define "hot" grid cells (90% density quantiles)
# Both global and depth-binned (500m bins)
# Set grid coordinates for clustering
hotspots <- bioPredict %>%
  dplyr::select(x, y, depthReal, ps118_slope, total_density) %>%
  mutate(
    hot_global = total_density >= quantile(total_density, 0.9, na.rm=T),
    depth_cat = cut(
      depthReal,
      breaks = seq(0, max(depthReal, na.rm=T), by=500),
      include.lowest=T)) %>%
  group_by(depth_cat) %>%
  mutate(
    hot_depth = total_density >= quantile(total_density, 0.9, na.rm=T)) %>%
  ungroup() %>%
  mutate(
    xr = round(x / 500) * 500,
    yr = round(y / 500) * 500)

## Hotspot assignment 
for(type in c("hot_depth", "hot_global")) {
  grid <- hotspots %>%
    group_by(xr, yr) %>%
    summarise(
      frac_hot     = mean(.data[[type]], na.rm=T),
      mean_density = mean(total_density, na.rm=T),
      mean_depth   = mean(depthReal, na.rm=T),
      mean_slope   = mean(ps118_slope, na.rm=T),
      .groups = "drop") %>%
    filter(frac_hot >0.5) %>%  
    mutate(
      density_z    = as.numeric(scale(mean_density)),
      hot_strength = frac_hot * density_z)
  
  cl_coarse <- dbscan(as.matrix(grid[, c("xr","yr")]), eps=1500, minPts=4)
  grid$cluster_coarse <- cl_coarse$cluster
  
  grid <- grid %>%
    split(.$cluster_coarse) %>%
    lapply(function(d) {
      if(nrow(d) < 4){
        d$cluster_fine <- NA_integer_
      } else {
        d$cluster_fine <- dbscan(as.matrix(d[, c("xr","yr")]), eps=500, minPts=3)$cluster
      }
      d
    }) %>%
    bind_rows() %>%
    filter(cluster_coarse > 0)
 
  # all-cell assignment
  knn_all <- FNN::get.knnx(
    data  = as.matrix(grid[, c("xr","yr")]),
    query = as.matrix(hotspots[, c("xr","yr")]),
    k = 1)
  hotspots[[paste0("cluster_all_", type)]]    <- grid$cluster_coarse[knn_all$nn.index[,1]]
  hotspots[[paste0("subcluster_all_", type)]] <- grid$cluster_fine[knn_all$nn.index[,1]]
  hotspots[[paste0("strength_all_", type)]]   <- grid$hot_strength[knn_all$nn.index[,1]]

  # hot-only assignment
  hotspots[[paste0("cluster_", type)]]    <- 0L
  hotspots[[paste0("subcluster_", type)]] <- NA_integer_
  hotspots[[paste0("strength_", type)]]   <- NA_real_
 idx <- which(hotspots[[type]] %in% TRUE)
  
  if(length(idx) > 0){
    knn_hot <- FNN::get.knnx(
      data  = as.matrix(grid[, c("xr","yr")]),
      query = as.matrix(hotspots[idx, c("xr","yr")]),
      k = 1)
    hotspots[[paste0("cluster_", type)]][idx]    <- grid$cluster_coarse[knn_hot$nn.index[,1]]
    hotspots[[paste0("subcluster_", type)]][idx] <- grid$cluster_fine[knn_hot$nn.index[,1]]
    hotspots[[paste0("strength_", type)]][idx]   <- grid$hot_strength[knn_hot$nn.index[,1]]
  }
}

# Determine quantiles
quantile(hotspots$strength_hot_depth, probs=c(0.5, 0.75, 0.9, 0.95), na.rm=T)
quantile(hotspots$strength_hot_global, probs=c(0.5, 0.75, 0.9, 0.95), na.rm=T)

## FILTER accordingly; min 6000 cells
# Sort / rename by strength, slope, depth, location
hottest <- hotspots %>%
  filter(hot_global==TRUE & hot_depth==TRUE) %>%
  filter(strength_hot_global >0.84 & strength_hot_depth >0.88) %>% # global 0.63 & depth > 0.7
  mutate(Hotspot = dbscan(as.matrix(cbind(x, y)), eps=5000, minPts=3)$cluster) %>%
  filter(Hotspot > 0) %>%
  group_by(Hotspot) %>%
  mutate(
    n_cells  = n(),
    strength = mean(strength_hot_global, na.rm=T)) %>%
  filter(n_cells >6000) %>%
  ungroup() %>%
  left_join(
    {df <- .
      cent <- df %>%
        group_by(Hotspot) %>%
        summarise(
          cx = mean(x, na.rm =T),
          cy = mean(y, na.rm =T),
          .groups = "drop")
      east_ids <- cent %>%
        arrange(desc(cx), Hotspot) %>%
        slice_head(n = 4) %>%
        pull(Hotspot)
     map_east <- cent %>%
        filter(Hotspot %in% east_ids) %>%
        arrange(desc(cy), Hotspot) %>%   #north -> south
        mutate(hotspot_id = row_number())
     map_rest <- cent %>%
        filter(!Hotspot %in% east_ids) %>%
        arrange(desc(cx), Hotspot) %>%  # east -> west for remaining
        mutate(hotspot_id = row_number() + nrow(map_east))
  bind_rows(map_east, map_rest) %>%
        select(Hotspot, hotspot_id)
},  by="Hotspot") 
  
####################################

# Taxon enrichment + hotspot contribution
hottest.enrch <- hotspots %>%
  select(x, y, hot_global, hot_depth) %>%
  left_join(hottest %>% select(x, y, Hotspot), by=c("x","y")) %>%
  left_join(bioPredict %>% select(x, y, matches("cor|demo|glass|ophi|sea")), by=c("x","y")) %>%
  mutate(zone = case_when(
    !is.na(Hotspot) ~"Hotspot",
    hot_global & hot_depth ~"High-density cluster",
    TRUE ~ "Non-Hotspot")) %>%
  group_by(zone) %>%
  summarise(across(matches("cor|demo|glass|ophi|sea"), ~ median(.x, na.rm=T)), .groups="drop") %>%
  pivot_longer(matches("cor|demo|glass|ophi|sea"), names_to="variable", values_to="median") %>%
  pivot_wider(names_from = zone, values_from = median) %>%
  mutate(
    EF_hotspot = Hotspot / (`Non-Hotspot` + 1e-6),
    contrib_hotspot = Hotspot / (sum(Hotspot, na.rm=T) + 1e-12)) %>%
  arrange(desc(EF_hotspot)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Canyon metrics 
hottest.summary <- hottest %>%
  group_by(hotspot_id) %>%
  summarise(
    n_cells = n(),
    area_km2 = n_cells * 625 / 1e6,
    mean_depth = mean(depthReal, na.rm=T),
    mean_slope = mean(ps118_slope, na.rm=T),
    mean_strength = mean(strength_hot_global, na.rm=T),
    median_density = median(total_density, na.rm=T),
    min_x = min(x, na.rm=T),
    max_x = max(x, na.rm=T),
    min_y = min(y, na.rm=T),
    max_y = max(y, na.rm=T),
    width_km  = (max_x - min_x) / 1000,
    height_km = (max_y - min_y) / 1000,
    ratio = width_km / (height_km + 1e-9), .groups="drop") %>%
  arrange(ratio) %>%
  select(-min_x, -min_y, -max_x, -max_y)

# Spatial extents 
hotspots %>%
  left_join(hottest %>% select(x, y, Hotspot), by = c("x","y")) %>%
  mutate(zone = case_when(
    !is.na(Hotspot) ~"Hotspot",
    hot_global & hot_depth ~"High-density cluster",
    TRUE ~"Non-Hotspot")) %>%
  group_by(zone) %>%
  summarise(
    area_km2 = (n() * 625) / 1e6,
    area_fraction = (n() / nrow(hotspots)) *100,
    total_km2 = (nrow(hotspots) *625) / 1e6) %>%
  arrange(desc(area_km2))
# Hotspots: 33 sqkm, 0.5% of study area 
# High-density clusters: 7%
# Non-hotspots 92%
# Total area 7400 km2

####################################

## EXPORT / EXPLORATORY PLOTS
# Load base raster
ps118_bathy <- rast("PS118_bathy_25m.tif")

# Create rasters for hotspots
hottest.sf <- st_as_sf(
  hottest,
  coords = c("x","y"),
  crs = crs(ps118_bathy))

# Define ecoregion centroids
hottest.centroids <- hottest.sf %>%
  group_by(hotspot_id) %>%
  summarise(
    mean_depth = mean(depthReal, na.rm=T),
    mean_density = mean(total_density, na.rm=T),
    n_points = n(),
    geometry = st_centroid(st_union(geometry)))

# Hottest coordinates
hottest.coords <- st_coordinates(hottest.centroids)

# Export hotspot points for ArcGis
st_write(hottest.sf, "~/Rdata/HotspotPoints.gpkg", delete_dsn=T)

# Raster for high-density clusters
clust.rast <- rast(ps118_bathy)
clust.rast[] <- NA
xy <- hotspots %>%
  filter(hot_global==T & hot_depth==T) %>%
  select(x, y) %>%
  as.matrix()
cells <- cellFromXY(clust.rast, xy)
cells <- unique(cells[!is.na(cells)])   
clust.rast[cells] <- 1

# Raster for hotspots
hot.rast <- rast(ps118_bathy)
xy <- as.matrix(hottest[, c("x","y")])
cells <- cellFromXY(hot.rast, xy)
hot.rast[] <- NA
hot.rast[cells] <- 1

# Plot points + centroids
plot(ps118_bathy, col = gray.colors(100, 0.9, 0.1), legend=T)
plot(hot.rast, add=T, cex=5, col = c("transparent","deeppink3"), legend=F)
plot(clust.rast, add=T, col = c("transparent","orange"), legend=F)
text(x = hottest.coords[,1], 
     y = hottest.coords[,2], 
     labels = hottest.centroids$hotspot_id, 
     col="white", cex=0.8, font = 2)   
