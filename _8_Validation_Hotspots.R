
#################################################################
 ## HOTSPOT VERIFICATION
#################################################################

## Do predicted hotspot cells indeed host greater densities?
## FIRST: annotate random new images predicted as hot v not
# Load complete image info
imagesAll <- read.table(
  "./envData/imagesComplete.txt", header=T, sep="\t", check.names=F) %>%
  select(-type)

# Load 32 newly annotated images 
veri.raw <- read.table(
  "./envData/hotspotVeri.txt", header=T, sep="\t", check.names=F) 

# Extract unique metadata (area and hotspot/non-hotspot type)
# Raw counts per image+organism 
# Sum demosponges + glass sponges
# Normalize (density per m2)
veri.counts <- veri.raw %>%
  filter(labels != "laser point") %>%
  mutate(area = image_m2 * rectify) %>%
  group_by(image_filename, labels, type) %>%
  summarise(raw_count = n(), area = first(area), .groups="drop") %>%
  pivot_wider(names_from = labels, values_from = raw_count, values_fill = 0) %>%
  mutate(
    obs_corals = corals / area,
    obs_sponges = (demosponges + `glass sponges`) / area,
    obs_total = obs_corals + obs_sponges) %>%
  select(image_filename, type, area, obs_corals, obs_sponges, obs_total) %>%
  left_join(imagesAll, by="image_filename")

# Wilcoxon test per taxon
veri.counts %>%
  dplyr::select(-area, -date, -datetime, -id, -image_URL, -image_txtfile) %>%
  pivot_longer(cols = -c(image_filename, type), names_to="taxon", values_to ="density") %>%
  mutate(type = factor(type, levels = c("Non-Hotspot","Hotspot"))) %>%
  group_by(taxon) %>%
  filter(length(unique(type)) == 2) %>% 
  pairwise_wilcox_test(density ~ type, p.adjust.method="BH") %>%
  ungroup()

# Plot
veri.counts %>%
  dplyr::select(-area, -date, -datetime, -id, -image_URL, -image_txtfile) %>%
  pivot_longer(cols = -c(image_filename, type), names_to="taxon", values_to ="density") %>%
  filter(taxon %in% c("obs_corals","obs_sponges")) %>%
  ggplot(aes(x = taxon, y = density, fill = type)) +
  geom_boxplot(outliers=F, alpha=0.8, width=0.5) +
  scale_fill_manual(values = hotspotCol) + 
  ylab(expression(bold("Density (individuals / ") ~ bold(m^2) ~ bold(")"))) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth=0.22, color="gray84"),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face="bold", size = 11),
    legend.position = "right", 
    axis.title.x = element_blank(),
    axis.text.x = element_text(face="bold", size = 10)) +
  ggpubr::stat_compare_means(
    method = "wilcox.test", label = "p.signif", 
    label.x = 1.5, hide.ns = FALSE)

#############################

## SECOND: Do observed counts match predicted densities from upscaling?
# Match image coordinates to nearest prediction cell
veri.xy <- veri.counts %>%
  terra::vect(geom = c("lon","lat"), crs="EPSG:4326") %>%
  terra::project(crs) %>%
  as.data.frame(geom="XY")

# Set nearest-neighbor
knn.hot <- FNN::get.knnx(
  data  = as.matrix(bioPredict[, c("x","y")]),
  query = as.matrix(veri.xy[, c("x","y")]),
  k = 1)

# Match data
veri.pred <- veri.xy %>% 
  bind_cols(bioPredict[knn.hot$nn.index[, 1], c("corals","demosponges","glassSponges","total_density")] %>%
  rename(pred_corals = corals, pred_total = total_density)) %>%
  mutate(pred_sponges = demosponges + glassSponges, nn_dist_m = knn.hot$nn.dist[, 1])

# Compact summary table
veri.stats <- bind_rows(
  veri.pred %>% group_by(type) %>% summarise(
    metric = "Observed median density",
    n = n(),
    total = median(obs_total, na.rm=T),
    corals = median(obs_corals, na.rm=T),
    sponges = median(obs_sponges, na.rm=T),
    pred_total = median(pred_total, na.rm=T),
    .groups="drop"),
  tibble(
    type = "Hotspot vs Non-Hotspot",
    metric = "Effect size",
    n = nrow(veri.pred),
    total = veri.pred %>% wilcox_effsize(obs_total ~ type) %>% pull(effsize),
    corals = veri.pred %>% wilcox_effsize(obs_corals ~ type) %>% pull(effsize),
    sponges = veri.pred %>% wilcox_effsize(obs_sponges ~ type) %>% pull(effsize),
    pred_total = NA_real_),
  tibble(
    type = "Predicted vs observed",
    metric = "Spearman rho",
    n = nrow(veri.pred),
    total = cor(veri.pred$pred_total, veri.pred$obs_total, method="spearman"),
    corals = cor(veri.pred$pred_corals, veri.pred$obs_corals, method="spearman"),
    sponges = cor(veri.pred$pred_sponges, veri.pred$obs_sponges, method="spearman"),
    pred_total = NA_real_),
  tibble(
      type = "Spatial matching",
      metric = "Nearest prediction cell distance",
      n = sum(!is.na(veri.pred$nn_dist_m)), # Use count of non-NA values
      total = median(veri.pred$nn_dist_m, na.rm = TRUE),
      corals = min(veri.pred$nn_dist_m, na.rm = TRUE),
      sponges = max(veri.pred$nn_dist_m, na.rm = TRUE),
      pred_total = NA_real_)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Plot relation for coral+sponges
veri.pred %>%
  pivot_longer(
    cols = c(pred_sponges, obs_sponges, pred_corals, obs_corals),
    names_to = c(".value", "taxon"),
    names_sep = "_") %>%
  ggplot(aes(x = pred, y = obs, fill = type)) +
  geom_point(shape = 21, size = 3, color = "black", alpha = 0.85) +
  facet_wrap(~ taxon, scales = "free") +
  scale_fill_manual(values=hotspotCol) +
  theme_classic() +
  labs(
    x = expression("Predicted density (individuals " * m^-2 * ")"),
    y = expression("Observed density (individuals " * m^-2 * ")"),
    fill = NULL)


#################################################################
 ## HOTSPOT ROBUSTNESS vs THRESHOLDS
#################################################################

# Function for threshold comparison
hotCheck <- function(
    q = 0.90, strengthGlobal = 0.90, strengthDepth = 0.90,
    cells = 6000, eps = 5000, minPts_hot = 3) {
  min_cells <- cells
  eps_hot <- eps
  hotspots_tmp <- bioPredict %>%
    dplyr::select(x, y, depthReal, ps118_slope, total_density) %>%
    mutate(
      hot_global = total_density >= quantile(total_density, q, na.rm=T),
      depth_cat = cut(
        depthReal,
        breaks = seq(0, max(depthReal, na.rm=T), by = 500),
        include.lowest=T)) %>%
    group_by(depth_cat) %>%
    mutate(
      hot_depth = total_density >= quantile(total_density, q, na.rm=T)) %>%
    ungroup() %>%
    mutate(
      xr = round(x / 500) * 500,
      yr = round(y / 500) * 500)
  for(type in c("hot_depth", "hot_global")) {
    grid <- hotspots_tmp %>%
      group_by(xr, yr) %>%
      summarise(
        frac_hot     = mean(.data[[type]], na.rm=T),
        mean_density = mean(total_density, na.rm=T),
        mean_depth   = mean(depthReal, na.rm=T),
        mean_slope   = mean(ps118_slope, na.rm=T),
        .groups = "drop") %>%
      filter(frac_hot > 0.5) %>%
      mutate(
        density_z    = as.numeric(scale(mean_density)),
        hot_strength = frac_hot * density_z)
    if(nrow(grid) == 0) next
    cl_coarse <- dbscan::dbscan(as.matrix(grid[, c("xr", "yr")]), eps = 1500, minPts = 4)
    grid$cluster_coarse <- cl_coarse$cluster
    grid <- grid %>%
      split(.$cluster_coarse) %>%
      lapply(function(d) {
        if(nrow(d) < 4 || unique(d$cluster_coarse) == 0) {
          d$cluster_fine <- NA_integer_} else {
            d$cluster_fine <- dbscan::dbscan(as.matrix(d[, c("xr", "yr")]), eps = 500, minPts = 3)$cluster} 
        d}) %>%
      bind_rows() %>%
      filter(cluster_coarse > 0)
    hotspots_tmp[[paste0("strength_", type)]] <- NA_real_
    idx <- which(hotspots_tmp[[type]] %in% TRUE)
    if(length(idx) > 0 && nrow(grid) > 0) {
      knn_hot <- FNN::get.knnx(
        data  = as.matrix(grid[, c("xr", "yr")]),
        query = as.matrix(hotspots_tmp[idx, c("xr", "yr")]), k = 1)
      hotspots_tmp[[paste0("strength_", type)]][idx] <-
        grid$hot_strength[knn_hot$nn.index[, 1]]}}
  g_thr <- quantile(hotspots_tmp$strength_hot_global, strengthGlobal, na.rm=T)
  d_thr <- quantile(hotspots_tmp$strength_hot_depth,  strengthDepth,  na.rm=T)
  candidate <- hotspots_tmp %>%
    filter(hot_global==T & hot_depth==T) %>%
    filter(strength_hot_global > g_thr, strength_hot_depth  > d_thr)
  if(nrow(candidate) == 0) {
    stats <- tibble(
      q = q,
      strengthGlobal = strengthGlobal, strengthDepth = strengthDepth,
      global_strength_thr = g_thr, depth_strength_thr = d_thr,
      cells = cells, eps = eps, n_hotspots = 0, n_cells = 0,
      area_km2 = 0, area_percent = 0,
      median_density = NA_real_, mean_density = NA_real_)
    return(list(stats = stats, cells = candidate))}
  hottest_tmp <- candidate %>% mutate(
    Hotspot = dbscan::dbscan(as.matrix(cbind(x, y)),
                             eps = eps_hot, minPts = minPts_hot)$cluster) %>%
    filter(Hotspot > 0) %>%
    group_by(Hotspot) %>%
    mutate(
      n_cells = n(),
      strength = mean(strength_hot_global, na.rm=T)) %>%
    filter(n_cells > min_cells) %>%
    ungroup()
  stats <- tibble(
    q = q, strengthGlobal = strengthGlobal, strengthDepth = strengthDepth,
    global_strength_thr = g_thr, depth_strength_thr = d_thr,  cells = cells,
    eps = eps, n_hotspots = n_distinct(hottest_tmp$Hotspot),
    n_cells = nrow(hottest_tmp), area_km2 = nrow(hottest_tmp) * 625 / 1e6,
    area_percent = nrow(hottest_tmp) / nrow(hotspots_tmp) * 100,
    median_density = median(hottest_tmp$total_density, na.rm=T),
    mean_density = mean(hottest_tmp$total_density, na.rm=T))
  return(list(stats = stats, cells = hottest_tmp))}

# Define scenarios to test
hotcheck.scenarios <- tibble::tribble(
  ~scenario, ~q, ~strengthGlobal, ~strengthDepth, ~cells, ~eps,
  "Original: density Q90 + strength Q90 + eps 5 km + min 6000 cells", 0.90, 0.90, 0.90, 6000, 5000,
  "Density Q90 + strength Q90 + eps 3 km + min 6000 cells", 0.90, 0.90, 0.90, 6000, 3000,
  "Density Q90 + strength Q95 + eps 5 km + min 6000 cells", 0.90, 0.95, 0.95, 6000, 5000,
  "Density Q95 + strength Q90 + eps 5 km + min 6000 cells", 0.95, 0.90, 0.90, 6000, 5000,
  "Density Q99 + strength Q50 + eps 5 km + min 2000 cells", 0.99, 0.50, 0.50, 2000, 5000)  %>%
  mutate(scenario_label = str_replace_all(scenario, " \\+ ", "\n"))

# Calculate
hotcheck.runs <- hotcheck.scenarios %>%
  mutate(res = pmap(
    list(q, strengthGlobal, strengthDepth, cells, eps),
    ~ hotCheck(q = ..1, strengthGlobal = ..2, strengthDepth = ..3, cells = ..4, eps = ..5)))

# Plot numerical summary
hotcheck.runs %>%
  mutate(stats = map(res,"stats")) %>%
  select(scenario, stats) %>%
  tidyr::unnest(stats) %>%
  mutate(
    plot_area = as.numeric(area_km2),
    plot_density = ifelse(is.na(mean_density), 0, mean_density),
    is_hotspot = n_hotspots > 0) %>%
  ggplot(aes(
    x = plot_area, y = reorder(str_wrap(scenario, width=40), plot_area))) +
  geom_point(
    data = . %>% filter(!is_hotspot), shape=4, size=3, stroke=2) +
  geom_point(
    data = . %>% filter(is_hotspot), aes(fill = plot_density),
    shape=21, size=12, color="black", alpha=0.9) +
  geom_text(
    data = . %>% filter(is_hotspot), aes(label = round(area_percent, 2)), 
    color="white", fontface="bold", size=3.5) +
  geom_segment(
    data = . %>% filter(is_hotspot), aes(x=0, xend = plot_area - 1.7, 
    yend = reorder(str_wrap(scenario, width=40), plot_area))) +
  scale_fill_gradient(
    low="deepskyblue4", high="gold2", 
    limits = c(27.5, 29.5), breaks = c(27.5, 28.5, 29.5),
    oob = scales::squish) +
  theme_classic() +
  scale_x_continuous(expand=c(0.07,1)) +
  labs(
    x = expression("Hotspot area (" * km^2 * ")"), 
    y=NULL, fill="Mean Density") +
  theme(
    axis.text.y = element_text(lineheight = 1.1), 
    axis.ticks.y = element_blank())

# Summary table
hotcheck.runs %>%
  mutate(stats = map(res, "stats")) %>%
  select(scenario, stats) %>%
  tidyr::unnest(stats) %>%
  mutate(across(c(
    global_strength_thr, depth_strength_thr, area_km2,
    area_percent, median_density, mean_density), ~ round(.x, 3)))

# Prepare map plotting
hotcheck.plot <- hotcheck.runs %>%
  mutate(cells_df = map(res, "cells")) %>%
  select(scenario_label, cells_df) %>%
  tidyr::unnest(cells_df) %>%
  filter(n_cells > 0) %>%
  mutate(scenario_label = factor(scenario_label,levels = hotcheck.scenarios$scenario_label))

# Create base-plot
ps118_bathy <- rast("./envData/PS118_bathy_25m.tif")

# Plot map 
ggplot() +
  geom_spatraster(data = ps118_bathy) +
  scale_fill_gradient(low="white", high="black", name="Depth") +
  geom_point(data = hotcheck.plot, aes(x=x, y=y), color="red", size=0.5) +
  facet_wrap(~scenario_label) +
  theme_void()

