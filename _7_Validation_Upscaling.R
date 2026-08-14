
#################################################################
 ## UPSCALING UNCERTAINTY
#################################################################

## What is the prediction uncertainty in biogeographic upscaling? 
# Subsample prediction grid
set.seed(1)
uncert.grid <- bioPredict %>%
  inner_join(hot.lookup, by = c("x","y")) %>%
  slice_sample(n = min(200000, nrow(.))) %>%
  mutate(slope = slope_scaled, depth = depth_scaled) %>%
  dplyr::select(x, y, slope, depth, zone)

bootPred <- function(count_data, model_formula, grid, B=50) {
  preds <- replicate(B, {
    dat_b <- count_data[sample(seq_len(nrow(count_data)), replace=T), ]
    mod_b <- try(
      gam(model_formula, data = dat_b, family = nb()),
      silent=T)
    if(inherits(mod_b, "try-error")) {
      rep(NA_real_, nrow(grid))
    } else {
      predict(mod_b, newdata = grid, type="response")}})
  preds}

# Bootstrapping per taxon
set.seed(2)
uncert.taxa <- list(
  corals = bootPred(counts.coral, Corals ~ slope + depth, uncert.grid, B=50),
  demosponges = bootPred(counts.demo, Demosponges ~ slope + depth, uncert.grid, B=50),
  glassSponges = bootPred(counts.glass, GlassSponges ~ slope + depth, uncert.grid, B=50),
  ophiuroids = bootPred(counts.ophi, Ophiuroids ~ slope + depth, uncert.grid, B=50),
  seaPens = bootPred(counts.pens, SeaPens ~ slope + depth, uncert.grid, B=50))  %>%
  reduce(`+`)

# Summarise uncertainties
uncert.map <- uncert.grid %>%
  dplyr::select(x, y, zone) %>%
  bind_cols(tibble(
    median = apply(uncert.taxa, 1, median, na.rm=T),
    low95  = apply(uncert.taxa, 1, quantile, probs = 0.025, na.rm=T),
    high95 = apply(uncert.taxa, 1, quantile, probs = 0.975, na.rm=T),
    range95 = high95 - low95,
    uncertainty = range95 / (median + 1e-9))) 

# Overall summary
# Uncertainty map plotted with terrain-matching below
cat(sprintf(
  "Median: %s billion\n95%% CI: %s – %s billion\nRel. Range: %s%% / +%s%%",
  round(median(bs <- colSums(uncert.taxa, na.rm=T) * 625 * nrow(bioPredict) / nrow(uncert.grid) / 1e9), 1),
  round(quantile(bs, 0.025), 1),
  round(quantile(bs, 0.975), 1),
  round(100 * (quantile(bs, 0.025) / median(bs) - 1), 0),
  round(100 * (quantile(bs, 0.975) / median(bs) - 1), 0)))

# Comparable uncertainties in hotspots vs non-hotspots?
uncert.map %>%
  filter(is.finite(uncertainty)) %>%
  wilcox_effsize(
    uncertainty ~ zone,
    comparisons = list(
      c("Hotspot", "High-density cluster"),
      c("Hotspot", "Non-Hotspot"),
      c("High-density cluster","Non-Hotspot")))


#################################################################
 ## TERRAIN REPRESENTATION 
#################################################################

## Do our 196 images represent the environmental space across the flank?
# Load whole-grid bathymetry (slope+depth already loaded)
ps118_vrm   <- terra::rast("./envData/PS118_VRM_25m.tif")
ps118_east  <- terra::rast("./envData/PS118_eastness_25m.tif")
ps118_north <- terra::rast("./envData/PS118_northness_25m.tif")
set.seed(1)

# Terrain attributes for 196 images
env.img <- meta %>%
  tibble::rownames_to_column("image") %>%
  transmute(
    image, dive, depth = abs(depth),  slope = slope_25m,
    vrm = VRM_25m, east = aspect_E_25m, north = aspect_N_25m) %>%
  drop_na(image, dive, depth, slope, vrm, east, north)

# Regional terrain sample
env.grid <- ps118_combined %>%
  inner_join(hot.lookup, by = c("x","y")) %>%
  dplyr::select(x, y, zone, depth, slope = ps118_slope) %>%
  mutate(
    depth = abs(depth),
    vrm = terra::extract(ps118_vrm, terra::vect(., geom = c("x","y"), crs=crs))[, 2],
    east = terra::extract(ps118_east, terra::vect(., geom = c("x","y"), crs=crs))[, 2],
    north = terra::extract(ps118_north, terra::vect(., geom = c("x","y"), crs=crs))[, 2]) %>%
  drop_na(depth, slope, vrm, east, north) %>%
  slice_sample(n = min(500000, nrow(.)))

# Image-derived tolerances
# Depth remains fixed at ±100 m.
# All other terrain variables use ±1 SD across OFOBS images.
terrain.tol <- env.img %>%
  summarise(
    slope = sd(slope, na.rm=T),
    vrm = sd(vrm, na.rm=T),
    east = sd(east, na.rm=T),
    north = sd(north, na.rm=T)) %>%
  unlist()

# Function for cell-image matching
terrainMatch <- function(grid, img, tolerances) {
  n_images_per_cell <- integer(nrow(grid))
  n_cells_per_image <- integer(nrow(img))
  dives <- sort(unique(img$dive))
  dive_support <- matrix(
    FALSE, nrow = nrow(grid),
    ncol = length(dives), dimnames = list(NULL, dives))
  for(i in seq_len(nrow(img))) {
    inside <- rep(TRUE, nrow(grid))
  for(v in names(tolerances)) {
    inside <- inside & abs(grid[[v]] - img[[v]][i]) <= tolerances[[v]]}
    inside[is.na(inside)] <- FALSE
    n_images_per_cell <- n_images_per_cell + as.integer(inside)
    n_cells_per_image[i] <- sum(inside)
    dive_support[inside, img$dive[i]] <- TRUE}
  list(
    cell_support = n_images_per_cell,
    image_coverage = tibble(
      image = img$image,
      dive = img$dive,
      matched_grid_cells = n_cells_per_image,
      matched_grid_percent = 100 * n_cells_per_image / nrow(grid)),
   dive_support = as_tibble(dive_support) %>%
      mutate(n_matching_dives = rowSums(dive_support)))}

# Prepare depth + slope
match.sd <- terrainMatch(
  grid = env.grid,
  img = env.img,
  tolerances = c(
    depth = 100,
    slope = unname(terrain.tol["slope"])))

# Prepare depth + slope + VRM + eastness + northness
match.5var <- terrainMatch(
  grid = env.grid,
  img = env.img,
  tolerances = c(
    depth = 100,
    slope = unname(terrain.tol["slope"]),
    vrm = unname(terrain.tol["vrm"]),
    east = unname(terrain.tol["east"]),
    north = unname(terrain.tol["north"])))

# Final cell-level object
env.grid <- env.grid %>%
  mutate(
    n_match_sd = match.sd$cell_support,
    n_match_5var = match.5var$cell_support,
    match_sd = n_match_sd > 0,
    match_5var = n_match_5var > 0) %>%
  bind_cols(match.5var$dive_support %>% rename(
    dive_39_1 = `39-1`, dive_69_1 = `69-1`, dive_81_1 = `81-1`))

# Terrain-match summary (plotted below with dive info)
bind_rows(
  tibble(
    comparison = "Depth-slope",
    matched_fraction = 100 * mean(env.grid$n_match_sd > 0),
    matched_by_2plus_images = 100 * mean(env.grid$n_match_sd >= 2)),
  tibble(
    comparison = "Depth-slope-N-E-VRM",
    matched_fraction = 100 * mean(env.grid$n_match_5var > 0),
    matched_by_2plus_images = 100 * mean(env.grid$n_match_5var >= 2))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# EFFECT SIZE image coverage hotspot vs non (depth-slope)
env.grid %>% wilcox_effsize(
  n_match_sd ~ zone, comparisons = list(
    c("Hotspot", "High-density cluster"),
    c("Hotspot", "Non-Hotspot"),
    c("High-density cluster", "Non-Hotspot")))

# EFFECT SIZE image coverage hotspot vs non 5vars
env.grid %>% wilcox_effsize(
  n_match_5var ~ zone, comparisons = list(
    c("Hotspot", "High-density cluster"),
    c("Hotspot", "Non-Hotspot"),
    c("High-density cluster", "Non-Hotspot")))  


#################################################################
 ## DEPTH-RESOLVED SUMMARY + STANDING STOCK
#################################################################

env.coverage <- env.grid %>%
  mutate(depth_bin = cut(depth, breaks = seq(0, 3500, 500),include.lowest=T)) %>%
  group_by(depth_bin) %>%
  summarise(
    sampled_cells = n(),
    matched_sd_percent = 100 * mean(n_match_sd > 0),
    matched_5var_percent = 100 * mean(n_match_5var > 0),
    matched_2plus_sd_percent = 100 * mean(n_match_sd >= 2),
    matched_2plus_5var_percent = 100 * mean(n_match_5var >= 2),
    median_matched_5var = median(n_match_5var[n_match_5var > 0], na.rm=T), .groups="drop") %>%
  left_join(
    bioPredict %>%
      mutate(depth_bin = cut(depthReal, breaks = seq(0, 3500, 500), include.lowest=T)) %>%
      group_by(depth_bin) %>%
      summarise(area_km2 = n() * 625 / 1e6, total_density = sum(total_density * 625, na.rm=T), .groups="drop") %>%
      mutate(standingStock_percent = 100 * total_density / sum(total_density, na.rm=T)),  by="depth_bin") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Plot
env.coverage %>% rename(
  "Standing stock" = standingStock_percent,
  "Depth-slope: ≥1 image" = matched_sd_percent,
  "Depth-slope: ≥2 images" = matched_2plus_sd_percent,
  "Depth-slope-N-E-VRM: ≥1 image" = matched_5var_percent,
  "Depth-slope-N-E-VRM: ≥2 images" = matched_2plus_5var_percent) %>%
  pivot_longer(cols = c(
    "Standing stock", 
    "Depth-slope: ≥1 image","Depth-slope: ≥2 images",  
    "Depth-slope-N-E-VRM: ≥1 image", 
    "Depth-slope-N-E-VRM: ≥2 images"),
    names_to = "metric", values_to = "percent") %>%
  mutate(
    depth_bin = fct_rev(depth_bin),
    metric = factor(metric, levels = c(
      "Depth-slope: ≥2 images", 
      "Depth-slope: ≥1 image", "Depth-slope-N-E-VRM: ≥1 image",
      "Depth-slope-N-E-VRM: ≥2 images","Standing stock"))) %>%
  arrange(metric) %>%
  ggplot(aes(
    x = depth_bin,y = percent)) +
  geom_line(aes(group = metric, color = metric), linewidth = 1.5) +
  # Draw points second (top) to ensure they are never obscured
  geom_point(aes(group = metric, color = metric), size = 3) +
  coord_flip() +
  scale_x_discrete(labels = c(
    "3000–3500 m","2500–3000 m","2000–2500 m","1500–2000 m",
    "1000–1500 m","500–1000 m","0–500 m")) +
  scale_y_continuous(
    position = "right", limits = c(0, 100)) +
  scale_color_manual(values = c(
    "Standing stock"="saddlebrown",
    "Depth-slope: ≥1 image"="plum4", 
    "Depth-slope: ≥2 images"="plum", 
    "Depth-slope-N-E-VRM: ≥1 image"="skyblue3", 
    "Depth-slope-N-E-VRM: ≥2 images"="skyblue2")) +
  labs(
    x = NULL,  y = "Percent", color = NULL, shape = NULL) +  
  theme_classic() +
  theme(panel.grid.major = element_line(linewidth = 0.2, color="grey84"))


#################################################################
 ## COMBINED PLOT
#################################################################

# Plotted together with uncertainty map
# Same extent for all maps
map.x <- range(env.grid$x, na.rm = TRUE)
map.y <- range(env.grid$y, na.rm = TRUE)

# Uncertainty map
plot8 <- uncert.map %>%
  st_as_sf(coords = c("x", "y"), crs = 3031) %>%
  ggplot() +
  geom_sf(aes(color = uncertainty), size = 0.15, shape = 15) +
  scale_color_viridis_c(
    option = "magma",
    name = "Relative uncertainty",
    limits = quantile(uncert.map$uncertainty, c(0.02, 0.98), na.rm = TRUE),
    oob = scales::squish,
    guide = guide_colorbar(
      direction = "horizontal",
      barwidth = 10,
      barheight = 0.5,
      title.position = "bottom",
      title.hjust = 0.5)) +
  coord_sf(xlim = map.x, ylim = map.y, expand=F) +
  labs(x=NULL, y=NULL, title="Model reliability") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face="bold"),
    plot.margin = margin(5, 5, 5, 5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "bottom")

# Terrain matching maps
plot9 <- bind_rows(
  env.grid %>% transmute(x, y, comparison="Depth-slope", n_images = n_match_sd),
  env.grid %>% transmute(x, y, comparison="Depth-slope-N-E-VRM", n_images = n_match_5var)) %>%
  mutate(
    image_support = case_when(
      n_images == 0 ~ "None",
      n_images == 1 ~ "1",
      n_images <= 5 ~ "2–5",
      TRUE ~ ">5"),
    image_support = factor(image_support, levels = c("None","1","2–5",">5")),
    comparison = factor(comparison, levels = c("Depth-slope","Depth-slope-N-E-VRM"))) %>%
  st_as_sf(coords = c("x", "y"), crs=3031) %>%
  ggplot() +
  geom_sf(aes(color = image_support), size = 0.12, shape = 15) +
  facet_wrap(~comparison, nrow = 1) +
  scale_color_manual(
    values = c("None"="grey75","1"="mediumpurple2","2–5"="purple4",">5"="grey4"),
    name = "Images with matching terrain",
    guide = guide_legend(
      override.aes = list(size = 5),
      nrow = 1,
      title.position = "bottom",
      title.hjust = 0.5)) +
  coord_sf(xlim = map.x, ylim = map.y, expand=F) +
  labs(x = NULL, y = NULL, title = "Terrain matches") +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(5, 5, 5, 5),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9))

# Combine directly
plot_grid(
  plot8,
  plot9,
  nrow = 1,
  rel_widths = c(1, 2),
  align = "h",
  axis = "tb")


#################################################################
 ## DENSITY COMPARISON OFOBS-GRID 
#################################################################

# Add observed densities to env.img
env.bio <- env.img %>%
  left_join(counts.coral %>% select(image, obs_corals = Corals)) %>%
  left_join(counts.demo %>% select(image, obs_demosponges = Demosponges)) %>%
  left_join(counts.glass %>% select(image, obs_glass = GlassSponges)) %>%
  left_join(counts.ophi %>% select(image, obs_ophiuroids = Ophiuroids)) %>%
  left_join(counts.pens %>% select(image,obs_seapens = SeaPens)) %>%
  mutate(obs_sponges = obs_demosponges + obs_glass)

# Function for calculation
terrainDensity <- function(grid, img, tolerances, density_vars) {
  n_matches <- integer(nrow(grid))
  density_sum <- matrix(
    0, nrow = nrow(grid), ncol = length(density_vars), 
    dimnames = list(NULL, density_vars))
  for(i in seq_len(nrow(img))) {
    inside <- rep(TRUE, nrow(grid))
  for(v in names(tolerances)) {
    inside <- inside & abs(grid[[v]] - img[[v]][i]) <= tolerances[[v]]}
    inside[is.na(inside)] <- FALSE
    n_matches <- n_matches + as.integer(inside)
  for(v in density_vars) {
    density_sum[inside, v] <- density_sum[inside, v] + img[[v]][i]}}
    out <- tibble(n_matching_images = n_matches)
  for(v in density_vars) {
    out[[paste0(v, "_analogue_mean")]] <- ifelse(
    n_matches > 0, density_sum[, v] / n_matches, NA_real_)}
  return(out)}

# Create dataframe
env.bio <- env.grid %>% bind_cols(
  terrainDensity(
    grid = ., img = env.bio,
    tolerances = c(
      depth = 100, slope = unname(terrain.tol["slope"]),
      vrm = unname(terrain.tol["vrm"]), east = unname(terrain.tol["east"]), north = unname(terrain.tol["north"])),
      density_vars = c("obs_corals","obs_sponges","obs_ophiuroids","obs_seapens"))) %>%
  dplyr::select(x, y, depth, slope, vrm, east, north, everything()) %>%
  left_join(
    bioPredict %>% transmute(
      x, y, 
      pred_corals = corals, 
      pred_sponges = demosponges + glassSponges,
      pred_ophiuroids = ophiuroids, 
      pred_seapens = seaPens), by = c("x","y")) %>%
  mutate(
    support = case_when(
      n_matching_images == 0 ~"None", n_matching_images == 1 ~"1 image",
      n_matching_images <= 5 ~"2–5 images", TRUE ~ ">5 images"), 
    support = factor(support, levels = c("None","1 image","2–5 images",">5 images")))

# Calculate correlations
env.cor <- tibble(
  taxon = c("Corals","Sponges","Ophiuroids","Sea pens"),
  rho = c(
    cor(env.bio$pred_corals, env.bio$obs_corals_analogue_mean, method="spearman", use="complete.obs"),
    cor(env.bio$pred_sponges, env.bio$obs_sponges_analogue_mean, method="spearman", use="complete.obs"),
    cor(env.bio$pred_ophiuroids, env.bio$obs_ophiuroids_analogue_mean, method="spearman", use="complete.obs"),
    cor(env.bio$pred_seapens, env.bio$obs_seapens_analogue_mean, method="spearman", use="complete.obs"))) %>%
  mutate(label = paste("rho ==", round(rho, 2)))

# Plot by taxon
set.seed(12)
env.bio %>%
  filter(n_matching_images > 0) %>%
  transmute(
    support,
    Corals_predicted = pred_corals,
    Corals_observed = obs_corals_analogue_mean,
    Sponges_predicted = pred_sponges,
    Sponges_observed = obs_sponges_analogue_mean,
    Ophiuroids_predicted = pred_ophiuroids,
    Ophiuroids_observed = obs_ophiuroids_analogue_mean,
    `Sea pens_predicted` = pred_seapens,
    `Sea pens_observed` = obs_seapens_analogue_mean) %>%
  pivot_longer(
    cols = -support,
    names_to = c("taxon", ".value"),
    names_pattern = "(.+)_(predicted|observed)") %>%
  filter(predicted > 0, observed > 0) %>%
  slice_sample(n = 3000, by = c("support","taxon")) %>%
  ungroup() %>%
  arrange(support) %>%
  ggplot(aes(x = predicted, y = observed, color = support)) +
  geom_point(size=0.8, alpha=0.5) +
  geom_abline(intercept=0, slope=1, linewidth=0.8, linetype=2, color="red") +
  geom_text(
    data = . %>% 
      group_by(taxon) %>% 
      summarize(
        x = 10^(log10(min(predicted))), 
        y = 10^(log10(max(observed)))) %>%
      left_join(env.cor, by = "taxon"),
    aes(x = x, y = y, label = label),
    parse=T, hjust = -0.05, vjust=1, size=3, color="black", inherit.aes=F) +
  facet_wrap(~ taxon, scales = "free") +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  labs(
    x = expression("Predicted density (individuals " * m^-2 * ")"),
    y = expression("Observed density (individuals " * m^-2 * ")"))  +
  scale_color_manual(
    values = c("1 image"="mediumpurple2","2–5 images"="purple3",">5 images"="gray22"),
    name = "Matching terrain",
    guide = guide_legend( override.aes = list(size = 5, alpha=1))) +
  theme_classic() 


#################################################################
 ## COVERAGE BY DIVE
#################################################################

# Function to calculate
terrainDive <- function(grid, img, tolerances) {
    matched <- rep(FALSE, nrow(grid))
 for(i in seq_len(nrow(img))) {
    inside <- rep(TRUE, nrow(grid))
 for(v in names(tolerances)) {
      inside <- inside &
        abs(grid[[v]] - img[[v]][i]) <= tolerances[[v]]}
      inside[is.na(inside)] <- FALSE
      matched <- matched | inside} 
    matched}

# Summary
purrr::map_dfr(
  unique(env.img$dive), \(d) {
    img.dive <- env.img %>% filter(dive == d)
    matched <- terrainDive(
      grid = env.grid,
      img = img.dive,
      tolerances = c(
        depth = 100,
        slope = unname(terrain.tol["slope"]),
        vrm = unname(terrain.tol["vrm"]),
        east = unname(terrain.tol["east"]),
        north = unname(terrain.tol["north"])))
   tibble(
      dive = d,
      n_images = nrow(img.dive),
      matched_cells = sum(matched),
      matched_percent = 100 * mean(matched),
      matched_percent_perImage = 100 * mean(matched) / nrow(img.dive))} %>%
   mutate(
    dive = factor(dive, levels = c("39-1","69-1","81-1")),
    across(where(is.numeric), ~ round(.x, 3))))

# Plot
bind_rows(
  tibble(
    comparison = rep(c("Depth–slope", "Depth–slope-E-N-VRM"), each = 2),
    result = "Overall", support = rep(c("≥1 image", "≥2 images"), 2),
    dive = NA_character_, n_images = NA_real_,
    matched_percent = c(
      100 * mean(env.grid$n_match_sd >= 1),
      100 * mean(env.grid$n_match_sd >= 2),
      100 * mean(env.grid$n_match_5var >= 1),
      100 * mean(env.grid$n_match_5var >= 2))),
  purrr::map_dfr(c("Depth–slope", "Depth–slope-E-N-VRM"), \(v) {
    purrr::map_dfr(unique(env.img$dive), \(d) {
      img.dive <- env.img %>% filter(dive == d)
      matched <- terrainDive(
        grid = env.grid,
        img = img.dive,
        tolerances = if(v == "Depth–slope") {
          c(depth = 100, slope = unname(terrain.tol["slope"]))} else {
          c(depth = 100, slope = unname(terrain.tol["slope"]), vrm = unname(terrain.tol["vrm"]), east = unname(terrain.tol["east"]), north = unname(terrain.tol["north"]))})
      tibble(
        comparison = v, result = "Dive",
        support = NA_character_,  dive = d, n_images = nrow(img.dive),
        matched_percent = 100 * mean(matched))})})) %>%
mutate(
    comparison = factor(comparison, levels = c("Depth–slope-E-N-VRM","Depth–slope")),
    support = factor(support, levels = c("≥1 image", "≥2 images")),
    dive = factor(dive, levels = c("39-1", "69-1", "81-1")),
    y = as.numeric(comparison),
    y_plot = case_when(
      result=="Overall" & support =="≥1 image" ~ y + 0.20,
      result=="Overall" & support =="≥2 images" ~ y - 0.20,
      result == "Dive" ~ y)) %>%
ggplot() +
geom_segment(
  data = \(x) x %>% filter(result == "Overall"),
  aes(x = 0, xend = matched_percent, y = y_plot,
  yend = y_plot, linetype = support), linewidth=0.8, color="grey22") +
geom_point(
  data = \(x) x %>% filter(result == "Overall"),
  aes(x = matched_percent, y = y_plot), shape=15, size = 3.8, color = "grey15") +
geom_text(
  data = \(x) x %>% filter(result == "Overall"),
  aes(x = matched_percent,  y = y_plot, label = paste0(round(matched_percent, 0), "%")),
  nudge_x = 1.3,  hjust = 0, fontface = "bold", size = 3.2) +
geom_point(
  data = \(x) x %>% filter(result == "Dive"),
  aes(x = matched_percent, y = y_plot, fill = dive),
  shape = 21, color = "grey20",  size=3, stroke = 0.5, alpha = 0.95) +
scale_y_continuous(
  breaks = c(1, 2),
  labels = c("Depth–slope–E–N–VRM", "Depth–slope"),
  limits = c(0.7, 2.3), # Changed limits to create the gap at the bottom
  expand = expansion(mult = c(0, 0.02))) +
scale_x_continuous(
  limits = c(0, 70),
  breaks = seq(0, 70, 10),
  expand = expansion(mult = c(0, 0.02))) +
scale_linetype_manual(values = c(
  "≥1 image"="solid","≥2 images"="dashed")) +
scale_fill_manual(values = c(
  "39-1"="chartreuse","69-1"="olivedrab4","81-1"="seagreen2"), name="Dive") +
labs(x = "Regional grid cells with matching terrain (%)", y = NULL) +
theme_classic() +
theme(
  plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
  legend.position = "top",
  legend.justification = "left",
  legend.box = "vertical",
  legend.box.just = "left",
  axis.text.x = element_text(margin = margin(t = 2)), # Minimal push
  axis.text.y = element_text(face = "bold", margin = margin(r = 8)),
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank()) +
guides(
    fill = guide_legend(title = "Dive", nrow = 1, byrow=T),
    linetype = guide_legend(title = "Matches", nrow = 1, byrow=T),
    shape = "none")
