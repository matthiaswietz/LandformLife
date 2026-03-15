
#################################################################
 ## SUMMARY STATS
#################################################################

# How many landform images in total?
meta %>% group_by(landform) %>% summarise(n_images = n()) 

# How many images per dive?
meta %>% group_by(dive) %>% summarise(n_images = n()) 

# Dive summary
overview.dive <- meta %>%
  group_by(dive) %>%
  summarise(
    n_images = n(),
    depth_mean = mean(depth, na.rm=T),
    depth_sd   = sd(depth, na.rm=T),
    depth_min  = min(depth, na.rm=T),
    depth_max  = max(depth, na.rm=T),
    lat_mean   = mean(lat, na.rm=T),
    lon_mean   = mean(lon, na.rm=T),
    .groups = "drop")

# Landform  summary
overview.lf <- meta %>%
  group_by(landform) %>%
  summarise(
    n_images = n(),
    depth_median = median(depth),
    depth_IQR = IQR(depth),
    depth_min = min(depth),
    depth_max = max(depth),
    .groups = "drop")


###############################################################
 ## IMAGE CLASSIFICATION ACCURACY
###############################################################

# Prepare data
confMatrix <- data.frame(
  Actual = factor(c(
    "Flat areas", "Broad slopes","Steep slopes", "Ridge", 
    "Depression", "Sand ripples","Terrace"), levels = c(
      "Flat areas","Broad slopes","Steep slopes", "Ridge", 
      "Depression","Sand ripples","Terrace")),
  `Flat areas` = c(70, 5, 1, 1, 0, 1, 0),
  `Broad slopes` = c(6, 70, 2, 1, 0, 0, 0),
  `Steep slopes` = c(0, 1, 52, 0, 0, 0, 3),
  Ridge = c(0, 0, 0, 18, 0, 0, 0),
  Depression = c(1, 1, 0, 0, 44, 0, 0),
  `Sand ripples` = c(0, 0, 0, 0, 0, 31, 0),
  Terrace = c(0, 4, 4, 0, 0, 0, 62),
  check.names = FALSE)

# Prepare for calculation accuracy
confValues <- as.matrix(confMatrix[, -1])
rownames(confValues) <- confMatrix$Actual
correct <- sum(diag(confValues))
total <- sum(confValues)

# Calculate
correct / total

# Print for SI table
confMatrix


###############################################################
 ## EFFECT OF LANDFORM AND RESOLUTION
###############################################################

## Which parameters differ between landforms?
## At which resolution is the effect strongest?
dunn.geo <- meta %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "variable",
    values_to = "value") %>%
  mutate(
    resolution = str_extract(variable, "(10cm|20cm|1m|2m)$"),
    metric   = str_remove(variable, "_(10cm|20cm|1m|2m)$")) %>%
  filter(!is.na(resolution), !variable %in% c("area")) %>%
  group_by(variable) %>%
  dunn_test(value~landform, p.adjust.method="BH") %>%
  ungroup() %>%
  mutate(
    resolution = str_extract(variable, "(10cm|20cm|1m|2m)$"),
    metric = str_remove(variable, "_(10cm|20cm|1m|2m)$"),
    sig = p.adj < 0.05)

# Summarize per variable + resolution
dunn.summary <- dunn.geo %>%
  group_by(metric, resolution) %>%
  summarise(
    n_sig     = sum(sig),
    n_total   = n(),
    frac_sig  = n_sig / n_total,
    median_p  = if_else(n_sig > 0, median(p.adj[sig]), NA_real_),
    .groups   = "drop") %>%
  arrange(desc(frac_sig), median_p)

# 20 cm confirmed as best resolution
# Print outcome and copy to Table SX_ResTest
dunn.summary


#################################################################
  ## LANDFORM-TERRAIN VARIABILITY: PCA
#################################################################

# Select variables
geoVars <- meta %>%
  dplyr::select(lat, lon, depth, tidyselect::ends_with("_20cm")) %>%
  drop_na() 

# Perform PCA
pca.geo <- PCA(geoVars, scale.unit=T, graph=F)

# Extract scores
scores.pca <- as.data.frame(pca.geo$ind$coord) %>%
  rownames_to_column("image") %>%
  left_join(meta %>% rownames_to_column("image") %>% dplyr::select(image, landform, dive)) %>%
  column_to_rownames("image")

# Extract coordinates; scale to make arrows more visible
# Slight offset beyond arrow tips for labels
max_r <- max(abs(c(scores.pca$Dim.1, scores.pca$Dim.2)))
arrow_scale <- max_r *0.5

# Adjust loadings for plot 
loadings <- as.data.frame(pca.geo$var$coord) %>%
  rownames_to_column("varnames") %>%
  mutate(
    Dim.1 = Dim.1 * arrow_scale,
    Dim.2 = Dim.2 * arrow_scale,
    label_x = Dim.1 *1.05,
    label_y = Dim.2 *1.05)

# PLOT -- FIG 3A
# combined with RDA / Fig 3b
plot1 = ggplot() +
  geom_point(
    data = scores.pca,
    aes(x = Dim.1, y = Dim.2, color = landform),
    size=2.44, alpha=0.8) +
  geom_segment(
    data = loadings,
    aes(x=0, y=0, xend = Dim.1, yend = Dim.2),
    arrow = arrow(length = unit(0.25,"cm"), type="closed"),
    color="black") +
  geom_text_repel(
    data = loadings,
    aes(x = label_x, y = label_y, label = str_remove(varnames, "_.*")),
    size=4, color="black",
    min.segment.length=0, max.overlaps=Inf) +
  labs(
    x = paste0("PC1 (", round(pca.geo$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(pca.geo$eig[2, 2], 1), "%)")) +
  scale_color_manual(values = landformCol) +
  scale_shape_manual(values = c(15, 17, 19)) +
  theme_classic() +
  theme(legend.position="right")

# Determine variable contributions; export
as.data.frame(pca.geo$var$contrib) %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to="PC", values_to="Contribution (%)") %>%
  filter(PC %in% c("Dim.1","Dim.2")) %>%
  arrange(PC, desc(`Contribution (%)`)) %>%
  write.table("./results/PCA.txt", sep="\t", row.names=F, quote=F)


###############################################################
 ## LANDFORM-TERRAIN VARIABILITY: BOXPLOT
###############################################################

# FIGURE 3C
# Export size 6.5 x 3
meta %>% 
  pivot_longer(cols = where(is.numeric), names_to="variable", values_to="value") %>%
  filter(variable %in% names(geoVars) & !variable %in% c("lat","lon","area","depth")) %>%
  mutate(
    landform = factor(landform, levels = c(
      "Terrace","Steep slope","Ridges","Broad slope",
      "Depression","Flat area","Sand ripples")),
    variable = str_remove(variable, "_20cm"),
    variable = factor(variable, levels = c(
      "slope","roughness","aspect",
      "profileCurvature","planCurvature"))) %>%
  ggplot(aes(x = landform, y = value, fill = landform)) +
  geom_boxplot(outlier.shape=NA, alpha=1) +
  geom_jitter(
    aes(x=landform, y=value, fill=landform),
    shape=21, color="gray44", stroke=0.05,
    size=1.6, alpha=0.8, width=0.2, show.legend=F) +
  facet_grid(variable ~ ., scales="free") +
  scale_fill_manual(values = landformCol) + 
  scale_color_manual(values = landformCol) +
  scale_y_continuous(n.breaks=3.6) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=10),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color="gray16", linewidth=0.2), 
    panel.spacing = unit(1, "lines"), 
    legend.position="none",
    panel.grid.major = element_line(color="gray87", linewidth=0.2),
    axis.text.x = element_text(angle=45, hjust=1))

