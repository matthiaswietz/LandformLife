
#################################################################
 ## SUMMARY STATS
#################################################################

# Dive summary
meta %>%
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

# Landform summary
meta %>%
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
 ## ATTRIBUTE STATISTICS
###############################################################

## CONSISTENCY ACROSS RESOLUTIONS
attributes.res <- meta %>%
  select(matches("_(10cm|20cm|1m|2m|25m)$")) %>%
  drop_na() %>% {dat <- .
    vars <- colnames(dat)
      expand.grid(var1 = vars[grepl("_25m$", vars)], var2 = vars, stringsAsFactors=F) %>%
      filter(var1 != var2) %>%
      filter(sub("_[^_]+$", "", var1) == sub("_[^_]+$", "", var2)) %>%
  purrr::pmap_dfr(function(var1, var2) {
    test <- cor.test(dat[[var1]], dat[[var2]], method="spearman", exact=F)
    data.frame(
      variable = sub("_[^_]+$", "", var1),
      resolution = paste0(sub(".*_", "", var1), " vs ", sub(".*_", "", var2)),
      spearman = as.numeric(test$estimate),
      p = as.numeric(test$p.value))})} %>%
  mutate(pAdj = p.adjust(p, method="BH")) %>%
  arrange(desc(abs(spearman))) %>%
  select(variable, resolution, spearman, pAdj)

## COLLINEARITY
attributes.vif <- meta %>%
  select(matches("_(10cm|20cm|1m|2m|25m)$")) %>%
  drop_na() %>% {dat <- . 
  res_list <- c("10cm", "20cm", "1m", "2m", "25m")
    purrr::map_dfr(res_list, function(r) {
      cols <- grep(paste0("_", r, "$"), colnames(dat), value = TRUE)
      temp_df <- dat %>% select(all_of(cols))
      clean_names <- sub(paste0("_", r, "$"), "", colnames(temp_df))
      colnames(temp_df) <- clean_names
      fit <- lm(seq_len(nrow(temp_df)) ~ ., data = temp_df)
      v_vals <- car::vif(fit)
      data.frame(
        resolution = r,
        variable = names(v_vals),
        VIF = as.numeric(v_vals))}) } %>%
  arrange(desc(VIF))

# DISTINCTION OF LANDFORMS
attributes.sel <- attributes.res %>%
  filter(grepl("25m vs", resolution)) %>%
  mutate(res_only = sub("25m vs ", "", resolution)) %>%
  left_join(attributes.vif, by = c("variable"="variable","res_only"="resolution")) %>%
  select(variable, res_only, spearman, VIF) %>%
  rename(
    resolution = res_only,
    spearman_25m = spearman,
    attribute = variable,
    #landform_sep_frac = frac_sig,
    #landform_sep_p = median_p
    ) %>%
  arrange(desc(spearman_25m))

# Print outcome
attributes.sel

####################################################

## BASED ON THIS: 20cm best; roughness removed
## Which attributes differ between landforms?
## Which resolution shows the strongest effect?
dunn.geo <- meta %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "variable",
    values_to = "value") %>%
  mutate(
    resolution = str_extract(variable, "(10cm|20cm|1m|2m)$"),
    metric = str_remove(variable, "_(10cm|20cm|1m|2m)$")) %>%
  filter(!is.na(resolution), !variable %in% c("area")) %>%
  group_by(variable) %>%
  dunn_test(value~landform, p.adjust.method="BH") %>%
  ungroup() %>%
  mutate(
    resolution = str_extract(variable, "(10cm|20cm|1m|2m)$"),
    metric = str_remove(variable, "_(10cm|20cm|1m|2m)$"),
    sig = p.adj < 0.05)

# Summarize per variable + resolution
attributes.lf <- dunn.geo %>%
  group_by(metric, resolution) %>%
  summarise(
    n_sig     = sum(sig),
    n_total   = n(),
    frac_sig  = n_sig / n_total,
    median_p  = if_else(n_sig > 0, median(p.adj[sig]), NA_real_),
    .groups   = "drop") %>%
  arrange(desc(frac_sig), median_p)

# Print outcome >> Table SX
print(attributes.lf, n=Inf)


#################################################################
  ## LANDFORM-TERRAIN VARIABILITY: PCA
#################################################################

# Select variables
geoVars <- meta %>%
  dplyr::select(lat, lon, depth, tidyselect::contains("_20cm")) %>%
  drop_na() 

# Perform PCA
pca.geo <- PCA(geoVars, scale.unit=T, graph=F)

# Extract scores
scores.pca <- as.data.frame(pca.geo$ind$coord) %>%
  rownames_to_column("image") %>%
  left_join(meta %>% rownames_to_column("image") %>% select(image, landform, dive)) %>%
  column_to_rownames("image")

# Extract coordinates; scale to make arrows more visible
# Slight offset beyond arrow tips for labels
arrowScale <- max(abs(c(scores.pca$Dim.1, scores.pca$Dim.2))) *0.5

# Adjust loadings for plot; remove vars thjat contribute <0.5% in both dimensions
loadings <- as.data.frame(pca.geo$var$coord) %>%
  rownames_to_column("varnames") %>%
  mutate(
    contrib_total = pca.geo$var$contrib[, "Dim.1"] + pca.geo$var$contrib[, "Dim.2"]) %>%
  filter(contrib_total >= 0.5) %>% 
  mutate(
    Dim.1 = Dim.1 * arrowScale,
    Dim.2 = Dim.2 * arrowScale,
    label_x = Dim.1 * 1.05,
    label_y = Dim.2 * 1.05)

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
    arrow = arrow(length = unit(0.3,"cm"), type="closed"),
    color="black") +
  geom_text_repel(
    data = loadings,
    aes(x = label_x, y = label_y, label = str_remove(varnames, "_[0-9]+(cm|m)")),
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
      "slope","VRM","roughness","aspect_E","aspect_N","aspect",
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

