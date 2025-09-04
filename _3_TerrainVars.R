
###############################################################
  ## EFFECT OF LANDFORM AND RESOLUTION
###############################################################

## Which parameters differ between landforms?
## At which resolution is the effect strongest?

# Run Dunn's test
dunn.env <- meta %>%
  pivot_longer(cols = where(is.numeric), names_to="variable", values_to="value") %>%
  mutate(
    # Extract resolution at the END
    resolution = str_extract(variable,"(10cm|20cm|1m|2m)$"),
    # Remove resolution from variable name
    var_base = str_remove(variable,"_(10cm|20cm|1m|2m)$")) %>%
  filter(!is.na(resolution)) %>%
  filter(!variable %in% c("area")) %>%
  group_by(variable) %>%
  dunn_test(value~landform, p.adjust.method="BH") %>%
  ungroup() %>%
  mutate(
    resolution = str_extract(variable,"(10cm|20cm|1m|2m)$"),
    var_base = str_remove(variable,"_(10cm|20cm|1m|2m)$"),
    sig = p.adj < 0.05)

# Summarize for each var + resolution
dunn.env <- dunn.env %>%
  group_by(var_base, resolution) %>%
  summarise(
    n_sig = sum(p.adj < 0.05),
    n_total = n(),
    median_padj = if_else(n_sig > 0, median(p.adj[p.adj < 0.05]), NA_real_),
    min_padj = if_else(n_sig > 0, min(p.adj[p.adj < 0.05]), NA_real_),
    frac_sig = n_sig / n_total, .groups="drop")

# Determine best resolution per variable:
# Max number of significant pairs (n_sig) and smallest median p.adj 
bestRes <- dunn.env %>%
  group_by(var_base) %>%
  arrange(median_padj, desc(n_sig)) %>%  # select by better pvalues
  #arrange(desc(n_sig), median_padj) %>%  # select by no. correlations
  slice(1) %>%
  ungroup() %>%
  mutate(best=T)
# After downstream modelling: best to prioritize pvalues 
# Higher model significance; lower VIF

# Join for plotting
dunn.res <- dunn.env %>%
  left_join(bestRes %>% dplyr::select(var_base, bestResolution = resolution), by="var_base") %>%
  mutate(isBest = resolution == bestResolution) %>%
  dplyr::select(var_base, resolution, n_sig, n_total, frac_sig, median_padj, min_padj, isBest) 

# Plot summary; gold diamonds for bestRes
# Export size 3.5 x 5.5
dunn.res %>% 
  filter(!is.na(n_sig), !is.na(median_padj)) %>% 
  mutate(var_base = factor(var_base, levels = c(
    "roughness","slope","aspect","planCurvature","profileCurvature")),
  resolution = factor(resolution, levels = c(
    "10cm","20cm","1m","2m"))) %>%
  ggplot(aes(x = var_base, y = resolution)) +
  geom_point(
    data = . %>% filter(!isBest),
    aes(size = rev(median_padj)),
    shape = 21, fill="grey64", color="black", stroke=0.2, alpha=0.96) +
  geom_point(
    data = . %>% filter(isBest),
    aes(size = rev(median_padj)),
    shape = 23, fill="gold", color="black", stroke=1) +
  scale_size(
    range = c(4, 7),
    breaks = c(0.0001, 0.005, 0.01, 0.05),
    name="Median p-adj") +
    labs(x=NULL, y="Resolution") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major = element_line(color="grey85"),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right")

###################################################

# Final test: extract vars at bestRes
# omit area: no ecological meaning, and may confound
dunnVars <- c(
  dunn.res %>%
    filter(isBest) %>%
    mutate(variable = paste(var_base, resolution, sep="_")) %>%
    pull(variable),
  meta %>%
    dplyr::select(where(is.numeric), -area) %>%
    names() %>%
    .[!str_detect(., "_")])

# Run Dunn
dunn.env <- meta %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "variable",
    values_to = "value") %>%
  filter(variable %in% dunnVars) %>%
  mutate(
    resolution = str_extract(variable, "(10cm|20cm|1m|2m)$"),
    var_base = str_remove(variable, "_(10cm|20cm|1m|2m)$")) %>%
  group_by(variable) %>%
  dunn_test(value ~ landform, p.adjust.method="BH") %>%
  ungroup() %>%
  filter(p.adj < 0.05)

# Plot variables per landform -- Figure X
# Export size 6.5 x 3
meta %>% 
  pivot_longer(cols = where(is.numeric), names_to="variable", values_to="value") %>%
  filter(variable %in% dunnVars & !variable %in% c("lat","lon","area","depth")) %>%
  mutate(
    landform = factor(landform, levels = c(
      "Terrace","Steep slope","Ridges","Broad slope",
      "Depression","Flat area","Sand ripples")),
    variable = str_remove(variable, "(_10cm|_20cm|_1m|_2m)"),
    variable = factor(variable, levels=c(
      "slope","roughness","aspect","planCurvature",
      "profileCurvature","depth"))) %>%
  ggplot(aes(x = landform, y = value, fill = landform)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  geom_jitter(
    aes(x=landform, y=value, fill=landform),
    shape=21, color="gray16", stroke=0.04,
    size=1.3, alpha=0.8, width=0.1, show.legend=F) +
  facet_grid(variable ~ ., scales = "free") +
  scale_fill_manual(values = landformCol) +
  scale_color_manual(values = landformCol) +
  scale_y_continuous(n.breaks = 3.4) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=10),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color="gray16", linewidth=0.2), 
    panel.spacing = unit(1, "lines"), 
    legend.position = "none",
    panel.grid.major = element_line(color="gray87", linewidth=0.2),
    axis.text.x = element_text(angle = 45, hjust = 1))


###############################################################
  ## LANDFORM-TERRAIN VARIABLITY: PCA
###############################################################

# Select numeric vars
geoVars <- meta %>%
  dplyr::select(all_of(dunnVars)) %>%
  drop_na()

# Perform PCA
pca.geo <- PCA(geoVars, scale.unit=T, graph=F)

# Extract scores
scores <- as.data.frame(pca.geo$ind$coord) %>%
  rownames_to_column("image") %>%
  left_join(meta %>% rownames_to_column("image") %>% dplyr::select(image, landform, dive)) %>%
  column_to_rownames("image")

# Extract coordinates; scale to make arrows more visible
# Slight offset beyond arrow tips for labels
max_r <- max(abs(c(scores$Dim.1, scores$Dim.2)))
loadings <- as.data.frame(pca.geo$var$coord) %>%
  rownames_to_column("varnames") %>%
  mutate(
    Dim.1 = Dim.1 * max_r * 0.7,
    Dim.2 = Dim.2 * max_r * 0.7,
    label_x = Dim.1 * 1.05,
    label_y = Dim.2 * 1.05)

# Plot -- Export size 5x6
ggplot() +
  geom_point(
    data = scores,
    aes(x = Dim.1, y = Dim.2, color=landform, shape=dive),
    size = 4, alpha = 0.8) +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    arrow = arrow(length = unit(0.25,"cm"), type="closed"),
    color = "black") +
  geom_text_repel(
    data = loadings,
    aes(x = label_x, y = label_y, label = str_remove(varnames, "_.*")),
    size = 4, color = "black",
    min.segment.length = 0, max.overlaps = Inf) +
  labs(
    x = paste0("PC1 (", round(pca.geo$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(pca.geo$eig[2, 2], 1), "%)")) +
  scale_color_manual(values = landformCol) +
  theme_classic() +
  theme(legend.position = "right")

# Determine variable contributions; export  (%)
as.data.frame(pca.geo$var$contrib) %>%
  rownames_to_column("variable") %>%
  pivot_longer(-variable, names_to="PC", values_to="Contribution (%)") %>%
  filter(PC %in% c("Dim.1","Dim.2")) %>%
  arrange(PC, desc(`Contribution (%)`)) %>%
  write.table("res_PCA.txt", sep="\t", row.names=F, quote=F)

