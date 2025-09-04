####################################################

## ALPHA-DIVERSITY AND RAREFACTION

# Remove taxa with any NA
# Keep samples with >1 observed taxa
# Keep columns with sum >0 
# Round counts: 0 < x < 1 → 1
counts.corr %>%
  filter(!if_any(everything(), is.na)) %>%
  dplyr::select(where(~ sum(. > 0, na.rm=T) > 1)) %>%
  mutate(across(everything(), ~ ifelse(. > 0.001 & . < 1, 1, round(.)))) %>%
  dplyr::select(where(~ sum(.) > 0)) %>%
  iNEXT(q=0, datatype="abundance", nboot=50) -> inext

# Compile outcome
richness <- inext$AsyEst[
  inext$AsyEst$Diversity=="Species richness",] %>%
  arrange(Assemblage) 
simpson <- inext$AsyEst[
  inext$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Assemblage) 
shannon <- inext$AsyEst[
  inext$AsyEst$Diversity=="Shannon diversity",] %>%
  arrange(Assemblage) 

# Compile; calculate evenness
AlphaDiv <- data.frame(
  data = richness$Assemblage,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Observed) %>%
  mutate(evenness = shannon/log(richness)) %>%
  left_join(
    meta %>% rownames_to_column("image"),
    by = c("data"="image")) %>%
  dplyr::rename(image=data)

# Test for significance
AlphaDivSign <- AlphaDiv %>%
  dunn_test(
    shannon ~ landform,
    p.adjust.method = "BH") %>%
  ungroup()

## PLOT WITH CLD
# Make symmetric p-value table
matSym <- bind_rows(
  AlphaDivSign %>% dplyr::select(group1, group2, p.adj = p.adj),
  AlphaDivSign %>% rename(group1 = group2, group2 = group1) %>% 
    dplyr::select(group1, group2, p.adj = p.adj)) %>%
  distinct()

# Get unique landforms from matSym for pivoting
landforms <- union(unique(AlphaDivSign$group1), unique(AlphaDivSign$group2))

# Make wide matrix 
matWide <- matSym %>%
  pivot_wider(names_from = group2, values_from = p.adj) %>%
  complete(group1 = landforms) %>%
  arrange(factor(group1, levels = landforms)) %>%
  dplyr::select(group1, all_of(landforms))

# Set correct order
matLevels <- c(
  "Terrace","Steep slope","Ridges","Broad slope",
  "Depression","Flat area","Sand ripples")

# Convert to matrix
matAll <- matWide %>%
  column_to_rownames("group1") %>%
  as.matrix()

# Reorder by factor levels
matAll <- matAll[matLevels, matLevels]

# Symmetrize matrix
matAll[upper.tri(matAll)] <- t(matAll)[upper.tri(matAll)]
matAll[is.na(matAll)] <- 1  # Treat missing as non-significant
diag(matAll) <- 0          # Diagonal = 0 (same group)

# Calculate compact letter display
sig.div <- multcompLetters(matAll)$Letters %>%
  enframe(name="landform", value="cld")

# Prepare labels: bit below min-y
labels.div <- AlphaDiv %>%
  group_by(landform) %>%
  summarise(y = min(shannon, na.rm=T) - 0.05 * diff(
    range(shannon))) %>%
  left_join(sig.div, by="landform")

# Plot; later combined with density and composition
plot1 = AlphaDiv %>%
  mutate(landform=factor(landform, levels=rev(c(
    "Sand ripples","Flat area","Depression",
    "Broad slope","Ridges","Steep slope","Terrace")))) %>%
  ggplot(aes(x = landform, y = shannon)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(
    width = 0.15, shape = 21, 
    fill = "grey8", color = "black", 
    alpha = 0.5, size = 1.6, stroke = 0.05) +
  geom_text(data = labels.div, aes(
    x = landform, y = y, label = cld), 
    color="darkred", size=5, inherit.aes=F) +
  labs(y="Shannon diversity") +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.ticks = element_blank(),
    axis.text.x = element_blank())


###############################################################
 ## Density differences between landforms?
###############################################################

## FIRST: Kruskal-Wallis test

# Prepare data
counts.image <- counts.corr %>%
  rownames_to_column("taxon") %>%
  pivot_longer(-taxon, names_to="image", values_to="density_m2") %>%
  left_join(meta %>% rownames_to_column("image"), by="image")

# Calculate density; scale to 10m2 (for later plotting)
counts.density <- counts.image %>%
  mutate(landform=factor(landform, levels=rev(c(
    "Sand ripples","Flat area","Depression","Broad slope",
    "Ridges","Steep slope","Terrace")))) %>%
  group_by(image, landform) %>%
  summarise(
    total.density_m2 = sum(density_m2, na.rm=T), 
    total.density_10m2 = sum(density_m2, na.rm=T)*10, .groups="drop") 

# Kruskal-Wallis test
kruskal.test(
  total.density_10m2 ~ landform, data=counts.density)

# Dunn's test 
dunn.landform = dunn_test(
  counts.density, total.density_10m2 ~ landform, p.adjust.method="BH")

## PLOT WITH CLD
# Make symmetric p-value table
matSym <- bind_rows(
  dunn.landform %>% dplyr::select(group1, group2, p.adj = p.adj),
  dunn.landform %>% rename(group1 = group2, group2 = group1) %>% 
    dplyr::select(group1, group2, p.adj = p.adj)) %>%
  distinct()

# Get unique landforms from matSym for pivoting
landforms <- union(unique(matSym$group1), unique(matSym$group2))

# Make wide matrix 
matWide <- matSym %>%
  pivot_wider(names_from = group2, values_from = p.adj) %>%
  complete(group1 = landforms) %>%
  arrange(factor(group1, levels = landforms)) %>%
  dplyr::select(group1, all_of(landforms))

# Set correct order
matLevels <- c(
  "Terrace","Steep slope","Ridges","Broad slope",
  "Depression","Flat area","Sand ripples")

# Convert to matrix
matAll <- matWide %>%
  column_to_rownames("group1") %>%
  as.matrix()

# Reorder by factor levels
matAll <- matAll[matLevels, matLevels]

# Symmetrize matrix
matAll[upper.tri(matAll)] <- t(matAll)[upper.tri(matAll)]
matAll[is.na(matAll)] <- 1  # Treat missing as non-significant
diag(matAll) <- 0          # Diagonal = 0 (same group)

# Calculate compact letter display
sig.dens <- multcompLetters(matAll)$Letters %>%
  enframe(name="landform", value="cld")

# Prepare labels: bit below min-y
label.dens <- counts.density %>%
  group_by(landform) %>%
  summarise(y = min(total.density_10m2, na.rm=T) - 0.05 * diff(
    range(total.density_10m2))) %>%
  left_join(sig.dens, by="landform")

## PLOT 
## later combined with diversity and composition
plot2 = counts.density %>%
  ggplot(aes(x = landform, y = total.density_10m2)) +
  geom_boxplot(fill="white", outlier.shape = NA) +
  geom_jitter(
    width = 0.15, shape = 21, 
    fill = "grey8", color = "black", 
    alpha = 0.5, size = 1.6, stroke = 0.05) +
  geom_text(data = label.dens, aes(
    x = landform, y = y, label = cld), 
    color="darkred", size=5, inherit.aes=F) +
  labs(y = expression("Density per 10" ~ m^{2})) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.ticks = element_blank(),
    axis.text.x = element_blank())

####################################################

## MODELLING OF TERRAIN VS DENSITIES
# Comparing landform+geoVars & geoVars only

## First check: correlations between geo-vars?
varCorr <- meta %>%
  dplyr::select(matches("^(aspect|planCurvature|profileCurvature|roughness|slope)_(10cm|20cm|1m|2m|25m)$")) %>%
  drop_na() %>%
  cor(method="spearman") %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to="var2", values_to="correlation") %>%
  filter(var1 < var2) %>%
  # Here extract resolution part from var1 and var2 accordingly:
  filter(
    sub(".*_(10cm|20cm|1m|2m|25m)$", "\\1", var1) !=
    sub(".*_(10cm|20cm|1m|2m|25m)$", "\\1", var2)) %>%
  arrange(desc(abs(correlation)))
# slope+roughness highly correlated
# VIF is OK (see below) -- will keep both with caution

############################

## PREPARE DENSITIY DATA
# Try different transformations
counts.image$logDensity <- log(counts.image$density_m2 + 1)
counts.image$sqrtDensity <- sqrt(counts.image$density_m2 + 1)

# Plotting shows: log-transform better 
# also coherent with earlier log-transform
counts.image %>%
  rename(`raw density`=density_m2) %>%
  pivot_longer(cols = c(
    `raw density`, sqrtDensity, logDensity), 
    names_to="transform", values_to="value") %>%
  mutate(transform = factor(transform, levels=c(
    "raw density","logDensity","sqrtDensity"))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill="blueviolet", color="black") +
  scale_y_continuous(expand=c(0.01, 0.01)) +
  facet_wrap(~transform, scales="free") +
  labs(x="Density", y="Frequency") +
  theme_minimal() +
  theme(axis.title = element_text(size=10))

############################

# Prepare model data
modelData <- counts.image %>%
  dplyr::select(logDensity, landform, any_of(dunnVars), where(is.numeric)) %>%
  dplyr::select(-matches("_"), any_of(dunnVars), -any_of(c("lat","lon","area","sqrtDensity"))) %>%
  mutate(landform = fct_relevel(factor(landform), "Flat area"))

# Prepare predictors
modelVars <- modelData %>%
  dplyr::select(where(is.numeric), -logDensity) %>%
  names()

# Check collinearity for terrain-only 
vif(lm(logDensity ~ ., data = modelData %>% dplyr::select(logDensity, all_of(modelVars))))

# --- 4a. GAM: terrain + landform ---
density.lf <- gam(
  reformulate(
    c("landform", paste0("s(", modelVars, ")")),
    response = "logDensity"),
  data = modelData,
  method = "REML")
summary(density.lf)

# --- 4b. GAM: terrain-only ---
density.geo <- gam(
  reformulate(paste0("s(", modelVars, ")"), response="logDensity"),
  data = modelData,
  method = "REML")
summary(density.geo)

############################

# Compute AIC 
AIC(density.lf, density.geo)
# lf 3984 / geo 4026

# Summary tables 
bind_rows(
  broom::tidy(density.lf, parametric=T) %>%
    filter(term != "(Intercept)") %>%
    mutate(term_group = "landform"),
  broom::tidy(density.lf, smooths=T) %>%
    mutate(term_group = "terrain")) %>%
  mutate(
    signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""),
    edf = ifelse(is.na(edf), "-", as.character(edf))) %>%
  dplyr::select(term, term_group, edf, statistic, p.value, signif) %>%
  arrange(p.value) %>%
  add_column(model="landform+terrain") %>% 
  write.table("res_gamLandform.txt", sep="\t", row.names=F, quote=F)

broom::tidy(density.geo) %>%
  filter(term != "(Intercept)") %>%
  mutate(signif = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "")) %>%
  dplyr::select(term, edf, statistic, p.value, signif) %>%
  arrange(p.value) %>%
  add_column(model="terrain-only", .before=1) %>%
  write.table("res_gamGeo.txt", sep="\t", row.names=F, quote=F)


###############################################################
 ## Compositional differences between landforms?
###############################################################

# Calculate distance matrix
dist <- counts.corr %>%
  t() %>% 
  as.data.frame() %>% 
  drop_na() %>%
  decostand(method = "log") %>%
  vegdist(method="bray")

# Log-transform counts
counts.log <- t(counts.corr) %>% as.data.frame() %>% 
  drop_na() %>% as.matrix() %>% 
  decostand(method="log")

# Extract original PCA scores 
pcaScores <- as.data.frame(pca.geo$ind$coord) %>% 
  dplyr::select(Dim.1, Dim.2)

# Vector to match all data
shared <- Reduce(intersect, list(
  rownames(pcaScores),
  rownames(as.matrix(dist)),
  rownames(counts.log)))

# Match samples with distance matrix
pcaScores <- pcaScores[shared, , drop=F]
dist <- as.dist(as.matrix(dist)[shared, shared])
counts.log <- counts.log[shared, , drop=F]

# PERMANOVA: community vs PCA axes
adonis2(as.dist(as.matrix(dist)[
  shared, shared]) ~ Dim.1 + Dim.2, 
  data = pcaScores[shared, , drop=F])

############################

## STANDARD RDA

# Prepare env variables
# Set 'flat area' as reference level
env.rda <- meta %>%
  dplyr::select(all_of(dunnVars), landform) %>%
  mutate(landform = fct_relevel(factor(landform), "Flat area"))

# Dummy-code landforms
landforms <- model.matrix(~ landform, data = env.rda)[, -1]  # remove intercept

# Combine
env.mat <- cbind(env.rda[dunnVars], landforms) %>%
  drop_na()

# Check if cols have zero variance
varZero <- apply(env.mat, 2, var, na.rm=T)
print(varZero)
# none -- no filtering needed

# Match to counts
counts.sub <- counts.log[rownames(
  counts.log) %in% rownames(env.mat), , drop=F]
env.mat <- env.mat[rownames(counts.sub), , drop=F]

# Run RDA
rda <- rda(counts.sub ~ ., data = as.data.frame(env.mat))

# Significance testing
anova(rda)
anova(rda, by="terms")
summary(rda)

# Extract sites + env vectors
scores.rda <- scores(rda, display="sites", scaling=2) %>% as.data.frame()
vectors.rda <- scores(rda, display="bp", scaling=2) %>% as.data.frame()
scores.rda$Sample <- rownames(scores.rda)
scores.rda <- left_join(scores.rda, meta %>% mutate(Sample=rownames(meta)), by="Sample")

# Extract stats
eig.rda <- summary(rda)$cont$importance[2, 1:2] * 100

# Plot -- export size 4x5
ggplot(scores.rda, aes(x = RDA1, y = RDA2, color = landform)) +
  geom_point(size=4, alpha=0.7) +
  scale_color_manual(values=landformCol) +
  labs(x = paste0("RDA1 (", round(eig.rda[1], 1), "%)"),
       y = paste0("RDA2 (", round(eig.rda[2], 1), "%)")) +
  theme_classic() +
  theme(axis.ticks = element_blank())

####################################

## dbRDA

# Bray-Curtis distance
counts.dist <- vegdist(counts.sub, method="bray")

# dbRDA with env vectors
dbrda <- capscale(
  counts.dist ~ ., data = as.data.frame(env.mat))

# ANOVA
anova(dbrda)
anova(dbrda, by="terms")
summary(dbrda)

# Extract sites + env vectors
scores.dbrda <- scores(dbrda, display="sites", scaling=2) %>% as.data.frame()
vectors.dbrda <- scores(dbrda, display="bp", scaling=2) %>% as.data.frame()
scores.dbrda$Sample <- rownames(scores.dbrda)
scores.dbrda <- left_join(scores.dbrda, meta %>% mutate(Sample=rownames(meta)), by="Sample")

# Stats
eig.dbrda <- summary(dbrda)$cont$importance[2, 1:2] * 100

# Plot
ggplot(scores.dbrda, aes(x = CAP1, y = CAP2, color = landform)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=landformCol) +
  labs(x = paste0("RDA1 (", round(eig.dbrda[1], 1), "%)"),
       y = paste0("RDA2 (", round(eig.dbrda[2], 1), "%)")) +
  theme_classic() +
  theme(axis.ticks = element_blank())

############################

## Conclusion: standard RDA better explanatory power
# Export full results for SI table
anova(rda, by="terms") %>%                       
  broom::tidy() %>%                                     
  filter(term !="Residual") %>%                        
  mutate(signif = case_when(                            
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "")) %>%
  drop_na(p.value) %>%
  arrange(p.value) %>%                                 
  write.table("res_rda.txt", sep="\t", row.names=F, quote=F)

####################################

# VARIANCE PARTITIONING

# Remove zero-variance taxa
vars <- apply(counts.sub, 2, var, na.rm=T)
counts.clean <- counts.sub[, vars > 0]

# Prepare Depth variable
depthVar <- env.mat %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, depth) %>%
  filter(sample %in% rownames(counts.clean)) %>%
  column_to_rownames("sample")

# Prepare Geomorphology variables (remove depth, lat, lon)
geomorphVars <- geoVars %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::select(-depth, -lat, -lon) %>%
  filter(sample %in% rownames(counts.clean)) %>%
  column_to_rownames("sample")

# Remove zero-variance columns in geomorphVars
geomorphVars <- geomorphVars[, apply(geomorphVars, 2, var, na.rm=T) > 0]

# Prepare Landform variables as dummy matrix
landformVars <- meta %>%
  rownames_to_column("sample") %>%
  filter(sample %in% rownames(counts.clean)) %>%
  column_to_rownames("sample") %>%
  { model.matrix(~ landform - 1, data = .) }

# Remove zero-variance columns in landformVars
landformVars <- landformVars[, apply(landformVars, 2, var, na.rm=T) > 0]

# Prepare Spatial variables (lat, lon)
spatialVars <- env.mat %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, lat, lon) %>%
  filter(sample %in% rownames(counts.clean)) %>%
  column_to_rownames("sample")

# PCA on Geomorphology, Landform, and Spatial to reduce collinearity
pca.geo <- prcomp(geomorphVars, scale.=T)$x[, 1:2]
pca.lf <- prcomp(landformVars, scale.=T)$x[, 1:2]
pca.spatial <- prcomp(spatialVars, scale.=T)$x[, 1:2]

# Run variance partitioning
varpart <- varpart(
  counts.clean,
  depthVar,
  pca.geo,
  pca.lf,
  pca.spatial)

# Plot -- export size 5 x 8
plot(
  varpart,
  bg = c("darkviolet","darkorange","mediumseagreen","skyblue"),
  Xnames = c("Depth","Geomorphology","Landform","Space"))


###############################################################
 ##  Which taxa differ between landforms?
###############################################################

# Kruskal-Wallis test
kruskal.tax <- counts.image %>%
  group_by(taxon) %>%
  kruskal_test(density_m2 ~ landform) %>%
  ungroup() %>%
  filter(p < 0.05)

# Dunn's test
dunn.tax <- counts.image %>%
  filter(taxon %in% kruskal.tax$taxon) %>%
  group_by(taxon) %>%
  dunn_test(density_m2 ~ landform, p.adjust.method="BH") %>%
  ungroup() %>%
  filter(p.adj < 0.05)

# Extract results
sig.pairs <- dunn.tax %>%
  filter(p.adj < 0.05) %>%
  dplyr::select(taxon, landformA = group1, landformB = group2)
sig.tax <- unique(c(sig.pairs$taxon))
sig.landform <- unique(c(sig.pairs$landformA, sig.pairs$landformB))

## PLOT
# Later combined with diversity & density
plot3 = counts.image %>%
  filter(taxon %in% sig.tax, landform %in% sig.landform) %>%
  group_by(landform, taxon) %>%
  summarise(value = mean(density_m2, na.rm=T), .groups="drop") %>%
  mutate(
    value10 = value*10+1e-3,  # small offset accounting for true zero counts
    taxon = factor(taxon, levels = rev(c(
      "Ophiuroids","Corals","Demosponges","Sea pens",
      "Bryozoans","Glass sponges","Asteroids","Anemones",
      "Echinoids","Crinoids","Sea cucumbers"))),
    landform = factor(landform, levels = rev(c(
      "Sand ripples","Flat area","Depression",
      "Broad slope","Ridges","Steep slope","Terrace")))) %>%
  ggplot(aes(x = landform, y = taxon, fill = value10)) +
  geom_tile() +
  #facet_grid(dive~.)+
  labs(fill = "Density per 10 m²") +
  scale_color_manual(values = c(
    "FALSE"="black","TRUE"="white"), 
    guide = "none") +
  scale_fill_gradientn(
    trans = "log10",
    colors = c("floralwhite","lemonchiffon","aliceblue","lightsteelblue2","lightsteelblue4","gray8"),
    breaks = c(1, 5, 10, 25, 70, 160),
    limits = c(1, 160),
    na.value = "white",
    oob = scales::squish) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "right",
    axis.ticks = element_blank()) 

## FIGURE X: diversity, densities and composition
# Export size 7 x 5.5
plot_grid(
  plot1, plot2, plot3,
  rel_heights = c(0.8,0.8,1.2),
  align="v",
  ncol=1)

# Sponge correlation?
cor(
  as.numeric(counts.corr["Demosponges", ]),
  as.numeric(counts.corr["Glass sponges", ]),
  method = "spearman",
  use = "complete.obs")

####################################

## SI FIgure: dunn.env + dunn.taxon details all.vs.all
# Define priority landform pair (manual order)
sitePrio <- c("Terrace","Steep slope")
siteAll <- sort(unique(c(dunn.tax$group1, dunn.tax$group2, dunn.env$group1, dunn.env$group2)))
siteLevels <- c(sitePrio, setdiff(siteAll, sitePrio))

# Process TAXA
allComp.tax <- dunn.tax %>%
  mutate(
    landformA = if_else(match(group1, siteLevels) <= match(group2, siteLevels), group1, group2),
    landformB = if_else(landformA == group1, group2, group1),
    landform_pair = paste(landformA, landformB, sep = " vs "))

# Process ENV 
allComp.env <- dunn.env %>%
  mutate(
    landformA = if_else(match(group1, siteLevels) <= match(group2, siteLevels), group1, group2),
    landformB = if_else(landformA == group1, group2, group1),
    landform_pair = paste(landformA, landformB, sep = " vs "))

# Create full ordered landform pair list from taxon data
allComp.order <- allComp.tax %>%
  distinct(landformA, landformB, landform_pair) %>%
  arrange(factor(landformA, levels = siteLevels), factor(landformB, levels = siteLevels)) %>%
  pull(landform_pair)

# Apply same factor levels to both datasets
allComp.tax <- allComp.tax %>%
  mutate(
    landform_pair_fct = factor(landform_pair, levels = allComp.order),
    sig_label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ ""))

allComp.env <- allComp.env %>%
  filter(!variable %in% c("lat","lon","area","depth")) %>%
  mutate(
    variable = str_remove(variable, "(_10cm|_20cm|_1m|_2m)"),
    landform_pair_fct = factor(landform_pair, levels = allComp.order),
    variable = factor(variable, levels = rev(c(
      "slope","roughness","aspect","profileCurvature","planCurvature"))))

# Fill in missing combinations for plotting alignment
allComp.tax <- allComp.tax %>%
  complete(landform_pair_fct, taxon) %>%
  filter(landform_pair!="NA vs NA")
allComp.env <- allComp.env %>%
  complete(landform_pair_fct, variable) %>%
  drop_na(landform_pair)  %>%
  filter(landform_pair!="NA vs NA")

# Compute global max
maxVal <- max(c(
  -log10(allComp.env$p.adj), -log10(allComp.tax$p.adj)), na.rm=T)

# plot ENV
plot4 = ggplot(allComp.env, aes(x = landform_pair_fct, y = variable)) +
  geom_point(aes(size = -log10(p.adj), color = -log10(p.adj)), alpha = 0.8) +
  scale_size_continuous(name = expression(-log[10](adj~p)), range = c(2, 7), limits = c(0, maxVal)) +
  scale_color_viridis_c(option = "magma", begin = 0, end = 0.9, direction = -1, limits = c(0, maxVal)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "gray30", linewidth = 0.2),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray87", linewidth = 0.2))

# plot taxa
plot5 = allComp.tax %>%
  filter(landform_pair_fct %in% allComp.env$landform_pair_fct) %>%
  ggplot(aes(x = landform_pair_fct, y = factor(taxon, levels = rev(c(
    "Corals","Asteroids","Demosponges","Glass sponges",
    "Ophiuroids","Sea pens","Bryozoans","Anemones","Echinoids","Crinoids"))))) +
  geom_point(aes(size = -log10(p.adj), color = -log10(p.adj)), alpha = 0.8) +
  scale_size_continuous(name = expression(-log[10](adj~p)), range = c(2, 7), limits = c(0, maxVal)) +
  scale_color_viridis_c(option="magma", begin = 0, end = 0.9, direction = -1, limits = c(0, maxVal)) +
  theme_minimal()  +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.line = element_line(color = "gray30", linewidth = 0.2),
    panel.grid.major = element_line(color = "gray87", linewidth = 0.2),
    legend.position = "right")

## FIGURE SX: all-vs-all ENV+TAX
# Export size 6 x 7 
plot_grid(
  plot4, plot5, 
  rel_heights = c(0.4,1.2),
  align="v",
  ncol=1)

########################################

## PAIRWISE TAXA_TERRAIN CORRELATIONS
# report some (no figure)

# Reshape counts.corr (taxa x images -> long) ---
cor.tax <- counts.corr %>%
  as.data.frame() %>%
  tibble::rownames_to_column("taxon") %>%
  pivot_longer(-taxon, names_to = "image", values_to = "density") %>%
  left_join(meta %>% rownames_to_column("image"),
            by = "image")


# --- Step 5: run only for best-resolution vars ---
# Example: suppose you have dunnVars defined somewhere:
# dunnVars <- c("slope_10cm","roughness_10cm","planCurvature_10cm","profileCurvature_10cm")
# And add the 25 m vars:
terrainVars <- meta %>%
  dplyr::select(matches("25m$"), all_of(dunnVars))

cor.res <- cor.tax %>%
  group_by(taxon) %>%
  summarise(across(.cols = names(terrainVars),.fns  = ~{
    if (all(is.na(.x)) | all(is.na(density))) return(NA_real_)
    tryCatch(cor(.x, density, method = "spearman", use = "pairwise.complete.obs"),
             error = function(e) NA_real_)},
    .names = "{.col}"), .groups = "drop") %>%
  pivot_longer(-taxon, names_to="variable", values_to="spearman_rho")

########################################
# remove temp-data
rm()
