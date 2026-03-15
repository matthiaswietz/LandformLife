
#################################################################
 ## LINKS BETWEEN BATHY RESOLUTIONS 
#################################################################

## PAIRWISE CORRELATIONS 
# Define order
varOrder <- c("slope","roughness","aspect","profileCurvature","planCurvature")
resOrder <- c("2m","1m","20cm","10cm","25m")

# Set helper function
ord_levels <- function(x) {
  base <- stringr::str_extract(x, paste(varOrder, collapse = "|"))
  res  <- stringr::str_extract(x, paste(resOrder, collapse = "|"))
 idx <- base::order(
    factor(base, levels = varOrder),
    factor(res,  levels = resOrder))
    lev <- base::unique(x[idx])
  factor(x, levels = lev)}

resCorr <- meta %>%
  select(matches("^(aspect|planCurvature|profileCurvature|roughness|slope)_(10cm|20cm|1m|2m|25m)$")) %>%
  drop_na() %>%
  cor(method="spearman") %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to="var2", values_to="correlation") %>%
  filter(var1 < var2) %>%
  filter(
    sub(".*_(10cm|20cm|1m|2m|25m)$", "\\1", var1) !=
    sub(".*_(10cm|20cm|1m|2m|25m)$", "\\1", var2)) %>%
 arrange(desc(abs(correlation))) %>%
  filter(xor(stringr::str_detect(var1,"_25m$"), stringr::str_detect(var2,"_25m$"))) %>%
  filter(abs(correlation) > 0.10) %>%
  mutate(
    hiRes = dplyr::if_else(stringr::str_detect(var1,"_25m$"), var2, var1),
    loRes = dplyr::if_else(stringr::str_detect(var1,"_25m$"), var1, var2),
    shape = factor(dplyr::if_else(correlation >= 0,"pos","neg"), levels = c("pos","neg")),
    hiRes_ord = ord_levels(hiRes),
    loRes_ord = ord_levels(loRes)) %>%
  arrange(desc(abs(correlation)))

# Export size 5x6
# omit roughness: colinearity with slope
resCorr %>% 
  filter(!grepl("roughness", var1)) %>%
  filter(!grepl("roughness", var2)) %>%
ggplot(aes(x = loRes_ord, y = hiRes_ord)) +
  geom_point(aes(
    shape = shape,
    color = abs(correlation),
    fill = abs(correlation),  
    size = abs(correlation)), 
  color="gray22", stroke = 0.4) +
  scale_fill_gradientn(
    colors = c("aliceblue","skyblue2","dodgerblue4","midnightblue","gray8"),
    limits = c(0, 0.8)) +
  scale_shape_manual(
    values = c("pos" = 24, "neg" = 25),
    labels = c("Positive", "Negative")) +
  scale_size(range = c(1.5, 5)) +
  guides(size="none") +
  theme_classic(base_size=11) +
  theme(
    panel.grid.major = element_line(linewidth=0.001, color="gray92"),
    axis.text.x = element_text(angle=45, hjust=1),
    axis.text = element_text(size=10),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())


#################################################################
 ## DENSITY MODELLING ACROSS BATHY RESOLUTIONS 
#################################################################

## Formula for extracting stats
mvabund_stats <- function(model, model_res, model_name="Model") {
  coefs <- coef(model)
  uni_p <- model_res$uni.p
  taxa <- colnames(uni_p)
  predictors <- rownames(uni_p)
  purrr::map_dfr(taxa, function(taxon) {
    tibble(
      Taxon = stringr::str_replace_all(taxon, "\\.", " "), 
      Predictor = predictors,
      EffectSize = coefs[, taxon],
      Pvalue = uni_p[, taxon],
      Model = model_name)}) %>%
    mutate(
      Significance = case_when(
        Pvalue < 0.001 ~ "***",
        Pvalue < 0.01  ~ "**",
        Pvalue < 0.05  ~ "*",
        Pvalue < 0.1   ~ ".",
        TRUE ~ ""),
      Significant = Pvalue < 0.05)}

# First: high-res OFOBS bathy (scaled predictors)
predict.hires <- meta %>%
  select(depth, lat, lon, tidyselect::ends_with("_20cm")) %>%
  drop_na() %>%
  rename_with(~ str_remove(.x, "_.*")) %>%
  scale() %>% 
  as.data.frame() %>%
  dplyr::select(-c("roughness"))

# Prepare & round counts
counts.glm <- t(counts.corr)[rownames(predict.hires), , drop=F] |>
  as.data.frame() |>
  tibble::rownames_to_column("row") |>
  mutate(across(-row, ~ ifelse(. > 0 & . < 1, 1, round(.)))) |>
  tibble::column_to_rownames("row")

# Align predictors to counts
predict.hires <- predict.hires[rownames(counts.glm), , drop=F]

# Set data
Y <- mvabund(counts.glm) 
X <- predict.hires      

# Run model (log-link by default)
mvglm.hires <- manyglm(
  Y ~ ., data = X, 
  family="negative.binomial")

# Assess model
summary(mvglm.hires)
mvglm.hires.out <- anova(
  mvglm.hires, p.uni="adjusted", resamp="pit.trap")

# Extract output
stats.hires <- mvabund_stats(
  model = mvglm.hires,
  model_res = mvglm.hires.out) %>%
  add_column(data="OFOBS")

##########################

# Ship-bathy only
predict.lowres <- meta %>%
  select(depth, lat, lon, tidyselect::ends_with("_25m")) %>%
  drop_na() %>%
  rename_with(~ str_remove(.x, "_.*")) %>%
  scale() %>% 
  as.data.frame() %>%
  dplyr::select(-c("roughness"))

# Prepare & round counts
counts.glm <- t(counts.corr)[rownames(predict.lowres), , drop=F] |>
  as.data.frame() |>
  tibble::rownames_to_column("row") |>
  mutate(across(-row, ~ ifelse(. > 0 & . < 1, 1, round(.)))) |>
  tibble::column_to_rownames("row")

# Align predictors to counts
predict.lowres <- predict.lowres[rownames(counts.glm), , drop=F]

# Set data
Y <- mvabund(counts.glm) 
X <- predict.lowres      

# Run model (log-link by default)
mvglm.lowres <- manyglm(
  Y ~ ., data = X, 
  family="negative.binomial")

# Assess model
summary(mvglm.lowres)
mvglm.lowres.out <- anova(
  mvglm.lowres, p.uni="adjusted", resamp="pit.trap")

stats.lowres <- mvabund_stats(
  model = mvglm.lowres,
  model_res = mvglm.lowres.out)  %>%
  add_column(data="Ship")

##########################

# PLOT: effect sizes for slope
# Keep only Taxa with sigf at both resolutions 
rbind(stats.hires, stats.lowres) %>%
  filter(Significant==TRUE & Predictor=="slope") %>%
  inner_join(
    rbind(stats.hires, stats.lowres) %>%
      filter(Significant==TRUE, data=="OFOBS") %>%
      select(Taxon, Predictor, SignOFOBS = EffectSize),
    by = c("Taxon","Predictor")) %>%
  group_by(Taxon, Predictor) %>%
  filter(n() > 1) %>% 
  ungroup() %>%
mutate(
  Taxon = factor(Taxon, levels = rev(c(
  "Corals","Demosponges","Glass sponges","Ophiuroids","Sea pens")))) %>%
  filter(sign(EffectSize) == sign(SignOFOBS)) %>%
  ggplot(aes(
    x = EffectSize, y = Taxon,
    color = data, shape = data)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(
    expand=c(0.1, 0.1), n.breaks = 3.58) +
  scale_color_manual(values = c(
    "OFOBS" = "deepskyblue4",
    "Ship" ="lightcoral"),
    name="Bathymetry") +
  scale_shape_manual(values = c(
    "OFOBS"=16,"Ship"=17), name="Bathymetry") +
  theme_bw() + theme(
    strip.text = element_text(face="bold"),
    strip.background = element_rect(fill=NA, color="black"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(linewidth=0.2))

## Report in text:
# Number of significance; avg effect size
rbind(stats.hires, stats.lowres) %>%
  filter(Significant == "TRUE" & Predictor=="slope") %>%
  select(Taxon, Predictor, data, EffectSize) %>%
  tidyr::pivot_wider(names_from = data, values_from = EffectSize) %>%
  drop_na() %>%
  mutate(
    diff = Ship - OFOBS,
    abs_diff = abs(diff)) %>%
  arrange(abs_diff)


#################################################################
 ## UPSCALING -- DENSITIES vs SLOPE+DEPTH
#################################################################

# Load bathy data
ps118_bathy <- rast("ps118_bathy_25m.tif")
ps118_slope <- rast("ps118_slope_25m.tif")

# Set coordinate system
crs <- crs(ps118_slope) 

# Combine bathy, depth, slope
ps118_combined <- as.data.frame(ps118_slope, xy=T, na.rm=T) %>%
  rename(ps118_slope = ps118_slope_25m) %>%
  mutate(depth = terra::extract(
  ps118_bathy, vect(., geom = c("x","y"), crs=crs))[, 2])

####################################################

# Scale for GAM 
predGam <- meta %>%
  select(c("slope_25m","depth","lat","lon")) %>%
  rename_with(~ str_remove(.x, "_.*")) %>%
  drop_na() %>%
  scale() %>% 
  as.data.frame()

# Prepare & round counts; join slope
# Calculate taxa with slope and/or depth significance OFOBS+ship
counts.coral <- predGam %>%
  rownames_to_column("image") %>%
  left_join(t(counts.corr)[ , "Corals", drop=F] %>% 
  as.data.frame() %>%
  rownames_to_column("image") %>%
  mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image")
counts.demo <- predGam %>%
  rownames_to_column("image") %>%
  left_join(t(counts.corr)[ , "Demosponges", drop=F] %>% 
  as.data.frame() %>%
  rownames_to_column("image") %>%
  mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image")
counts.glass <- predGam %>%
  rownames_to_column("image") %>% 
  left_join(t(counts.corr)[,"Glass sponges", drop=F] %>% 
  as.data.frame() %>% rownames_to_column("image") %>%
  mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image") %>% 
  rename(GlassSponges = `Glass sponges`)
counts.ophi <- predGam %>%
  rownames_to_column("image") %>%
  left_join(t(counts.corr)[ , "Ophiuroids", drop=F] %>% 
  as.data.frame() %>%
  rownames_to_column("image") %>%
  mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image")
counts.pens <- predGam %>%
  rownames_to_column("image") %>%
  left_join(t(counts.corr)[ , "Sea pens", drop=F] %>% 
  as.data.frame() %>%
  rownames_to_column("image") %>%
  mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image") %>% 
  rename(SeaPens = `Sea pens`)
counts.crinoi <- predGam %>%
  rownames_to_column("image") %>% 
  left_join(t(counts.corr)[ , "Crinoids", drop=F] %>% 
  as.data.frame() %>%
  rownames_to_column("image") %>%
  mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image")
counts.bryo <- predGam %>%
  rownames_to_column("image") %>%
  left_join(t(counts.corr)[ , "Bryozoans", drop=F] %>% 
  as.data.frame() %>%
  rownames_to_column("image") %>%
  mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image")

# Run models
gam.coral <- gam(
  Corals ~ slope + depth,
  data = counts.coral,
  family = nb())
gam.demo <- gam(
  Demosponges ~ slope + depth,
  data = counts.demo,
  family = nb())
gam.glass <- gam(
  GlassSponges ~ slope + depth,
  data = counts.glass,
  family = nb())
gam.ophi <- gam(
  Ophiuroids ~ slope + depth,
  data = counts.ophi,
  family = nb())
gam.pens <- gam(
  SeaPens ~ slope + depth,
  data = counts.pens,
  family = nb())
gam.crinoi <- gam(
  Crinoids ~ slope + depth,
  data = counts.crinoi,
  family = nb())
gam.bryo <- gam(
  Bryozoans ~ slope + depth,
  data = counts.bryo,
  family = nb())

####################################################

# Prepare scaling (link to predGam)
scaled <- meta %>%
  dplyr::select(c("slope_25m","depth")) %>%
  rename_with(~ str_remove(.x, "_.*")) %>%
  drop_na() %>%
  as.data.frame() %>%
  scale() 

# Extract scaling parameters
slopeMean <- attr(scaled, "scaled:center")["slope"]
slopeSD   <- attr(scaled, "scaled:scale")["slope"]
depthMean <- attr(scaled, "scaled:center")["depth"]
depthSD   <- attr(scaled, "scaled:scale")["depth"]

# Apply scaling; run predictive models 
# Also depth-responsive taxa (crinoids+bryozoans) were predicted
# Full pipeline showed similar hotspot extents, but weaker ENV links  
# Hence: focus on slope-responsive taxa
# All calculations documented here for reproducibility
bioPredict <- ps118_combined %>% 
  mutate(
    depthReal = abs(depth),
    slope = (ps118_slope - slopeMean) / slopeSD,
    depth = (depthReal - depthMean) / depthSD) %>%
  mutate(
    corals = predict(gam.coral, newdata = ., type="response"),
    demosponges = predict(gam.demo,  newdata = ., type="response"),
    glassSponges = predict(gam.glass, newdata = ., type="response"),
    ophiuroids = predict(gam.ophi, newdata = ., type="response"),
    seaPens = predict(gam.pens, newdata = ., type="response"),
    crinoids = predict(gam.crinoi, newdata = ., type="response"),
    bryozoans = predict(gam.bryo, newdata = ., type="response")) %>%
  mutate(
    total_density = corals + demosponges + seaPens +  
    glassSponges + ophiuroids) %>%
  mutate(
    total_density_all = corals + demosponges +  
      glassSponges + ophiuroids + seaPens +
      crinoids + bryozoans) %>%
  dplyr::rename(
     depth_scaled = depth, slope_scaled = slope)

## Total density: 90 billion
sum(bioPredict$total_density *625 / 1e9)
sum(bioPredict$total_density_all *625 / 1e9)

# Detailed summary
bioPredict %>%
  summarise(
    total = sum(total_density, na.rm=T) *625,
    median_m2 = median(total_density, na.rm=T),
    mean_m2 = mean(total_density, na.rm=T),
    max_m2 = max(total_density, na.rm =T)) %>%
  mutate(across(everything(), ~ format(round(.x, 2), scientific=F, big.mark=","))) %>%
  pivot_longer(everything())

## GEOTIFF EXPORT
writeRaster(c(
  rast(bioPredict[, c("x","y","ps118_slope")], type="xyz", crs=crs),
  rast(bioPredict[, c("x","y","depthReal")], type="xyz", crs=crs),
  rast(bioPredict[, c("x","y","total_density")], type="xyz", crs=crs)),
  "bioPredict.tif", filetype="GTiff", overwrite=T)


#################################################################
 ## MODEL PERFORMANCE METRICS ##
#################################################################

## TRAIN-TEST SPLITS
# Check model performance (for coral)
trainData <- counts.coral[counts.coral$image %in% rownames(meta)[meta$dive %in% c("39-1","69-1")], ]
testData  <- counts.coral[counts.coral$image %in% rownames(meta)[meta$dive %in% c("81-1")], ]

# Fit GAM on training data
gam.train <- bam(
  Corals ~ slope + depth,    
  data = trainData,
  family = nb(),
  select = TRUE)
# Predict on test data
predTest <- predict(gam.train, newdata = testData, type="response")

# Evaluate
sqrt(mean((testData$Corals - predTest)^2)) #rmse: 8.2
mean(abs(testData$Corals - predTest)) #mae: 4.7
cor(testData$Corals, predTest) #0.33

############################

## STANDARD ERROR
# Convert to DT
setDT(ps118_combined)

# Scale slope/depth
ps118_combined[, `:=`(
  slope = (ps118_slope - slopeMean) / slopeSD,
  depth = (abs(depth) - depthMean) / depthSD)]

# Set function
# Scale densities (ind m⁻²) to 25 × 25 m cells
stdErr <- function(model, dt_data) {
  p <- predict(
    model,
    newdata = as.data.frame(dt_data[, .(slope, depth)]),
    type="link",
    se.fit=T)
  fit_resp <- exp(p$fit) *625
  se_resp  <- exp(p$fit) * p$se.fit *625
  res <- c(
    sum_fit = sum(fit_resp, na.rm=T),
    sum_var = sum(se_resp^2, na.rm=T))
  return(res)
}

# Run analytics
gam.comb <- list(gam.coral, gam.demo, gam.glass, gam.ophi)
gam.final <- sapply(gam.comb, stdErr, dt_data = ps118_combined)

# Global SE: 4073307
sqrt(sum(gam.final["sum_var", ])) 

############################

## DEVIANCE EXPLAINED
data.frame(
  Taxon = c("Coral","Demosponges","Glass sponges","Bryozoans","Crinoids","Ophiuroids","Sea pens"),
  Deviance_Percent = c(
    round(summary(gam.coral)$dev.expl *100, 2),
    round(summary(gam.demo)$dev.expl *100, 2),
    round(summary(gam.glass)$dev.expl *100, 2),
    round(summary(gam.bryo)$dev.expl *100, 2),
    round(summary(gam.crinoi)$dev.expl *100, 2),
    round(summary(gam.ophi)$dev.expl *100, 2),
    round(summary(gam.pens)$dev.expl *100, 2)))
