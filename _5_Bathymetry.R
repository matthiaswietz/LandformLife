
#################################################################
 ## LINKS BETWEEN BATHY RESOLUTIONS 
#################################################################

## OFOBS vs SHIP BATHY -- pairwise correlations 

# Order variables
varOrder <- c("slope","roughness","aspect","profileCurvature","planCurvature")
resOrder <- c("2m","1m","20cm","10cm","25m")

# Helper function to order factors
order_var <- function(x) {
  base <- str_extract(x, paste(varOrder, collapse = "|"))
  res  <- str_extract(x, paste(resOrder, collapse = "|"))
  uniq_levels <- unique(x[order(factor(base, levels = varOrder), factor(res, levels = resOrder))])
  factor(x, levels = uniq_levels)}

## Plot blue and symbols
resCorr = varCorr %>%
  filter(xor(str_detect(var1,"_25m$"), str_detect(var2,"_25m$"))) %>%
  filter(abs(correlation) > 0.11) %>%
  mutate(
    hiRes = if_else(str_detect(var1,"_25m$"), var2, var1),
    loRes = if_else(str_detect(var1,"_25m$"), var1, var2),
    shape = factor(if_else(correlation >= 0, "pos", "neg"), levels = c("pos","neg")),
    hiRes_ord = order_var(hiRes),
    loRes_ord = order_var(loRes)) %>%
  arrange(desc(abs(correlation)))

# Export size 5x6
ggplot(resCorr, aes(x = loRes_ord, y = hiRes_ord)) +
  geom_point(aes(
    shape = shape,
    color = abs(correlation),
    fill = abs(correlation),   # Optional: used to fill triangle
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
  theme_classic(base_size = 11) +
  theme(
    panel.grid.major = element_line(linewidth = 0.001, color = "gray92"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank())

###########################################################

## MODELLING
# Linking high- and low-res bathymetry?

## Due to collinearity of slope & roughness:
## check whether one of them adds more variance
## If yes: only select this one
## Done for coral as most abundant taxon
varSelect <- counts.corr %>%
  rownames_to_column("taxa") %>%
  filter(taxa=="Corals") %>%
  column_to_rownames("taxa") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("image") %>%
  left_join(
    meta %>%
      rownames_to_column("image") %>%
      select(image, all_of(dunnVars)),
    by="image") 

# Tested by ANOVA
lm1 <- lm(Corals ~ scale(roughness_2m), data = varSelect)
lm2 <- lm(Corals ~ scale(slope_1m), data = varSelect)
lm3 <- lm(Corals ~ scale(roughness_2m) + scale(slope_1m), data=varSelect)
anova(lm1, lm3)
avPlots(lm3)
## Significant improvement with slope over roughness 
## Hence, only slope was included
## Also in line with model.geo: only slope significant

####################################

## Formula for extracting stats

mvabund_stats <- function(model, model_res, model_name="Model") {
  coefs <- coef(model)
  uni_p <- model_res$uni.p
  taxa <- colnames(uni_p)
  predictors <- rownames(uni_p)
  purrr::map_dfr(taxa, function(taxon) {
    tibble(
      Taxon = stringr::str_replace_all(taxon, "\\.", " "), # remove everything after "_"
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

####################################

## Separate models, then overlapped
# First: high-res OFOBS bathy (scaled predictors)
predict.hires <- meta %>%
  select(all_of(setdiff(dunnVars, c("lat","lon","depth")))) %>%
  drop_na() %>%
  rename_with(~ str_remove(.x, "_.*")) %>%
  scale() %>% 
  as.data.frame() %>%
  dplyr::select(-c("roughness"))

# round counts, match everything
#counts.sdm <- round(t(counts.corr)[rownames(predictors), , drop=F])
#predictors <- predictors[rownames(counts.sdm), , drop=F]

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
 dplyr::select(contains("_25m")) %>%
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
  add_column(data="Shipboard")

##########################

# Plot effect size comparison
rbind(stats.hires, stats.lowres) %>%
  filter(Significant == TRUE, abs(EffectSize) > 0.01) %>%
  inner_join(
    rbind(stats.hires, stats.lowres) %>%
      filter(Significant == TRUE, data == "OFOBS") %>%
      dplyr::select(Taxon, Predictor, SignOFOBS = EffectSize),
    by = c("Taxon", "Predictor")) %>%
  mutate(
    Predictor=factor(Predictor, levels=c("slope","aspect")),
    Taxon=factor(Taxon, levels=rev(c(
      "Corals","Demosponges","Glass sponges","Ophiuroids","Crinoids","Sea pens")))) %>%
  filter(sign(EffectSize) == sign(SignOFOBS)) %>%
  ggplot(aes(
    x = EffectSize, y = Taxon,
    color = data, shape = data)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  facet_wrap(~ Predictor, scales = "free_x") +
  scale_x_continuous(
    expand=c(0.1, 0.1), n.breaks = 3.58) +
  scale_color_manual(values = c(
    "OFOBS" = "deepskyblue4",
    "Shipboard" ="lightcoral"),
    name="Bathymetry") +
  scale_shape_manual(values = c(
    "OFOBS" = 16,"Shipboard" = 17),
    name="Bathymetry") +
  theme_bw() + theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = "black"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2))

# Plot correlation between effect sizes
bind_rows(stats.hires, stats.lowres) %>%
  filter(Significant == TRUE, abs(EffectSize) > 0.01) %>%
  filter(Predictor == "slope") %>%
  select(Taxon, EffectSize, data) %>%
  pivot_wider(names_from = data, values_from = EffectSize) %>%
  drop_na() %>%
  mutate(diff = abs(OFOBS - Shipboard)) %>% { 
    top_labels <- slice_max(., order_by = diff, n = 6)$Taxon
    ggplot(., aes(x = OFOBS, y = Shipboard)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method="lm", se=F, color="gray28") +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="grey40") +
      geom_text_repel(
        data = filter(., Taxon %in% top_labels),
        aes(label = Taxon),
        size = 4, box.padding = 0.35, point.padding = 0.3,
        segment.color="grey50") +
  theme_classic() + theme(
    panel.grid.major = element_line(linewidth = 0.2))}

## Report in text:
# Number of significance; avg effect size
rbind(stats.hires, stats.lowres) %>%
  group_by(data, Predictor) %>% summarize(
    mean_abs_effect = mean(abs(EffectSize[Significant]), na.rm=T),
    n_significant = sum(Significant), .groups="drop") %>%
  arrange(desc(mean_abs_effect))

## Report in text:
# Spearman correlation
cor(
  meta$slope_1m, meta$slope_25m, 
  method="spearman")


#################################################################
 ## UPSCALING -- CORAL vs SLOPE via GAM
#################################################################

# Load bathy data
ps118_bathy <- rast("../env/PS118_25m_bathy.img")
ps118_slope <- rast("../env/PS118_25m_slope.tif")

# Set coordinate system
crs <- crs(ps118_slope) 

# Combine points and slope
ps118_combined <- as.data.frame(ps118_slope, xy=T, na.rm=T)

# Convert to SpatVector
pts <- vect(ps118_combined, geom = c("x","y"), crs=crs)

# Extract values from bathy based on slope locations
bathyVals <- terra::extract(ps118_bathy, pts)

# Combine slope + bathy
ps118_combined$depth <- bathyVals[,2]

# Scale for GAM 
slope.pred <- meta %>%
  dplyr::select(c("slope_1m")) %>%
  rename_with(~ str_remove(.x, "_.*")) %>%
  drop_na() %>%
  scale() %>% 
  as.data.frame()

####################################################

# Prepare & round counts; join slope
counts.coral <- slope.pred %>%
rownames_to_column("image") %>%
left_join(t(counts.corr)[ , "Corals", drop=F] %>% 
as.data.frame() %>%
rownames_to_column("image") %>%
mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image")

counts.demo <- slope.pred %>%
rownames_to_column("image") %>%
left_join(t(counts.corr)[ , "Demosponges", drop=F] %>% 
as.data.frame() %>%
rownames_to_column("image") %>%
mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), by="image")

counts.glass <- slope.pred %>%
rownames_to_column("image") %>%
left_join(t(counts.corr)[ , "Glass sponges", drop=F] %>% as.data.frame() %>%
rownames_to_column("image") %>%
mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), 
by="image") %>% rename(GlassSponges = `Glass sponges`)

counts.crinoi <- slope.pred %>%
rownames_to_column("image") %>%
left_join(t(counts.corr)[ , "Crinoids", drop=F] %>% 
as.data.frame() %>%
rownames_to_column("image") %>%
mutate(across(-image, ~ ifelse(. > 0 & . < 1, 1, round(.)))), 
by="image")

############################

# Check model performance (for coral)
# Split into train (80%) and test (20%)
set.seed(123)
trainID <- sample(seq_len(nrow(counts.coral)), size = 0.8 * nrow(counts.coral))
trainData <- counts.coral[trainID, ]
testData  <- counts.coral[-trainID, ]

# Fit GAM on training data
gam.train <- bam(
  Corals ~ slope,     # Corals is the response
  data = trainData,
  family = nb(),
  select = TRUE)

# Predict on test data
pred <- predict(gam.train, newdata = testData, type="response")

# Evaluate
sqrt(mean((testData$Corals - pred)^2)) #rmse: 6.83
mean(abs(testData$Corals - pred)) #mae: 4.06
cor(testData$Corals, pred) #0.39

############################

# Run models
gam.coral <- gam(
  Corals ~ slope,
  data = counts.coral,
  family = nb())
gam.demo <- gam(
  Demosponges ~ slope,
  data = counts.demo,
  family = nb())
gam.glass <- gam(
  GlassSponges ~ slope,
  data = counts.glass,
  family = nb())
gam.crinoi <- gam(
  Crinoids ~ slope,
  data = counts.crinoi,
  family = nb())

##########################

# Prepare scaling (link to slope.pred)
meta %>%
  dplyr::select(c("slope_1m")) %>%
  rename_with(~ str_remove(.x, "_.*")) %>%
  drop_na() %>%
  as.data.frame() %>%
  scale() -> scaled

# Extract scaling parameters
slopeMean <- attr(scaled, "scaled:center")["slope"]
slopeSD   <- attr(scaled, "scaled:scale")["slope"]

# Apply scaling; run predictive models 
ps118_predict <- ps118_combined %>% 
  mutate(slope = (PS118_25m_slope - slopeMean) / slopeSD) %>%
  mutate(
    corals_predicted  = predict(gam.coral, newdata =., type="response"),
    demosponges_predicted = predict(gam.demo,  newdata =., type="response"),
    glassSponges_predicted = predict(gam.glass, newdata =., type="response"),
    crinoids_predicted = predict(gam.crinoi, newdata =., type="response")) %>%
  dplyr::rename(slope_scaled=slope) 

# Total number of individuals
ps118_predict %>%
  summarise(
    corals_sum = sum(corals_predicted, na.rm=T), # 61402086
    demosponges_sum = sum(demosponges_predicted, na.rm=T), #22825904
    glassSponges_sum = sum(glassSponges_predicted, na.rm=T), #11633960
    crinoids_sum = sum(crinoids_predicted, na.rm=T)) #19253688
    
####################################################

## GEOTIFF EXPORT

# Coral only
writeRaster(c(
  rast(ps118_predict[, c("x","y","PS118_25m_slope")], type="xyz", crs=crs),
  rast(ps118_predict[, c("x","y","depth")], type="xyz", crs = crs),
  rast(ps118_predict[, c("x","y","corals_predicted")], type="xyz", crs=crs)),
  "coralPredictions.tif", filetype="GTiff", overwrite=T)

# All taxa
writeRaster(c(
  rast(ps118_predict[, c("x","y","PS118_25m_slope")], type="xyz", crs=crs),
  rast(ps118_predict[, c("x","y","depth")], type="xyz", crs = crs),
  rast(ps118_predict[, c("x","y","corals_predicted")], type="xyz", crs=crs),
  rast(ps118_predict[, c("x","y","demosponges_predicted")], type="xyz", crs=crs),
  rast(ps118_predict[, c("x","y","glassSponges_predicted")], type="xyz", crs=crs),
  rast(ps118_predict[, c("x","y","crinoids_predicted")], type="xyz", crs=crs)),
  "allPredictions.tif", filetype="GTiff", overwrite=T)

###########################################################

# remove temp.data
rm(Y, X)
