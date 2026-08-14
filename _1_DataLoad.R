
#################################################################
 ## LOAD AND FORMAT OFOBS DATA
#################################################################

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(forcats)
library(stringr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(scales)
library(vegan)
library(rstatix)
library(FactoMineR)
library(factoextra)
library(multcompView)
library(data.table)
library(mvabund)
library(mgcv)
library(terra)
library(tidyterra)
library(tidync)
library(sf)
library(FNN)
library(dbscan)
library(caret)
library(car)
library(cluster)
library(broom)
library(iNEXT)

####################################################

## LOAD IMAGE AREAS
# Correct for partially counted images
areas <- read.table(
  "./envData/areas.txt", header=T, sep="\t", check.names=F) %>%
  mutate(area = image_m2 * rectify) %>%
  distinct(image_filename, .keep_all=T) %>%
  dplyr::select(image_filename, area, landform)

## LOAD BATHYMETRY
# Join area info and hydro data 
meta <- read.table(
  "./envData/bathymetry.txt", h=T, sep="\t", check.names=F) %>%
  left_join(areas) %>%
  mutate(image = paste0(str_replace_all(
    str_extract(image_filename, "\\d{4}_\\d{2}_\\d{2}"), "_", ""),"_",
    str_extract(image_filename, "IMG_\\d+"))) %>%
  group_by(image) %>%
  mutate(image = if (n() == 1) image else paste0(image, "_", letters[1:n()])) %>%
  ungroup() %>%
  column_to_rownames("image") %>%
  mutate(landform=as.factor(landform)) %>%
  drop_na(area) %>%
  mutate(landform=factor(landform, levels=c(
    "Sand ripples","Depression","Flat area",
    "Broad slope","Ridges","Steep slope","Terrace")))

## LOAD COUNTS
counts.raw <- read.table(
  "./envData/counts.txt", header=T, sep="\t", check.names=F) %>% 
  mutate(image = paste0(str_replace_all(
    str_extract(image_filename, "\\d{4}_\\d{2}_\\d{2}"), "_", ""),"_",
    str_extract(image_filename, "IMG_\\d+"))) %>%
  group_by(image) %>%
  mutate(image = if (n() == 1) image else paste0(image, "_", letters[1:n()])) %>%
  ungroup() %>%
  filter(image_filename %in% meta$image_filename) %>%
  column_to_rownames("image") %>%
  dplyr::select(-c(image_filename, landform)) %>%
  t() %>%
  as.data.frame()

# Normalize by corrected image areas 
areaVec <- setNames(meta$area, rownames(meta))

# Apply to raw counts (= scale to 1m2)
counts.corr <- counts.raw %>%
  dplyr::select(all_of(names(areaVec))) %>%
  mutate(across(everything(), ~ .x / areaVec[cur_column()]))


#################################################################
 ## REMOVE OUTLIERS
#################################################################

# Prepare data
outlier.test = meta %>%
  dplyr::select(c(matches("10cm|20cm|1m|2m")))

# Define resolutions and name patterns 
resolutions <- c("10cm","20cm","1m","2m")

# Function to select columns by name pattern
outlier_detect <- function(df, res) {
  cols <- grep(res, colnames(df), value=T)
  df[, cols, drop = F]}

outlier.list <- list()

for (res in resolutions) {
  cat("Processing resolution:", res, "\n")
  geo_res <- outlier_detect(outlier.test, res)
  
  # Skip if no columns for this resolution
  if (ncol(geo_res) == 0) {
    warning(paste("No columns for resolution", res))
    next
  }
  
  # Remove rows with NA to avoid PCA errors
  geo_res <- na.omit(geo_res)
  
  # Run PCA (no graph)
  pca_res <- PCA(geo_res, graph=F)
  
  # Extract individual scores (samples)
  scores_res <- as.data.frame(pca_res$ind$coord)
  
  # Calculate Mahalanobis distances on scores
  md_res <- mahalanobis(scores_res, colMeans(scores_res), cov(scores_res))
  
  # Threshold=50 (defined from exploratory plot)
  outliers_res <- rownames(scores_res)[md_res > 50]
  
  cat("Outliers at", res, "resolution:", length(outliers_res), "\n")
  
  # Store outliers
  outlier.list[[res]] <- outliers_res}

# Combine; count how often a sample is flagged as outlier
# Select outliers at 2 or more resolutions
outliers <- outlier.list %>%
  unlist() %>%
  tibble(sample = .) %>%
  count(sample) %>%
  filter(n >= 2) %>%
  pull(sample)

# Subset original data accordingly
meta <- meta[!rownames(meta) %in% outliers, , drop=F]
counts.corr <- counts.corr[, !colnames(counts.corr) %in% outliers, drop=F]

####################################################

# Set plotting colors
landformCol = c(
  "Sand ripples"="gray22",
  "Flat area"="yellow2",
  "Depression"="deeppink4",
  "Broad slope"="goldenrod2",
  "Ridges"="goldenrod4",
  "Steep slope"="mediumturquoise",
  "Terrace"="turquoise4")

hotspotCol = c(
  "Hotspot"="red",
  "Non-Hotspot"="gray")

