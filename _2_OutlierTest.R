
# Prepare data
outlier.test = meta %>%
dplyr::select(c(matches("10cm|20cm|1m|2m")))

# Define resolutions and name patterns 
resolutions <- c("10cm","20cm","1m","2m")

# Function to select columns by name pattern
outlier_detect <- function(df, res) {
  # Adjust this pattern depending on your exact column names
  cols <- grep(res, colnames(df), value=T)
  df[, cols, drop = F]}

# Run function
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

# Combine; count how many times each sample is flagged as outlier
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

# Remove temp-data
rm(resolutions, res, md_res, scores_res, pca_res, geo_res, outliers_res)
