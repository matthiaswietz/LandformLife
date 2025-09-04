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
library(tidync)
library(sf)
library(FNN)
library(dbscan)
library(caret)
library(car)
library(broom)
library(iNEXT)

setwd("/Users/mwietz/ownCloud - mwietz@owncloud.mpi-bremen.de/AWI_MPI/collaborations/meineMiao/LandformLife/Rstats")
setwd("Y:/AWI_MPI/AWI_MPI/collaborations/meineMiao/LandformLife/Rstats")
save.image("~/LandformLife.Rdata")
#load("~/LandformLife.Rdata")

####################################################

## LOAD DATA

# Load image areas
# Correct for partially counted images
areas <- read.table(
  "areas.txt", header=T, sep="\t", check.names=F) %>%
  mutate(area = image_m2 * rectify) %>%
  distinct(image_filename, .keep_all=T) %>%
  dplyr::select(image_filename, area, landform)

# Load bathy data
# Join area info and hydro data 
meta <- read.table(
  "bathy.txt", h=T, sep="\t", check.names=F) %>%
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

# Load counts, reformat names
counts.raw <- read.table(
  "counts.txt", header=T, sep="\t", check.names=F) %>% 
  dplyr::select(-c(
    "Polychaete","Crustaceans", "Fish",
    "Worms","Sea cucumbers","laser point")) %>%
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

# Normalize counts by corrected image areas 
area.vec <- setNames(meta$area, rownames(meta))

# Apply to raw counts 
counts.corr <- counts.raw %>%
  dplyr::select(all_of(names(area.vec))) %>%
  mutate(across(everything(), ~ .x / area.vec[cur_column()]))

# Set plotting colors
landformCol =c(
  "Sand ripples"="gray22",
  "Flat area"="yellow2",
  "Depression"="deeppink4",
  "Broad slope"="goldenrod3",
  "Ridges"="gray44",
  "Steep slope"="mediumturquoise",
  "Terrace"="turquoise4")
