rm(list=ls())

# Compute functional space axes using mFD package

library(tidyverse)
library(mFD)

#### Load data ####

# load data with complete functional traits for all your species
fishes <- readRDS("Data/Fish_exo_database.rds")


#### Obtain clean trait data ####

# One column with species name
# one column for each trait

##detritus is aliased with the other diet variables
fishes_ft <- fishes %>%
  # select all traits of interest
  dplyr::select(
    c(Species, EhBd, BlBd, PFiBd, MoBd, Length, detritus, nekton,
      plants, zoobenthos, zooplankton, RepGuild1, Amplitudetemp)) %>% 
  # remove potential duplicates
  distinct() %>%
  # categorical variables as factor
  mutate_if(is.character, as.factor) %>%
  # species to rownames
  column_to_rownames("Species")

colSums(is.na(fishes_ft))
str(fishes_ft)
summary(fishes_ft)


###### Clean traits ######

fishes_ft_trans <- fishes_ft

# order parental care
fishes_ft_trans$RepGuild1 <- factor(fishes_ft_trans$RepGuild1, 
                                    levels = c("nonguarders", "guarders","bearers"),
                                    ordered = T)

str(fishes_ft_trans)
levels(fishes_ft_trans$RepGuild1)

num_traits <- colnames(fishes_ft_trans %>% dplyr::select(where(is.numeric)))
ordered_traits <- colnames(fishes_ft_trans %>% dplyr::select(where(is.factor)))

traits_cat <- data.frame(
  trait_name = colnames(fishes_ft_trans)) %>%
  mutate(trait_type = case_when(
    trait_name %in% num_traits ~ "Q",
    trait_name %in% ordered_traits ~ "O"
  ))


##### Trait-based distances between species #########


sp_dist_fish <- mFD::funct.dist(
  sp_tr         = fishes_ft_trans,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

#### Functional space & quality ####

fspaces_quality_fish <- mFD::quality.fspaces(
  sp_dist             = sp_dist_fish,
  maxdim_pcoa         = 222,
  deviation_weighting = c("absolute", "squared"),
  fdist_scaling       = FALSE,
  fdendro             = "average")

round(fspaces_quality_fish$quality_fspaces, 3)

fspaces_quality_fish$details_fspaces$pc_eigenvalues
colSums(fspaces_quality_fish$details_fspaces$pc_eigenvalues)
fspaces_quality_fish$details_fspaces$pc_eigenvalues$Eigenvalues[1:10]/18.4432899

# contribution of each axis to variance 
# [1] 0.24733142 0.17275135 0.13757403 0.10185779 0.06980763 0.05601513
# [7] 0.04281421 0.03489633 0.02204993 0.01687940

mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_fish,
  quality_metric             = "rmsd", # or "mad"
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


#### correlation between traits & axis ####

sp_faxes_coord_fish <- fspaces_quality_fish$"details_fspaces"$"sp_pc_coord"

fish_tr_faxes_bin <- mFD::traits.faxes.cor(
  sp_tr          = fishes_ft_trans %>% select(where(is.integer)), 
  sp_faxes_coord = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

fish_tr_faxes_other <- mFD::traits.faxes.cor(
  sp_tr          = fishes_ft_trans %>% select(-where(is.integer)), 
  sp_faxes_coord = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

# plot correlation between traits and axes 
fish_tr_faxes_bin$tr_faxes_plot
fish_tr_faxes_other$tr_faxes_plot

fish_tr_faxes_other$"tr_faxes_stat"
# filter only significant corr
fish_tr_faxes_other$"tr_faxes_stat"[
  which(fish_tr_faxes_other$"tr_faxes_stat"$"p.value" < 0.05), ]


fish_tr_faxes_bin$"tr_faxes_stat"
# filter only significant corr
fish_tr_faxes_bin$"tr_faxes_stat"[
  which(fish_tr_faxes_bin$"tr_faxes_stat"$"p.value" < 0.05 &
          fish_tr_faxes_bin$"tr_faxes_stat"$"value" > 0.2), ]



stat_traits_axes_14 <- bind_rows(fish_tr_faxes_other$"tr_faxes_stat", 
                                 fish_tr_faxes_bin$"tr_faxes_stat")




#### PC1 to PC4 axes functional space ####

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  plot_ch         = TRUE,
  plot_vertices   = TRUE,
  plot_sp_nm      = NULL,
  check_input     = TRUE)

#big_plot$"patchwork"


coord_4d <- as.data.frame(sp_faxes_coord_fish[ , c("PC1", "PC2", "PC3", "PC4")])


saveRDS(coord_4d, "Output/01_PCoA_coordinates")

correlation <- cor(fishes_ft_trans %>% 
                     mutate(RepGuild1 = as.numeric(RepGuild1)), 
                   coord_4d)


saveRDS(correlation, "Output/01_PCoA_corr_traits")


cor_long <- as.data.frame(correlation) %>%
  rownames_to_column("trait") %>%
  pivot_longer(cols = PC1:PC4, names_to = "axis", values_to = "corr")

stat_traits_axes_14_all <- left_join(cor_long, stat_traits_axes_14,
                                     by = c("trait","axis")) %>%
  select(-c(test, stat)) %>%
  pivot_wider(
    names_from = axis,
    names_glue = "{axis}_{.value}",
    names_sort = T,
    values_from = c(corr, value, p.value)
  ) %>%
  mutate_if(is.numeric, round, 3)

colnames(stat_traits_axes_14_all)

write.csv2(stat_traits_axes_14_all, "Output/01_table_corr_traits_axis.csv")


