# hypervolume overlap and distance

rm(list=ls())

library(BAT)
library(hypervolume)
library(tidyverse)
library(parallel)

#---------------- Load data and set up species groups ---------------

# load PC axis 
#pc_axis <- readRDS("Output/01_PCA_coordinates_10D")
pc_axis <- readRDS("Output/01_PCoA_coordinates")
colnames(pc_axis) <- paste0("V", 1:4)


# load data with species groups
fishes <- readRDS("Data/Fish_exo_database.rds")

# set up species groups
fish_groups <- fishes %>% 
  mutate(
    Outside_cat = if_else(Outside==1, "Exists outside", "Only inside"),
    Nb_intro = if_else(occ.introduced > 4, "More than 4", "1 to 4")) %>%
  mutate(combi_gp = paste(Nb_intro, Outside_cat, sep = " & ")) %>%
  dplyr::select(Species, Outside_cat, Nb_intro, combi_gp)

fish_gp_hv <- left_join(
  fish_groups, 
  pc_axis %>% rownames_to_column("Species"),
  by = "Species")

colSums(is.na(fish_gp_hv))
str(fish_gp_hv)
summary(fish_gp_hv)

#### Compute hv Using BAT package -------------

comm <- fish_groups %>% select(Species, combi_gp) %>%
  mutate(val = 1) %>%
  mutate(combi_gp = case_when(
    combi_gp == "More than 4 & Only inside" ~ "I+4",
    combi_gp == "1 to 4 & Only inside" ~ "I-4",
    combi_gp == "1 to 4 & Exists outside" ~ "O-4",
    combi_gp == "More than 4 & Exists outside" ~ "O+4"
  )) %>%
  pivot_wider(names_from = Species, values_from = val, values_fill = 0) %>%
  column_to_rownames("combi_gp")

trait <- pc_axis

(rownames(trait) == colnames(comm))


hvlist = kernel.build(comm, trait)

alph = kernel.alpha(hvlist)
beta = kernel.beta(hvlist, comp = TRUE)
#dis = kernel.dispersion(hvlist)
#even = kernel.evenness(hvlist)
sim = kernel.similarity(hvlist)
#ori = kernel.originality(hvlist)

beta
#sim



##### Fig 3 beta density & centroid distance ####

# venn diagrams between the 4 groups
library(eulerr)
library(ggplotify)
library(ggrepel)
library(gridExtra)

ggp1 <- as.ggplot(plot(euler(c(
  "In>4" = beta$Unique_to_Cols[1],
  "In<4" = beta$Unique_to_Rows[1],
  "In>4&In<4" = beta$Shared[1]))))

ggp2 <- as.ggplot(plot(euler(c(
  "In>4" = beta$Unique_to_Cols[2],
  "Out<4" = beta$Unique_to_Rows[2],
  "In>4&Out<4" = beta$Shared[2]))))

ggp3 <- as.ggplot(plot(euler(c(
  "In>4" = beta$Unique_to_Cols[3],
  "Out>4" = beta$Unique_to_Rows[3],
  "In>4&Out>4" = beta$Shared[3]))))

ggp4 <- as.ggplot(plot(euler(c(
  "In<4" = beta$Unique_to_Cols[4],
  "Out<4" = beta$Unique_to_Rows[4],
  "In<4&Out<4" = beta$Shared[4]))))

ggp5 <- as.ggplot(plot(euler(c(
  "In<4" = beta$Unique_to_Cols[5],
  "Out>4" = beta$Unique_to_Rows[5],
  "In<4&Out>4" = beta$Shared[5]))))

ggp6 <- as.ggplot(plot(euler(c(
  "Out<4" = beta$Unique_to_Cols[6],
  "Out>4" = beta$Unique_to_Rows[6],
  "Out<4&Out>4" = beta$Shared[6]))))


# distances between centroids
cent <- as.data.frame(get_centroid(hvlist)[1:4,])
rownames(cent)
rownames(cent) = c("In>4", "In<4", "Out<4", "Out>4")

cent12 <- ggplot() +
  geom_point(data = pc_axis, aes(x = V1, y = V2), alpha = 0.5)  +
  geom_label_repel(data = cent, 
                   aes(x = V1, y = V2, label = rownames(cent),
                       color = rownames(cent))) +
  geom_point(data = cent, aes(x = V1, y = V2, color = rownames(cent)),
             size = 4)+
  theme_light()+
  theme(legend.position = "none")

cent34 <- ggplot() +
  geom_point(data = pc_axis, aes(x = V3, y = V4), alpha = 0.5)  +
  geom_label_repel(data = cent, 
                   aes(x = V3, y = V4, label = rownames(cent),
                       color = rownames(cent))) +
  geom_point(data = cent, aes(x = V3, y = V4, color = rownames(cent)),
             size = 4)+
  theme_light()+
  theme(legend.position = "none")

# total plot 
#pdf("Fig/Fig3_beta_dissim.pdf", 4, 5)
(cent12 / cent34) | (ggp1 / ggp2 / ggp3 / ggp4 / ggp5 / ggp6)
#dev.off()

# Export PDF table
pdf("table_beta_sim.pdf", 4,4)       
grid.table(
  data.frame(
    B_total = c(beta$Btotal),
    B_repl = c(beta$Brepl),
    B_rich = c(beta$Brich)))
dev.off()

df <- data.frame(
  B_total = c(beta$Btotal),
  B_repl = c(beta$Brepl),
  B_rich = c(beta$Brich))



#####---------------- Null models beta / sim ------------------

repet <- 1:999

trait <- pc_axis

fish_groups_to_sample <- fish_groups %>%
  select(Species, combi_gp) %>%
  mutate(combi_gp = case_when(
    combi_gp == "More than 4 & Only inside" ~ "I+4",
    combi_gp == "1 to 4 & Only inside" ~ "I-4",
    combi_gp == "1 to 4 & Exists outside" ~ "O-4",
    combi_gp == "More than 4 & Exists outside" ~ "O+4"
  ))

library(pbapply)

nm_hv_list <- pbmcapply::pbmclapply(repet, function(x){
  
  # sample comm
  comm_null <- fish_groups_to_sample %>%
    mutate(combi_gp_null = sample(combi_gp)) %>%
    select(-combi_gp) %>%
    mutate(val = 1) %>%
    pivot_wider(names_from = Species, values_from = val, values_fill = 0) %>%
    column_to_rownames("combi_gp_null")
  
  # compute hv
  hvlist_null = kernel.build(comm_null, trait)
  
  # compute metrics
  alpha = kernel.alpha(hvlist_null)
  beta = kernel.beta(hvlist_null, comp = TRUE)
  sim = kernel.similarity(hvlist_null)
  
  # store results
  result_list <- list(
    alpha = alpha,
    beta = beta,
    sim = sim)
  
  return(result_list)},
  mc.cores = 20) # takes 50min

Sys.time()

saveRDS(nm_hv_list, "Output/03_null_mod_hv_metrics_mfd")

#### extract results from list

nm_hv_list <- readRDS("Output/03_null_mod_hv_metrics_mfd")

names(nm_hv_list[[1]][["alpha"]])
comb <- combn(names(nm_hv_list[[1]][["alpha"]]), 2)
df_2by2 <- data.frame(
  hv1 = comb[1,],
  hv2 = comb[2,])

alpha_df <- data.frame()
tot_df <- data.frame()

for(i in 1:length(nm_hv_list)){
  alpha_temp <- as.data.frame(nm_hv_list[[i]][["alpha"]])
  beta_tot_temp <- as.data.frame(as.matrix(nm_hv_list[[i]][["beta"]]$Btotal))
  beta_repl_temp <- as.data.frame(as.matrix(nm_hv_list[[i]][["beta"]]$Brepl))
  dist_cent_temp <- as.data.frame(as.matrix(nm_hv_list[[i]][["sim"]]$Distance_centroids))
  
  colnames(alpha_temp) <- "alpha"
  alpha_temp$hv = rownames(alpha_temp)
  
  df_2by2_temp <- df_2by2 %>%
    mutate(beta_tot = numeric(nrow(df_2by2)),
           beta_repl = numeric(nrow(df_2by2)),
           dist_cent = numeric(nrow(df_2by2)))
  
  for(k in 1:nrow(df_2by2_temp)){
    hv1 = df_2by2_temp$hv1[k]
    hv2 = df_2by2_temp$hv2[k]
    
    df_2by2_temp$beta_tot[k] <- beta_tot_temp[hv1, hv2]
    df_2by2_temp$beta_repl[k] <- beta_repl_temp[hv1, hv2]
    df_2by2_temp$dist_cent[k] <- dist_cent_temp[hv1, hv2]
  }
  
  # store results in df
  alpha_df <- bind_rows(alpha_df, alpha_temp)
  tot_df <- bind_rows(tot_df, df_2by2_temp)
}

##### volume of each group compared to NM volume ----------

# observed value
alph_obs <- as.data.frame(alph)
alph_obs$hv = rownames(alph_obs)
# compare observed value with null model
ggplot() +
  geom_violin(data = alpha_df, aes(x=hv, y=alpha)) + 
  geom_boxplot(data = alpha_df, aes(x=hv, y=alpha)) +
  geom_point(data = alph_obs, aes(x = hv, y = alph), color="red", size=2)

##### diff metrics compared to NM -------------

# shape null model output
tot_df_plt <- tot_df  %>%
  mutate(
    hv1_ok = case_when(
      hv1 == "I+4" ~ "In>4",
      hv1 == "I-4" ~ "In<4",
      hv1 == "O-4" ~ "Out<4",
      hv1 == "O+4" ~ "Out>4"),
    hv2_ok = case_when(
      hv2 == "I+4" ~ "In>4",
      hv2 == "O+4" ~ "Out>4",
      hv2 == "I-4" ~ "In<4",
      hv2 == "O-4" ~ "Out<4")) %>%
  mutate(combi = paste0(hv1_ok, " / ", hv2_ok)) %>%
  select(combi, beta_tot, beta_repl, dist_cent) %>%
  pivot_longer(!combi, names_to = "metric", values_to = "value")

# shape observed output
beta_repl_obs <- as.data.frame(as.matrix(beta$Brepl))
beta_tot_obs <- as.data.frame(as.matrix(beta$Btotal))
dist_cent_obs <- as.data.frame(as.matrix(sim$Distance_centroids))

df_2by2_obs <- data.frame(
  hv1 = comb[1,],
  hv2 = comb[2,],
  beta_tot = numeric(ncol(comb)),
  beta_repl = numeric(ncol(comb)),
  dist_cent = numeric(ncol(comb)))

for(k in 1:nrow(df_2by2_obs)){
  hv1 = df_2by2_obs$hv1[k]
  hv2 = df_2by2_obs$hv2[k]
  df_2by2_obs$beta_tot[k] <- beta_tot_obs[hv1, hv2]
  df_2by2_obs$beta_repl[k] <- beta_repl_obs[hv1, hv2]
  df_2by2_obs$dist_cent[k] <- dist_cent_obs[hv1, hv2]
}

# reshape table for same look as NM
df_2by2_obs_plt <- df_2by2_obs %>%
  mutate(
    hv1_ok = case_when(
      hv1 == "I+4" ~ "In>4",
      hv1 == "I-4" ~ "In<4",
      hv1 == "O-4" ~ "Out<4"),
    hv2_ok = case_when(
      hv2 == "O+4" ~ "Out>4",
      hv2 == "I-4" ~ "In<4",
      hv2 == "O-4" ~ "Out<4")) %>%
  mutate(combi = paste0(hv1_ok, " / ", hv2_ok)) %>%
  select(combi, beta_tot, beta_repl, dist_cent) %>%
  pivot_longer(!combi, names_to = "metric", values_to = "value")


ggplot() +
  geom_violin(data = tot_df_plt, aes(x=combi, y=value), fill = "grey90") + 
  geom_boxplot(data = tot_df_plt, aes(x=combi, y=value), fill = "lightblue1", alpha = .50) +
  geom_point(data = df_2by2_obs_plt, aes(x = combi, y = value), 
             color="hotpink2", size=5, shape = 18) +
  facet_wrap(~metric, scales = "fixed") +
  theme_minimal()


#only distance to centroids
dist <- ggplot() +
  geom_violin(data = tot_df_plt %>% filter(metric=="dist_cent"), 
              aes(x=combi, y=value), fill = "grey90", color=NA) + 
  geom_boxplot(data = tot_df_plt %>% filter(metric=="dist_cent"), 
               aes(x=combi, y=value), fill = "lightblue2", alpha = .50) +
  geom_point(data = df_2by2_obs_plt %>% filter(metric=="dist_cent"),
             aes(x = combi, y = value), 
             color="orange2", size=5, shape = 18) +
  xlab("Combinaison of compared groups") +
  ylab("Distance between centroids") +
  theme_bw()

pdf("Fig/Suppl_dist_centroids.pdf", 6, 4)
dist
dev.off()

# count nb of values above observed

colnames(tot_df_plt)
colnames(df_2by2_obs_plt)

obs_vs_sim <- left_join(tot_df_plt %>% rename(value_sim = value),
                        df_2by2_obs_plt %>% rename(value_obs = value),
                        by = c("combi","metric")) %>%
  mutate(diff = value_obs-value_sim) %>%
  mutate(above = if_else(diff<0, "yes","no"))

colnames(obs_vs_sim)
obs_vs_sim_count <- obs_vs_sim %>% 
  group_by(combi, metric, above) %>%
  summarise(nb = n())

saveRDS(obs_vs_sim_count, "Output/03_obs_vs_sim_values_count")

obs_vs_sim_count <- readRDS("Output/03_obs_vs_sim_values_count")
