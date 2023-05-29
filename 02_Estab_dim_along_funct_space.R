# association between establishment dimensions and functional space axes

rm(list=ls())

library(tidyverse)
library(patchwork)


#---------------- Load data and set up species groups ---------------

# load PC axis 
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


fish_gp_c <- left_join(
  fish_groups, 
  pc_axis %>% rownames_to_column("Species"),
  by = "Species") %>% 
  mutate(OI = if_else(Outside_cat == "Only inside", "OI", "Ou"),
         Nb = if_else(Nb_intro == "1 to 4", "4-","4+"))

colSums(is.na(fish_gp_c))
str(fish_gp_c)
summary(fish_gp_c)



#---------------- Plot groups along axes ---------------



##################### biplots

cor <- readRDS("Output/01_PCoA_corr_traits")

colnames(cor) = colnames(pc_axis)

cor_keep <- as.data.frame(cor) %>%
  mutate(
    axe12 = if_else(-0.3<V1 & V1<0.3 & -0.3<V2 & V2<0.3, "no", "yes"),
    axe34 = if_else(-0.3<V3 & V3<0.3 & -0.3<V4 & V4<0.3, "no", "yes")
  ) %>%
  mutate(trait = c("EhBd", "BlBd", "PFiBd", "MoBd", "TL", "DiD", "DiN", 
                   "DiP", "DiZb", "DiZp", "PC", "TA")) %>%
  mutate_if(is.numeric, function(x){x/4})


v12 <- ggplot() +
  geom_point(data = pc_axis, aes(x = V1, y = V2), color = "gray30",
             alpha = 0.5) +
  # add arrows
  geom_segment(data = cor_keep %>% filter(axe12=="yes"), 
               aes(x = 0, y = 0, xend = 2*V1, yend = 2*V2), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "black")  +
  # add dashed arrows ends
  geom_segment(data = cor_keep %>% filter(axe12=="yes"), 
               aes(x = 0, y = 0, xend = -2*V1, yend = -2*V2), 
               lty = 5, colour = "darkgrey")+
  # add arrow labels
  geom_text(data = cor_keep %>% filter(axe12=="yes"), 
            aes(x = 2*V1, y = 2*V2, label = trait),
            size = 4, nudge_x = c(0, 0, 0, 0, -0.05),
            nudge_y = c(0.05, -0.05, -0.05, -0.05, 0.05)) +
  # axis labels - see comp_var
  labs(x = "PC1 (18.8%)", y = "PC2 (15.5%)") +
  theme_light()+
  theme(legend.position = "none")

v34 <- ggplot() +
  geom_point(data = pc_axis, aes(x = V3, y = V4), color = "gray30",
             alpha = 0.5) +
  # add arrows
  geom_segment(data = cor_keep %>% filter(axe34=="yes"), 
               aes(x = 0, y = 0, xend = 2*V3, yend = 2*V4), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  # add dashed arrows ends
  geom_segment(data = cor_keep %>% filter(axe34=="yes"), 
               aes(x = 0, y = 0, xend = -2*V3, yend = -2*V4), 
               lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text(data = cor_keep %>% filter(axe34=="yes"), 
            aes(x = 2*V3, y = 2*V4, label = trait),
            size = 4, nudge_x = c(0, 0, 0, 0, -0.05),
            nudge_y = c(0.05, -0.05, -0.05, -0.05, 0.05)) +
  # axis labels - see comp_var
  labs(x = "PC3 (12.5%)", y = "PC4 (10.0%)") +
  theme_light()+
  theme(legend.position = "none")



#################### density plots 

# on V1
v1doi <- ggplot()+
  geom_density(data = fish_gp_c, aes(x = V1, fill = OI, color = OI), alpha = 0.5)+
  scale_fill_manual(values=c("hotpink2", "hotpink4")) +
  scale_color_manual(values=c("hotpink2", "hotpink4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


v1nb <- ggplot()+
  geom_density(data = fish_gp_c, aes(x = V1, fill = Nb, color = Nb), alpha = 0.5)+
  scale_fill_manual(values=c("lightblue2", "lightblue4")) +
  scale_color_manual(values=c("lightblue2", "lightblue4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")


# on V2
v2doi <- ggplot()+
  geom_density(data = fish_gp_c, aes(x = V2, fill = OI, color = OI), alpha = 0.5)+
  scale_fill_manual(values=c("hotpink2", "hotpink4")) +
  scale_color_manual(values=c("hotpink2", "hotpink4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") + coord_flip()

v2nb <- ggplot()+
  #geom_density(data = fish_gp_c, aes(x = V2, fill = OI, color = OI), alpha = 0.5)+
  geom_density(data = fish_gp_c, aes(x = V2, fill = Nb, color = Nb), alpha = 0.5)+
  scale_fill_manual(values=c("lightblue2", "lightblue4")) +
  scale_color_manual(values=c("lightblue2", "lightblue4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") + coord_flip()

# on V3
v3doi <- ggplot()+
  geom_density(data = fish_gp_c, aes(x = V3, fill = OI, color = OI), alpha = 0.5)+
  scale_fill_manual(values=c("hotpink2", "hotpink4")) +
  scale_color_manual(values=c("hotpink2", "hotpink4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

v3nb <- ggplot()+
  #geom_density(data = fish_gp_c, aes(x = V3, fill = OI, color = OI), alpha = 0.5)+
  geom_density(data = fish_gp_c, aes(x = V3, fill = Nb, color = Nb), alpha = 0.5)+
  scale_fill_manual(values=c("lightblue2", "lightblue4")) +
  scale_color_manual(values=c("lightblue2", "lightblue4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

# on V4
v4doi <- ggplot()+
  geom_density(data = fish_gp_c, aes(x = V4, fill = OI, color = OI), alpha = 0.5)+
  scale_fill_manual(values=c("hotpink2", "hotpink4")) +
  scale_color_manual(values=c("hotpink2", "hotpink4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") + coord_flip()

v4nb <- ggplot()+
  #geom_density(data = fish_gp_c, aes(x = V4, fill = OI, color = OI), alpha = 0.5)+
  geom_density(data = fish_gp_c, aes(x = V4, fill = Nb, color = Nb), alpha = 0.5)+
  scale_fill_manual(values=c("lightblue2", "lightblue4")) +
  scale_color_manual(values=c("lightblue2", "lightblue4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") + coord_flip()



########### define the layout you want

layout2 <- "
AAA##FFF##
BBB##GGG##
CCCDEHHHIJ
CCCDEHHHIJ
CCCDEHHHIJ
"

########### final plot 

tot = v1doi + v1nb + v12 + v2nb + v2doi +
  v3doi + v3nb + v34 + v4nb + v4doi +
  plot_layout(design = layout2)

#pdf("Fig/Fig2_PCA_groups_axis.pdf", 10, 5)
tot
#dev.off()

# legend
leg_doi <- ggplot()+
  geom_density(data = fish_gp_c, aes(x = V1, fill = OI, color = OI), alpha = 0.5)+
  scale_fill_manual(values=c("hotpink2", "hotpink4")) +
  scale_color_manual(values=c("hotpink2", "hotpink4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
leg_doi

leg_nb <- ggplot()+
  geom_density(data = fish_gp_c, aes(x = V1, fill = Nb, color = Nb), alpha = 0.5)+
  scale_fill_manual(values=c("lightblue2", "lightblue4")) +
  scale_color_manual(values=c("lightblue2", "lightblue4")) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
leg_nb


#---------------- Statistical test distrib group axes ---------------


# test inside v. outside dimension

a1 <- ks.test(fish_gp_c %>% filter(OI=="OI") %>% pull(V1),
              fish_gp_c %>% filter(OI=="Ou") %>% pull(V1))

a2 <- ks.test(fish_gp_c %>% filter(OI=="OI") %>% pull(V2),
              fish_gp_c %>% filter(OI=="Ou") %>% pull(V2))

a3 <- ks.test(fish_gp_c %>% filter(OI=="OI") %>% pull(V3),
              fish_gp_c %>% filter(OI=="Ou") %>% pull(V3))

a4 <- ks.test(fish_gp_c %>% filter(OI=="OI") %>% pull(V4),
              fish_gp_c %>% filter(OI=="Ou") %>% pull(V4))

output_tbl_location <- data.frame(
  Dimension = "Location of establishment",
  Axis = c("PC1", "PC2","PC3","PC4"),
  D = round(c(a1$statistic, a2$statistic, a3$statistic, a4$statistic),3),
  p = round(c(a1$p.value, a2$p.value, a3$p.value, a4$p.value), 3)
)


# test -4 vs + 4 dimansion

b1 <- ks.test(fish_gp_c %>% filter(Nb=="4+") %>% pull(V1),
              fish_gp_c %>% filter(Nb=="4-") %>% pull(V1))

b2 <- ks.test(fish_gp_c %>% filter(Nb=="4+") %>% pull(V2),
              fish_gp_c %>% filter(Nb=="4-") %>% pull(V2))

b3 <- ks.test(fish_gp_c %>% filter(Nb=="4+") %>% pull(V3),
              fish_gp_c %>% filter(Nb=="4-") %>% pull(V3))

b4 <- ks.test(fish_gp_c %>% filter(Nb=="4+") %>% pull(V4),
              fish_gp_c %>% filter(Nb=="4-") %>% pull(V4))

output_tbl_nb_basins <- data.frame(
  Dimension = "Number of basins", 
  Axis = c("PC1", "PC2","PC3","PC4"),
  D = round(c(b1$statistic,b2$statistic, b3$statistic, b4$statistic),3),
  p = round(c(b1$p.value, b2$p.value, b3$p.value, b4$p.value), 3)
)


final_tab <- bind_rows(output_tbl_nb_basins, output_tbl_location)

write.csv2(final_tab, "Output/Suppl_table_KS_tests.csv", row.names = F)



