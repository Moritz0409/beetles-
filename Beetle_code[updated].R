setwd("C:/Users/PC/Documents/Beetles Case Study")

packages <- c("iNEXT", "ggpubr", "knitr", "ggplot2", "broom", "devtools", "reshape2", "vegan", "dendextend", "flextable", "officer","BiodiversityR", "dplyr", "tidyr", "BAT", "gridExtra", "spaa", "tinytex", "SpadeR")

install.packages(setdiff(packages, rownames(installed.packages()))) 
for (pkg in packages) {
  library(pkg, character.only = TRUE)
}

library("iNEXT")
library("ggpubr")
library("ggplot2")

# load data
raw_2017 <- read.csv("data_2017_species.csv")
data_2018_species <- read.csv("data_2018_species.csv")

#beetles category
species_taxonomy_hahn <-read.csv("species_taxonomy_hahn.csv")

taxonomy_2017 <- merge(x = raw_2017, y = species_taxonomy_hahn[ , c("Colkat", "xylobiont", "Rote.Liste", "Urwaldrelikt")], by.x ="FHL", by.y = "Colkat")
taxonomy_2018 <- merge(x = data_2018_species, y = species_taxonomy_hahn[ , c("Colkat", "xylobiont", "Rote.Liste", "Urwaldrelikt")], by.x ="FHL", by.y = "Colkat")
taxonomy_2018 <- taxonomy_2018 %>%
  mutate(Species = paste(Family, Genus, Art, sep = " "))

beetles_2017_xylobiont <- taxonomy_2017[which(taxonomy_2017$xylobiont == "j" & taxonomy_2017$Number != 0 ),]
beetles_2017_red <- taxonomy_2017[which(taxonomy_2017$Rote.Liste!="n"& taxonomy_2017$Number != 0),]
beetles_2017_urwaldrelikt <- taxonomy_2017[which(taxonomy_2017$Urwaldrelikt!="n"),]

### Alpha diversity
#tree species & xylobiont 2017
sp_per_trap_Rubra_xylobiont<-dcast( beetles_2017_xylobiont[beetles_2017_xylobiont$TreeSp=="Quercus rubra",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Rubra_xylobiont) <- sp_per_trap_Rubra_xylobiont[,1]
sp_per_trap_Rubra_xylobiont<-sp_per_trap_Rubra_xylobiont[,-1]
sp_per_trap_Rubra_xylobiont[sp_per_trap_Rubra_xylobiont>=1]<-1

sp_per_trap_Fraxinus_xylobiont<-dcast( beetles_2017_xylobiont[beetles_2017_xylobiont$TreeSp=="Fraxinus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Fraxinus_xylobiont) <- sp_per_trap_Fraxinus_xylobiont[,1]
sp_per_trap_Fraxinus_xylobiont<-sp_per_trap_Fraxinus_xylobiont[,-1]
sp_per_trap_Fraxinus_xylobiont[sp_per_trap_Fraxinus_xylobiont>=1]<-1

sp_per_trap_Quercus_xylobiont<-dcast( beetles_2017_xylobiont[beetles_2017_xylobiont$TreeSp=="Quercus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Quercus_xylobiont) <- sp_per_trap_Quercus_xylobiont[,1]
sp_per_trap_Quercus_xylobiont<-sp_per_trap_Quercus_xylobiont[,-1]
sp_per_trap_Quercus_xylobiont[sp_per_trap_Quercus_xylobiont>=1]<-1

sp_per_trap_Acer_xylobiont<-dcast( beetles_2017_xylobiont[beetles_2017_xylobiont$TreeSp=="Acer",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Acer_xylobiont) <- sp_per_trap_Acer_xylobiont[,1]
sp_per_trap_Acer_xylobiont<-sp_per_trap_Acer_xylobiont[,-1]
sp_per_trap_Acer_xylobiont[sp_per_trap_Acer_xylobiont>=1]<-1

sp_per_trap_Tilia_xylobiont<-dcast( beetles_2017_xylobiont[beetles_2017_xylobiont$TreeSp=="Tilia",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Tilia_xylobiont) <- sp_per_trap_Tilia_xylobiont[,1]
sp_per_trap_Tilia_xylobiont<-sp_per_trap_Tilia_xylobiont[,-1]
sp_per_trap_Tilia_xylobiont[sp_per_trap_Tilia_xylobiont>=1]<-1

sp_per_trap_Ulmus_xylobiont<-dcast( beetles_2017_xylobiont[beetles_2017_xylobiont$TreeSp=="Ulmus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Ulmus_xylobiont) <- sp_per_trap_Ulmus_xylobiont[,1]
sp_per_trap_Ulmus_xylobiont<-sp_per_trap_Ulmus_xylobiont[,-1]
sp_per_trap_Ulmus_xylobiont[sp_per_trap_Ulmus_xylobiont>=1]<-1 #turn into 0-1 values for absence/presence analysis

Forest_treesp_xylobiont= list("Fraxinus" = sp_per_trap_Fraxinus_xylobiont,"Quercus" = sp_per_trap_Quercus_xylobiont, "Tilia" = sp_per_trap_Tilia_xylobiont, "Quercus rubra" = sp_per_trap_Rubra_xylobiont, "Acer" = sp_per_trap_Acer_xylobiont, "Ulmus" = sp_per_trap_Ulmus_xylobiont)

tree_xylobiont_2 <- iNEXT(Forest_treesp_xylobiont, q=1, datatype ="incidence_raw", endpoint=2)

plot1<-ggiNEXT(tree_xylobiont_2, type=1) + ggtitle("Total number xylobiont species with rarefaction = 2 traps (2017)") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot2<-ggiNEXT(tree_xylobiont_2, type=2) + ggtitle("Coverage estimatation with rarefaction = 2 traps")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot20<-ggiNEXT(tree_xylobiont_2, type=3) + ggtitle("Species estimation based on sample coverage")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

Tree_Plot_xylo <- ggarrange(plot1, plot2, plot20, 
                            common.legend = TRUE, legend = "bottom",
                            ncol = 3, nrow = 1)
ggsave("Tree_Plot_xylo.png", width = 20, height = 10)

#extracting results from iNEXT: estimation  -> Source: chatGPT
#standardizing based on trap = 2
data_1_df <- bind_rows(tree_xylobiont_2$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))

#standardizing based on coverage = 0.6
data_1_sc <- bind_rows(data_1) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD      = approx(SC, qD, xout = 0.60)$y,
    qD.LCL  = approx(SC, qD.LCL, xout = 0.60)$y,
    qD.UCL  = approx(SC, qD.UCL, xout = 0.60)$y,
    SC      = 0.60
  ) %>%
  ungroup()

#red_list and tree species 2017
sp_per_trap_Rubra_red<-dcast( beetles_2017_red[beetles_2017_red$TreeSp=="Quercus rubra",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Rubra_red) <- sp_per_trap_Rubra_red[,1]
sp_per_trap_Rubra_red<-sp_per_trap_Rubra_red[,-1]
sp_per_trap_Rubra_red[sp_per_trap_Rubra_red>=1]<-1

sp_per_trap_Fraxinus_red<-dcast( beetles_2017_red[beetles_2017_red$TreeSp=="Fraxinus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Fraxinus_red) <- sp_per_trap_Fraxinus_red[,1]
sp_per_trap_Fraxinus_red<-sp_per_trap_Fraxinus_red[,-1]
sp_per_trap_Fraxinus_red[sp_per_trap_Fraxinus_red>=1]<-1

sp_per_trap_Quercus_red<-dcast( beetles_2017_red[beetles_2017_red$TreeSp=="Quercus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Quercus_red) <- sp_per_trap_Quercus_red[,1]
sp_per_trap_Quercus_red<-sp_per_trap_Quercus_red[,-1]
sp_per_trap_Quercus_red[sp_per_trap_Quercus_red>=1]<-1

sp_per_trap_Acer_red<-dcast( beetles_2017_red[beetles_2017_red$TreeSp=="Acer",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Acer_red) <- sp_per_trap_Acer_red[,1]
sp_per_trap_Acer_red<-sp_per_trap_Acer_red[,-1]
sp_per_trap_Acer_red[sp_per_trap_Acer_red>=1]<-1

sp_per_trap_Tilia_red<-dcast( beetles_2017_red[beetles_2017_red$TreeSp=="Tilia",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Tilia_red) <- sp_per_trap_Tilia_red[,1]
sp_per_trap_Tilia_red<-sp_per_trap_Tilia_red[,-1]
sp_per_trap_Tilia_red[sp_per_trap_Tilia_red>=1]<-1

sp_per_trap_Ulmus_red<-dcast( beetles_2017_red[beetles_2017_red$TreeSp=="Ulmus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Ulmus_red) <- sp_per_trap_Ulmus_red[,1]
sp_per_trap_Ulmus_red<-sp_per_trap_Ulmus_red[,-1]
sp_per_trap_Ulmus_red[sp_per_trap_Ulmus_red>=1]<-1

Forest_treesp_red= list("Fraxinus" = sp_per_trap_Fraxinus_red,"Quercus" = sp_per_trap_Quercus_red, "Tilia" = sp_per_trap_Tilia_red, "Quercus rubra" = sp_per_trap_Rubra_red, "Acer" = sp_per_trap_Acer_red, "Ulmus" = sp_per_trap_Ulmus_red)

tree_red_2_Shannon <- iNEXT(Forest_treesp_red, q=1, datatype ="incidence_raw", endpoint=2)

plot3<-ggiNEXT(tree_red_2_Shannon, type=1) + ggtitle("Total number red-list species with rarefaction = 2 traps (2017)")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot4<-ggiNEXT(tree_red_2_Shannon, type=2) + ggtitle("Coverage estimatation with rarefaction = 2 traps")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot21<-ggiNEXT(tree_red_2_Shannon, type=3) + ggtitle("Species estimation based on sample coverage")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))


Tree_Plot_red <- ggarrange(plot3, plot4, plot21, 
                           common.legend = TRUE, legend = "bottom",
                           ncol = 3, nrow = 1)
ggsave("Tree_Plot_red.png", width = 20, height = 10)

redlist2017_sc <- bind_rows(tree_red_2_Shannon$iNextEst) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD      = approx(SC, qD, xout = 0.50)$y,
    qD.LCL  = approx(SC, qD.LCL, xout = 0.50)$y,
    qD.UCL  = approx(SC, qD.UCL, xout = 0.50)$y,
    SC      = 0.50
  ) %>%
  ungroup()

redlist2017_trap <- bind_rows(tree_red_2_Shannon$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))

#total and tree species 2017
sp_per_trap_Rubra_total<-dcast( taxonomy_2017[taxonomy_2017$TreeSp=="Quercus rubra",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Rubra_total) <- sp_per_trap_Rubra_total[,1]
sp_per_trap_Rubra_total<-sp_per_trap_Rubra_total[,-1]
sp_per_trap_Rubra_total[sp_per_trap_Rubra_total>=1]<-1

sp_per_trap_Fraxinus_total<-dcast( taxonomy_2017[taxonomy_2017$TreeSp=="Fraxinus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Fraxinus_total) <- sp_per_trap_Fraxinus_total[,1]
sp_per_trap_Fraxinus_total<-sp_per_trap_Fraxinus_total[,-1]
sp_per_trap_Fraxinus_total[sp_per_trap_Fraxinus_total>=1]<-1

sp_per_trap_Quercus_total<-dcast( taxonomy_2017[taxonomy_2017$TreeSp=="Quercus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Quercus_total) <- sp_per_trap_Quercus_total[,1]
sp_per_trap_Quercus_total<-sp_per_trap_Quercus_total[,-1]
sp_per_trap_Quercus_total[sp_per_trap_Quercus_total>=1]<-1

sp_per_trap_Acer_total<-dcast( taxonomy_2017[taxonomy_2017$TreeSp=="Acer",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Acer_total) <- sp_per_trap_Acer_total[,1]
sp_per_trap_Acer_total<-sp_per_trap_Acer_total[,-1]
sp_per_trap_Acer_total[sp_per_trap_Acer_total>=1]<-1

sp_per_trap_Tilia_total<-dcast( taxonomy_2017[taxonomy_2017$TreeSp=="Tilia",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Tilia_total) <- sp_per_trap_Tilia_total[,1]
sp_per_trap_Tilia_total<-sp_per_trap_Tilia_total[,-1]
sp_per_trap_Tilia_total[sp_per_trap_Tilia_total>=1]<-1

sp_per_trap_Ulmus_total<-dcast( taxonomy_2017[taxonomy_2017$TreeSp=="Ulmus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Ulmus_total) <- sp_per_trap_Ulmus_total[,1]
sp_per_trap_Ulmus_total<-sp_per_trap_Ulmus_total[,-1]
sp_per_trap_Ulmus_total[sp_per_trap_Ulmus_total>=1]<-1

Forest_treesp_total= list("Fraxinus" = sp_per_trap_Fraxinus_total,"Quercus" = sp_per_trap_Quercus_total, "Tilia" = sp_per_trap_Tilia_total, "Quercus rubra" = sp_per_trap_Rubra_total, "Acer" = sp_per_trap_Acer_total, "Ulmus" = sp_per_trap_Ulmus_total)

tree_total_2 <- iNEXT(Forest_treesp_total, q=1, datatype ="incidence_raw", endpoint=2)

plot5<-ggiNEXT(tree_total_2, type=1) + ggtitle("Total species")+ ggtitle("Total number species with rarefaction = 2 traps (2017)") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot6<-ggiNEXT(tree_total_2, type=2) + ggtitle("Total species")+ ggtitle("Coverage estimatation with rarefaction = 2 traps") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot13<-ggiNEXT(tree_total_2, type=3) + geom_vline(xintercept = 0.65, linetype = "dashed") + ggtitle("Species estimation based on sample coverage") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

Tree_Plot_total <- ggarrange(plot5, plot6, plot13,  
                             common.legend = TRUE, legend = "bottom",
                             ncol = 3, nrow = 1)
ggsave("Tree_Plot_total.png", width = 20, height = 10)

total2017_sc <- bind_rows(tree_total_2$iNextEst) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD      = approx(SC, qD, xout = 0.60)$y,
    qD.LCL  = approx(SC, qD.LCL, xout = 0.60)$y,
    qD.UCL  = approx(SC, qD.UCL, xout = 0.60)$y,
    SC      = 0.60
  ) %>%
  ungroup()

total2017_trap <- bind_rows(tree_total_2$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))


#total beetles and tree species 2018
sp_per_trap_Rubra_total_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Quercus rubra",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Rubra_total_2018) <- sp_per_trap_Rubra_total_2018[,1]
sp_per_trap_Rubra_total_2018<-sp_per_trap_Rubra_total_2018[,-1]
sp_per_trap_Rubra_total_2018[sp_per_trap_Rubra_total_2018>=1]<-1

sp_per_trap_Fraxinus_total_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Fraxinus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Fraxinus_total_2018) <- sp_per_trap_Fraxinus_total_2018[,1]
sp_per_trap_Fraxinus_total_2018<-sp_per_trap_Fraxinus_total_2018[,-1]
sp_per_trap_Fraxinus_total_2018[sp_per_trap_Fraxinus_total_2018>=1]<-1

sp_per_trap_Quercus_total_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Quercus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Quercus_total_2018) <- sp_per_trap_Quercus_total_2018[,1]
sp_per_trap_Quercus_total_2018<-sp_per_trap_Quercus_total_2018[,-1]
sp_per_trap_Quercus_total_2018[sp_per_trap_Quercus_total_2018>=1]<-1

sp_per_trap_Acer_total_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Acer",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Acer_total_2018) <- sp_per_trap_Acer_total_2018[,1]
sp_per_trap_Acer_total_2018<-sp_per_trap_Acer_total_2018[,-1]
sp_per_trap_Acer_total_2018[sp_per_trap_Acer_total_2018>=1]<-1

sp_per_trap_Tilia_total_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Tilia",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Tilia_total_2018) <- sp_per_trap_Tilia_total_2018[,1]
sp_per_trap_Tilia_total_2018<-sp_per_trap_Tilia_total_2018[,-1]
sp_per_trap_Tilia_total_2018[sp_per_trap_Tilia_total_2018>=1]<-1

sp_per_trap_Ulmus_total_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Ulmus",], Species~Trap, value.var="Number")
rownames(sp_per_trap_Ulmus_total_2018) <- sp_per_trap_Ulmus_total_2018[,1]
sp_per_trap_Ulmus_total_2018<-sp_per_trap_Ulmus_total_2018[,-1]
sp_per_trap_Ulmus_total_2018[sp_per_trap_Ulmus_total_2018>=1]<-1

Forest_treesp_total_2018= list("Fraxinus" = sp_per_trap_Fraxinus_total_2018,"Quercus" = sp_per_trap_Quercus_total_2018, "Tilia" = sp_per_trap_Tilia_total_2018, "Quercus rubra" = sp_per_trap_Rubra_total_2018, "Acer" = sp_per_trap_Acer_total_2018, "Ulmus" = sp_per_trap_Ulmus_total_2018)

tree_total__2018 <- iNEXT(Forest_treesp_total_2018, q=1, datatype ="incidence_raw", endpoint=2)

plot7<-ggiNEXT(tree_total__2018, type=1) + ggtitle("Total species (2018)")+ ggtitle("Total number species with rarefaction = 2 traps (2018)")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot8<-ggiNEXT(tree_total__2018, type=2) + ggtitle("Total species (2018)")+ ggtitle("Coverage estimatation with rarefaction = 2 traps") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10)) ##Coverage
plot14<-ggiNEXT(tree_total__2018, type=3) + geom_vline(xintercept = 0.65, linetype = "dashed") + ggtitle("Species estimation based on sample coverage") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

Tree_Plot_total_2018 <- ggarrange(plot7, plot8, plot14,  
                                  common.legend = TRUE, legend = "bottom",
                                  ncol = 3, nrow = 1)
ggsave("Tree_Plot_total_2018.png", width = 20, height = 10)

total2018_sc <- bind_rows(tree_total__2018$iNextEst) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD      = approx(SC, qD, xout = 0.60)$y,
    qD.LCL  = approx(SC, qD.LCL, xout = 0.60)$y,
    qD.UCL  = approx(SC, qD.UCL, xout = 0.60)$y,
    SC      = 0.60
  ) %>%
  ungroup()

total2018_trap <- bind_rows(tree_total__2018$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))

#tree species & xylobiont 2018
sp_per_trap_Rubra_xylobiont_2018<-dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Quercus rubra" & taxonomy_2018$xylobiont == "j" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Rubra_xylobiont_2018) <- sp_per_trap_Rubra_xylobiont_2018[,1]
sp_per_trap_Rubra_xylobiont_2018<-sp_per_trap_Rubra_xylobiont_2018[,-1]
sp_per_trap_Rubra_xylobiont_2018[sp_per_trap_Rubra_xylobiont_2018>=1]<-1

sp_per_trap_Fraxinus_xylobiont_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Fraxinus" & taxonomy_2018$xylobiont == "j" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Fraxinus_xylobiont_2018) <- sp_per_trap_Fraxinus_xylobiont_2018[,1]
sp_per_trap_Fraxinus_xylobiont_2018<-sp_per_trap_Fraxinus_xylobiont_2018[,-1]
sp_per_trap_Fraxinus_xylobiont_2018[sp_per_trap_Fraxinus_xylobiont_2018>=1]<-1

sp_per_trap_Quercus_xylobiont_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Quercus"& taxonomy_2018$xylobiont == "j" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Quercus_xylobiont_2018) <- sp_per_trap_Quercus_xylobiont_2018[,1]
sp_per_trap_Quercus_xylobiont_2018<-sp_per_trap_Quercus_xylobiont_2018[,-1]
sp_per_trap_Quercus_xylobiont_2018[sp_per_trap_Quercus_xylobiont_2018>=1]<-1

sp_per_trap_Acer_xylobiont_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Acer" & taxonomy_2018$xylobiont == "j" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Acer_xylobiont_2018) <- sp_per_trap_Acer_xylobiont_2018[,1]
sp_per_trap_Acer_xylobiont_2018<-sp_per_trap_Acer_xylobiont_2018[,-1]
sp_per_trap_Acer_xylobiont_2018[sp_per_trap_Acer_xylobiont_2018>=1]<-1

sp_per_trap_Tilia_xylobiont_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Tilia" & taxonomy_2018$xylobiont == "j" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Tilia_xylobiont_2018) <- sp_per_trap_Tilia_xylobiont_2018[,1]
sp_per_trap_Tilia_xylobiont_2018<-sp_per_trap_Tilia_xylobiont_2018[,-1]
sp_per_trap_Tilia_xylobiont_2018[sp_per_trap_Tilia_xylobiont_2018>=1]<-1

sp_per_trap_Ulmus_xylobiont_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Ulmus" & taxonomy_2018$xylobiont == "j" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Ulmus_xylobiont_2018) <- sp_per_trap_Ulmus_xylobiont_2018[,1]
sp_per_trap_Ulmus_xylobiont_2018<-sp_per_trap_Ulmus_xylobiont_2018[,-1]
sp_per_trap_Ulmus_xylobiont_2018[sp_per_trap_Ulmus_xylobiont_2018>=1]<-1 #turn into 0-1 values for absence/presence analysis

Forest_treesp_xylobiont_2018= list("Fraxinus" = sp_per_trap_Fraxinus_xylobiont_2018,"Quercus" = sp_per_trap_Quercus_xylobiont_2018, "Tilia" = sp_per_trap_Tilia_xylobiont_2018, "Quercus rubra" = sp_per_trap_Rubra_xylobiont_2018, "Acer" = sp_per_trap_Acer_xylobiont_2018, "Ulmus" = sp_per_trap_Ulmus_xylobiont_2018)

tree_xylobiont__2018 <- iNEXT(Forest_treesp_xylobiont_2018, q=1, datatype ="incidence_raw", endpoint=2)

plot9<-ggiNEXT(tree_xylobiont__2018, type=1) + ggtitle("Total number xylobiont species with rarefaction = 2 traps (2018)")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot10<-ggiNEXT(tree_xylobiont__2018, type=2) + ggtitle("Coverage estimatation with rarefaction = 2 traps)")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot22<-ggiNEXT(tree_xylobiont__2018, type=3) + ggtitle("Species estimation based on sample coverage")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

Tree_Plot_xylo_2018 <- ggarrange(plot9, plot10, plot22, 
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 3, nrow = 1)
ggsave("Tree_Plot_xylo_2018.png", width = 20, height = 10)

xylo2018_sc <- bind_rows(tree_xylobiont__2018$iNextEst) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD      = approx(SC, qD, xout = 0.60)$y,
    qD.LCL  = approx(SC, qD.LCL, xout = 0.60)$y,
    qD.UCL  = approx(SC, qD.UCL, xout = 0.60)$y,
    SC      = 0.60
  ) %>%
  ungroup()

xylo2018_trap <- bind_rows(tree_xylobiont__2018$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))

#red_list and tree species 2018
sp_per_trap_Rubra_red_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Quercus rubra"& taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Rubra_red_2018) <- sp_per_trap_Rubra_red_2018[,1]
sp_per_trap_Rubra_red_2018<-sp_per_trap_Rubra_red_2018[,-1]
sp_per_trap_Rubra_red_2018[sp_per_trap_Rubra_red_2018>=1]<-1

sp_per_trap_Fraxinus_red_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Fraxinus" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Fraxinus_red_2018) <- sp_per_trap_Fraxinus_red_2018[,1]
sp_per_trap_Fraxinus_red_2018<-sp_per_trap_Fraxinus_red_2018[,-1]
sp_per_trap_Fraxinus_red_2018[sp_per_trap_Fraxinus_red_2018>=1]<-1

sp_per_trap_Quercus_red_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Quercus" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Quercus_red_2018) <- sp_per_trap_Quercus_red_2018[,1]
sp_per_trap_Quercus_red_2018<-sp_per_trap_Quercus_red_2018[,-1]
sp_per_trap_Quercus_red_2018[sp_per_trap_Quercus_red_2018>=1]<-1

sp_per_trap_Acer_red_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Acer" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Acer_red_2018) <- sp_per_trap_Acer_red_2018[,1]
sp_per_trap_Acer_red_2018<-sp_per_trap_Acer_red_2018[,-1]
sp_per_trap_Acer_red_2018[sp_per_trap_Acer_red_2018>=1]<-1

sp_per_trap_Tilia_red_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Tilia" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Tilia_red_2018) <- sp_per_trap_Tilia_red_2018[,1]
sp_per_trap_Tilia_red_2018<-sp_per_trap_Tilia_red_2018[,-1]
sp_per_trap_Tilia_red_2018[sp_per_trap_Tilia_red_2018>=1]<-1

sp_per_trap_Ulmus_red_2018<-dcast( taxonomy_2018[taxonomy_2018$TreeSp=="Ulmus" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number")
rownames(sp_per_trap_Ulmus_red_2018) <- sp_per_trap_Ulmus_red_2018[,1]
sp_per_trap_Ulmus_red_2018<-sp_per_trap_Ulmus_red_2018[,-1]
sp_per_trap_Ulmus_red_2018[sp_per_trap_Ulmus_red_2018>=1]<-1

Forest_treesp_red_2018= list("Fraxinus" = sp_per_trap_Fraxinus_red_2018,"Quercus" = sp_per_trap_Quercus_red_2018, "Tilia" = sp_per_trap_Tilia_red_2018, "Quercus rubra" = sp_per_trap_Rubra_red_2018, "Acer" = sp_per_trap_Acer_red_2018, "Ulmus" = sp_per_trap_Ulmus_red_2018)

tree_red_2018 <- iNEXT(Forest_treesp_red_2018, q=1, datatype ="incidence_raw", endpoint=2)

plot11<-ggiNEXT(tree_red_2018, type=1) + ggtitle("Total number red-list species with rarefaction = 2 traps (2018)")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot12<-ggiNEXT(tree_red_2018, type=2) + ggtitle("Coverage estimatation with rarefaction = 2 traps")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot23<-ggiNEXT(tree_red_2018, type=3) + ggtitle("Species estimation based on sample coverage")+
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

Tree_Plot_red_2018 <- ggarrange(plot11, plot12,plot23,  
                                common.legend = TRUE, legend = "bottom",
                                ncol = 3, nrow = 1)
ggsave("Tree_Plot_red_2018.png", width = 20, height = 10)

redlist2018_sc <- bind_rows(tree_red_2018$iNextEst) %>%  ##2018 not enough data for all tree species to have numbers
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD      = approx(SC, qD, xout = 0.60)$y,
    qD.LCL  = approx(SC, qD.LCL, xout = 0.60)$y,
    qD.UCL  = approx(SC, qD.UCL, xout = 0.60)$y,
    SC      = 0.60
  ) %>%
  ungroup()

redlist2018_trap <- bind_rows(tree_red_2018$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))

### Beetles uniqueness total 2017

beetles_2017_treesp <- taxonomy_2017[taxonomy_2017$TreeSp != "Ground", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number))

distinct_2017 <- beetles_2017_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n()) 

beetles_2017_b <- beetles_2017_treesp %>%
  left_join(distinct_2017, by = c("TreeSp", "FHL")) %>%
  group_by(TreeSp, n_appearance) %>%
  summarise(n_beetlesp = n_distinct(FHL))

beetels_2017_b <- dcast( beetles_2017_b, n_appearance~TreeSp, value.var="n_beetlesp")

### Beetles uniqueness redlist 2017

red_2017_treesp <- taxonomy_2017[taxonomy_2017$TreeSp != "Ground" & taxonomy_2017$Rote.Liste != "n" & taxonomy_2017$Number != "0", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number))

red_distinct_2017 <- red_2017_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n()) 

red_2017_b <- red_2017_treesp %>%
  left_join(red_distinct_2017, by = c("TreeSp", "FHL")) %>%
  group_by(TreeSp, n_appearance) %>%
  summarise(n_beetlesp = n_distinct(FHL))%>%
  
  red_2017_dcast <- dcast( red_2017_b, n_appearance~TreeSp, value.var="n_beetlesp")
red_2017_dcast[is.na(red_2017_dcast)] <- 0

### Beetles uniqueness total 2018

beetles_2018_treesp <- taxonomy_2018[taxonomy_2018$TreeSp != "Ground", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number))

distinct_2018 <- beetles_2018_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n()) 

beetles_2018_b <- beetles_2018_treesp %>%
  left_join(distinct_2018, by = c("TreeSp", "FHL")) %>%
  group_by(TreeSp, n_appearance) %>%
  summarise(n_beetlesp = n_distinct(FHL))

beetels_2018_b <- dcast( beetles_2018_b, n_appearance~TreeSp, value.var="n_beetlesp")

### Beetles uniqueness redlist 2018

red_2018_treesp <- taxonomy_2018[taxonomy_2018$TreeSp != "Ground" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != "0", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number))

red_distinct_2018 <- red_2018_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n()) 

red_2018_b <- red_2018_treesp %>%
  left_join(red_distinct_2018, by = c("TreeSp", "FHL")) %>%
  group_by(TreeSp, n_appearance) %>%
  summarise(n_beetlesp = n_distinct(FHL))

red_2018_dcast <- dcast( red_2018_b, n_appearance~TreeSp, value.var="n_beetlesp")
red_2018_dcast[is.na(red_2018_dcast)] <- 0



