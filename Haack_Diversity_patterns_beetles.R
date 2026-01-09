###Response of common and rare beetle species to tree species and vertical stratification in a floodplain forest

#R script supplement to paper

# load packages and set working directory
setwd("D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2021\\Beta diversity paper\\Paper 2 Nov20")

#install.packages("knitr")
library(knitr)
library(ggplot2)
#install.packages("broom")
library(broom)
library(devtools)

#install.packages("xlsx")
#library(xlsx)

#install.packages("reshape2")
library("reshape2")

#install.packages("vegan")
library("vegan")

#install.packages("dendextend")
library("dendextend")

#install.packages("flextable")
library(flextable)

#install.packages("officer")
library(officer)

#install.packages("BiodiversityR")
#library(BiodiversityR)

library(dplyr)

library(tidyr)

#install.packages("BAT")
library(BAT)

library(dplyr)

library(gridExtra)

#install.packages("spaa")
library(spaa)

#install.packages("tinytex")
#tinytex::install_tinytex() # install TinyTeX

library(SpadeR)

#############################

# load data

coleo2017 <- read.csv2("LAK_2017_Coleoptera_21.04.2020.csv")

coleo2017$date <- as.Date(with(coleo2017, paste(coleo2017$End.Tag, coleo2017$End.Monat, coleo2017$Jahr, sep="-")), "%d-%m-%Y")

coleo2017$Anzahl<- as.numeric(coleo2017$Anzahl)

#make subset with only determined species
beetles_2017<-coleo2017[ which(coleo2017$Art!='NA'
                               & coleo2017$Art!='spec.'& coleo2017$Art!='???? femoralis'), ]
beetles_2017<-beetles_2017 %>% filter(!(grepl("spec.", Art)))
beetles_2017<-beetles_2017 %>% filter(!(grepl("cf.", Art)))

beetles_2017<-beetles_2017 %>% unite("Artname", Gattung:Art, sep=" ")

beetles_2017$Falle<-as.numeric(beetles_2017$Falle)


#make subset with only un-determined species
undet_taxa_2017<-coleo2017 %>% filter((grepl(c("spec.|cf.|NA"), Art)))

undet_taxa_2017$Anzahl<-as.numeric(undet_taxa_2017$Anzahl)

beetles_traits<-read.csv2("Species_traits_2016_2017_complete.csv")
beetles_2017_traits<-merge(x = beetles_2017, y = beetles_traits[ , c("Artname", "xylobiont", "general.rarity")], by="Artname")
beetles_2017_xylo<-beetles_2017_traits[which(beetles_2017_traits$xylobiont=="y"),]

non_xylo_sp<-unique(beetles_2017_traits$Artname[which(beetles_2017_traits$xylobiont=="n")])
xylo_sp_vector<-unique(beetles_2017_xylo$Artname[which(beetles_2017_xylo$xylobiont=="y")])

#write.csv(beetles_2017_xylo, "beetles_2017_xylo.csv")

########################

#read environmental dataset with only stratum and treesp as a beginning
crane.env.strat<- read.csv(file="crane.env.2017.csv", header=TRUE, sep=";")
crane.env.strat$Fallennummer<-as.numeric(crane.env.strat$Fallennummer)

crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)


beetles_strata<-merge(x = beetles_2017_xylo, y = crane.env.strat[ , c("Fallennummer", "Stratum")], by.x="Falle", by.y = "Fallennummer")

beetles_2017_xylo_treesp <- beetles_2017_xylo[ which(beetles_2017_xylo$Falle!=31&beetles_2017_xylo$Falle!=32&beetles_2017_xylo$Falle!=21&beetles_2017_xylo$Falle!=22&beetles_2017_xylo$Falle!=27&beetles_2017_xylo$Falle!=28&beetles_2017_xylo$Falle!=29&beetles_2017_xylo$Falle!=30&beetles_2017_xylo$Falle!=39&beetles_2017_xylo$Falle!=40&beetles_2017_xylo$Falle!=41&beetles_2017_xylo$Falle!=42&beetles_2017_xylo$Falle!=43&beetles_2017_xylo$Falle!=44),]

beetles_tree<-merge(x = beetles_2017_xylo_treesp, y = crane.env.tree[ , c("Fallennummer", "TreeSp")], by.x="Falle", by.y = "Fallennummer")

######################

#assign rarity

## Red list rarity

# rarity according to red list
beetles_strata$rare<-NA
beetles_tree$rare<-NA

beetles_strata$rare<-ifelse((beetles_strata$general.rarity=="1"|beetles_strata$general.rarity=="2"|beetles_strata$general.rarity=="3"), "y", "n")

beetles_tree$rare<-ifelse((beetles_tree$general.rarity=="1"|beetles_tree$general.rarity=="2"|beetles_tree$general.rarity=="3"), "y", "n")

## Octave rarity
#install.packages("gambin")
library('gambin')

beetles_octaves<- aggregate(Anzahl~Artname, beetles_2017_xylo, FUN=sum)
beetles_octaves_num<-beetles_octaves[,2]

oct<-create_octaves(beetles_octaves_num, subsample = 0)

#barplot(oct$species, names=oct$octave)

beetles_octaves$rare<-"n"
beetles_octaves$rare[beetles_octaves$Anzahl<=64]<-"y"


beetles_octaves_strata<-beetles_strata[,c(1:17)]

beetles_octaves_tree<-beetles_tree[,c(1:17)]

beetles_octaves_strata<-merge(x = beetles_octaves_strata, y = beetles_octaves[ , c("Artname", "rare")], by.x="Artname", by.y = "Artname")

beetles_octaves_tree<-merge(x = beetles_octaves_tree, y = beetles_octaves[ , c("Artname", "rare")], by.x="Artname", by.y = "Artname")


################################

### Alpha diversity of xylobiont beetle communities


#red list

library(iNEXT)

#iNEXT estimation

#long to wide format of data

#strata & common sp
sp_per_trap_ng_common<-dcast( beetles_strata[beetles_strata$rare=="n"&beetles_strata$Stratum=="Understorey",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_ng_common) <- sp_per_trap_ng_common[,1]
sp_per_trap_ng_common<-sp_per_trap_ng_common[,-1]
sp_per_trap_ng_common[sp_per_trap_ng_common>=1]<-1

sp_per_trap_lc_common<-dcast( beetles_strata[beetles_strata$rare=="n"&beetles_strata$Stratum=="Lower",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_lc_common) <- sp_per_trap_lc_common[,1]
sp_per_trap_lc_common<-sp_per_trap_lc_common[,-1]
sp_per_trap_lc_common[sp_per_trap_lc_common>=1]<-1


sp_per_trap_uc_common<-dcast( beetles_strata[beetles_strata$rare=="n"&beetles_strata$Stratum=="Upper",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_uc_common) <- sp_per_trap_uc_common[,1]
sp_per_trap_uc_common<-sp_per_trap_uc_common[,-1]
sp_per_trap_uc_common[sp_per_trap_uc_common>=1]<-1


#strata &rare sp
sp_per_trap_ng_rare<-dcast( beetles_strata[beetles_strata$rare=="y"&beetles_strata$Stratum=="Understorey",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_ng_rare) <- sp_per_trap_ng_rare[,1]
sp_per_trap_ng_rare<-sp_per_trap_ng_rare[,-1]
sp_per_trap_ng_rare[sp_per_trap_ng_rare>=1]<-1


sp_per_trap_lc_rare<-dcast( beetles_strata[beetles_strata$rare=="y"&beetles_strata$Stratum=="Lower",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_lc_rare) <- sp_per_trap_lc_rare[,1]
sp_per_trap_lc_rare<-sp_per_trap_lc_rare[,-1]
sp_per_trap_lc_rare[sp_per_trap_lc_rare>=1]<-1

sp_per_trap_uc_rare<-dcast( beetles_strata[beetles_strata$rare=="y"&beetles_strata$Stratum=="Upper",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_uc_rare) <- sp_per_trap_uc_rare[,1]
sp_per_trap_uc_rare<-sp_per_trap_uc_rare[,-1]
sp_per_trap_uc_rare[sp_per_trap_uc_rare>=1]<-1

#tree & common sp
sp_per_trap_Tilia_common<-dcast( beetles_tree[beetles_tree$rare=="n"&beetles_tree$TreeSp=="Tilia",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Tilia_common) <- sp_per_trap_Tilia_common[,1]
sp_per_trap_Tilia_common<-sp_per_trap_Tilia_common[,-1]
sp_per_trap_Tilia_common[sp_per_trap_Tilia_common>=1]<-1

sp_per_trap_Quercus_common<-dcast( beetles_tree[beetles_tree$rare=="n"&beetles_tree$TreeSp=="Quercus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Quercus_common) <- sp_per_trap_Quercus_common[,1]
sp_per_trap_Quercus_common<-sp_per_trap_Quercus_common[,-1]
sp_per_trap_Quercus_common[sp_per_trap_Quercus_common>=1]<-1

sp_per_trap_Fraxinus_common<-dcast( beetles_tree[beetles_tree$rare=="n"&beetles_tree$TreeSp=="Fraxinus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Fraxinus_common) <- sp_per_trap_Fraxinus_common[,1]
sp_per_trap_Fraxinus_common<-sp_per_trap_Fraxinus_common[,-1]
sp_per_trap_Fraxinus_common[sp_per_trap_Fraxinus_common>=1]<-1

#tree &rare sp
sp_per_trap_Tilia_rare<-dcast( beetles_tree[beetles_tree$rare=="y"&beetles_tree$TreeSp=="Tilia",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Tilia_rare) <- sp_per_trap_Tilia_rare[,1]
sp_per_trap_Tilia_rare<-sp_per_trap_Tilia_rare[,-1]
sp_per_trap_Tilia_rare[sp_per_trap_Tilia_rare>=1]<-1

sp_per_trap_Quercus_rare<-dcast( beetles_tree[beetles_tree$rare=="y"&beetles_tree$TreeSp=="Quercus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Quercus_rare) <- sp_per_trap_Quercus_rare[,1]
sp_per_trap_Quercus_rare<-sp_per_trap_Quercus_rare[,-1]
sp_per_trap_Quercus_rare[sp_per_trap_Quercus_rare>=1]<-1

sp_per_trap_Fraxinus_rare<-dcast( beetles_tree[beetles_tree$rare=="y"&beetles_tree$TreeSp=="Fraxinus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Fraxinus_rare) <- sp_per_trap_Fraxinus_rare[,1]
sp_per_trap_Fraxinus_rare<-sp_per_trap_Fraxinus_rare[,-1]
sp_per_trap_Fraxinus_rare[sp_per_trap_Fraxinus_rare>=1]<-1

#### now the real species richness estimation takes place!
##########################

Forest_strata_common= list("near ground" = sp_per_trap_ng_common, "lower canopy" = sp_per_trap_lc_common, "upper canopy" = sp_per_trap_uc_common)

Forest_strata_rare= list("near ground" = sp_per_trap_ng_rare, "lower canopy" = sp_per_trap_lc_rare, "upper canopy" = sp_per_trap_uc_rare)

Forest_treesp_common= list("Tilia" = sp_per_trap_Tilia_common, "Quercus" = sp_per_trap_Quercus_common, "Fraxinus" = sp_per_trap_Fraxinus_common)

Forest_treesp_rare= list("Tilia" = sp_per_trap_Tilia_rare, "Quercus" = sp_per_trap_Quercus_rare, "Fraxinus" = sp_per_trap_Fraxinus_rare)

###### non-asymptotic analyses ######
#install.packages ("iNEXT")
library(iNEXT)
library (ggplot2)
library(ggpubr)

#for strata and common sp
out_strata_common <- iNEXT(Forest_strata_common, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_strata_common
# gives species richness AND diversity estimates!

#for strata and rare sp
out_strata_rare <- iNEXT(Forest_strata_rare, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_strata_rare
# gives species richness AND diversity estimates!

#for tree sp and common sp
out_tree_common <- iNEXT(Forest_treesp_common, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_tree_common
# gives species richness AND diversity estimates!

#for tree sp and rare sp
out_tree_rare <- iNEXT(Forest_treesp_rare, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_tree_rare
# gives species richness AND diversity estimates!

plot1<-ggiNEXT(out_strata_common, type=1) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot2<-ggiNEXT(out_strata_common, type=2) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot3<-ggiNEXT(out_strata_common, type=3) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

plot4<-ggiNEXT(out_strata_rare, type=1) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot5<-ggiNEXT(out_strata_rare, type=2) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot6<-ggiNEXT(out_strata_rare, type=3) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

plot7<-ggiNEXT(out_tree_common, type=1) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot8<-ggiNEXT(out_tree_common, type=2) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot9<-ggiNEXT(out_tree_common, type=3) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

plot10<-ggiNEXT(out_tree_rare, type=1) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot11<-ggiNEXT(out_tree_rare, type=2) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot12<-ggiNEXT(out_tree_rare, type=3) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

Strata_Plot <- ggarrange(plot1, plot2, plot4, plot5,
                         common.legend = TRUE, legend = "bottom",
                         ncol = 2, nrow = 2)

ggsave("Strata_Plot.png", width = 10, height = 10)

Tree_Plot <- ggarrange(plot7, plot8, plot10, plot11,
                       common.legend = TRUE, legend = "bottom",
                       ncol = 2, nrow = 2)
ggsave("Tree_Plot.png", width = 10, height = 10)


#test significance
Rarity<-c("common", "common", "common", "rare", "rare", "rare","common", "common", "common", "rare", "rare", "rare")
Factor<-c("Strata", "Strata", "Strata", "Strata", "Strata", "Strata", "Tree species", "Tree species", "Tree species", "Tree species", "Tree species", "Tree species")
Unit<-c("Near ground", "Lower canopy", "Upper canopy","Near ground", "Lower canopy", "Upper canopy", "Quercus", "Tilia", "Fraxinus" , "Quercus", "Tilia", "Fraxinus")
Estimate<-c(
  ((out_strata_common[["iNextEst"]][["near ground"]]$qD)[out_strata_common[["iNextEst"]][["near ground"]]$t==30]),
  ((out_strata_common[["iNextEst"]][["lower canopy"]]$qD)[out_strata_common[["iNextEst"]][["lower canopy"]]$t==30]),
  ((out_strata_common[["iNextEst"]][["upper canopy"]]$qD)[out_strata_common[["iNextEst"]][["upper canopy"]]$t==30]), 
  ((out_strata_rare[["iNextEst"]][["near ground"]]$qD)[out_strata_rare[["iNextEst"]][["near ground"]]$t==30]),
  ((out_strata_rare[["iNextEst"]][["lower canopy"]]$qD)[out_strata_rare[["iNextEst"]][["lower canopy"]]$t==30]),
  ((out_strata_rare[["iNextEst"]][["upper canopy"]]$qD)[out_strata_rare[["iNextEst"]][["upper canopy"]]$t==30]),
  ((out_tree_common[["iNextEst"]][["Quercus"]]$qD)[out_tree_common[["iNextEst"]][["Quercus"]]$t==30]),
  ((out_tree_common[["iNextEst"]][["Tilia"]]$qD)[out_tree_common[["iNextEst"]][["Tilia"]]$t==30]),
  ((out_tree_common[["iNextEst"]][["Fraxinus"]]$qD)[out_tree_common[["iNextEst"]][["Fraxinus"]]$t==30]),
  ((out_tree_rare[["iNextEst"]][["Quercus"]]$qD)[out_tree_rare[["iNextEst"]][["Quercus"]]$t==30]),
  ((out_tree_rare[["iNextEst"]][["Tilia"]]$qD)[out_tree_rare[["iNextEst"]][["Tilia"]]$t==30]),
  ((out_tree_rare[["iNextEst"]][["Fraxinus"]]$qD)[out_tree_rare[["iNextEst"]][["Fraxinus"]]$t==30])
)

Alpha_frame<-data.frame(Factor, Unit, Rarity, Estimate)

aov.res <- aov(Estimate~ Rarity + Unit, data=Alpha_frame)
summary(aov.res)

aov.res<-as.data.frame(summary(aov.res)[[1]])
names(aov.res)[2] <- "Sum.Sq"
names(aov.res)[3] <- "Mean.Sq"
names(aov.res)[4] <- "F.value"
names(aov.res)[5] <- "P"


#octaves

#iNEXT estimation

#long to wide format of data

#strata & common sp
sp_per_trap_ng_common_oct<-dcast( beetles_octaves_strata[beetles_octaves_strata$rare=="n"&beetles_octaves_strata$Stratum=="Understorey",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_ng_common_oct) <- sp_per_trap_ng_common_oct[,1]
sp_per_trap_ng_common_oct<-sp_per_trap_ng_common_oct[,-1]
sp_per_trap_ng_common_oct[sp_per_trap_ng_common_oct>=1]<-1

sp_per_trap_lc_common_oct<-dcast( beetles_octaves_strata[beetles_octaves_strata$rare=="n"&beetles_octaves_strata$Stratum=="Lower",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_lc_common_oct) <- sp_per_trap_lc_common_oct[,1]
sp_per_trap_lc_common_oct<-sp_per_trap_lc_common_oct[,-1]
sp_per_trap_lc_common_oct[sp_per_trap_lc_common_oct>=1]<-1


sp_per_trap_uc_common_oct<-dcast( beetles_octaves_strata[beetles_octaves_strata$rare=="n"&beetles_octaves_strata$Stratum=="Upper",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_uc_common_oct) <- sp_per_trap_uc_common_oct[,1]
sp_per_trap_uc_common_oct<-sp_per_trap_uc_common_oct[,-1]
sp_per_trap_uc_common_oct[sp_per_trap_uc_common_oct>=1]<-1


#strata &rare sp
sp_per_trap_ng_rare<-dcast( beetles_octaves_strata[beetles_octaves_strata$rare=="y"&beetles_octaves_strata$Stratum=="Understorey",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_ng_rare) <- sp_per_trap_ng_rare[,1]
sp_per_trap_ng_rare<-sp_per_trap_ng_rare[,-1]
sp_per_trap_ng_rare[sp_per_trap_ng_rare>=1]<-1


sp_per_trap_lc_rare<-dcast( beetles_octaves_strata[beetles_octaves_strata$rare=="y"&beetles_octaves_strata$Stratum=="Lower",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_lc_rare) <- sp_per_trap_lc_rare[,1]
sp_per_trap_lc_rare<-sp_per_trap_lc_rare[,-1]
sp_per_trap_lc_rare[sp_per_trap_lc_rare>=1]<-1

sp_per_trap_uc_rare_oct<-dcast( beetles_octaves_strata[beetles_octaves_strata$rare=="y"&beetles_octaves_strata$Stratum=="Upper",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_uc_rare_oct) <- sp_per_trap_uc_rare_oct[,1]
sp_per_trap_uc_rare_oct<-sp_per_trap_uc_rare_oct[,-1]
sp_per_trap_uc_rare_oct[sp_per_trap_uc_rare_oct>=1]<-1

#tree & common sp
sp_per_trap_Tilia_common_oct<-dcast( beetles_octaves_tree[beetles_octaves_tree$rare=="n"&beetles_octaves_tree$TreeSp=="Tilia",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Tilia_common_oct) <- sp_per_trap_Tilia_common_oct[,1]
sp_per_trap_Tilia_common_oct<-sp_per_trap_Tilia_common_oct[,-1]
sp_per_trap_Tilia_common_oct[sp_per_trap_Tilia_common_oct>=1]<-1

sp_per_trap_Quercus_common_oct<-dcast( beetles_octaves_tree[beetles_octaves_tree$rare=="n"&beetles_octaves_tree$TreeSp=="Quercus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Quercus_common_oct) <- sp_per_trap_Quercus_common_oct[,1]
sp_per_trap_Quercus_common_oct<-sp_per_trap_Quercus_common_oct[,-1]
sp_per_trap_Quercus_common_oct[sp_per_trap_Quercus_common_oct>=1]<-1

sp_per_trap_Fraxinus_common_oct<-dcast( beetles_octaves_tree[beetles_octaves_tree$rare=="n"&beetles_octaves_tree$TreeSp=="Fraxinus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Fraxinus_common_oct) <- sp_per_trap_Fraxinus_common_oct[,1]
sp_per_trap_Fraxinus_common_oct<-sp_per_trap_Fraxinus_common_oct[,-1]
sp_per_trap_Fraxinus_common_oct[sp_per_trap_Fraxinus_common_oct>=1]<-1

#tree &rare sp
sp_per_trap_Tilia_rare_oct<-dcast( beetles_octaves_tree[beetles_octaves_tree$rare=="y"&beetles_octaves_tree$TreeSp=="Tilia",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Tilia_rare_oct) <- sp_per_trap_Tilia_rare_oct[,1]
sp_per_trap_Tilia_rare_oct<-sp_per_trap_Tilia_rare_oct[,-1]
sp_per_trap_Tilia_rare_oct[sp_per_trap_Tilia_rare_oct>=1]<-1

sp_per_trap_Quercus_rare_oct<-dcast( beetles_octaves_tree[beetles_octaves_tree$rare=="y"&beetles_octaves_tree$TreeSp=="Quercus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Quercus_rare_oct) <- sp_per_trap_Quercus_rare_oct[,1]
sp_per_trap_Quercus_rare_oct<-sp_per_trap_Quercus_rare_oct[,-1]
sp_per_trap_Quercus_rare_oct[sp_per_trap_Quercus_rare_oct>=1]<-1

sp_per_trap_Fraxinus_rare_oct<-dcast( beetles_octaves_tree[beetles_octaves_tree$rare=="y"&beetles_octaves_tree$TreeSp=="Fraxinus",], Artname~Falle, value.var="Anzahl")
rownames(sp_per_trap_Fraxinus_rare_oct) <- sp_per_trap_Fraxinus_rare_oct[,1]
sp_per_trap_Fraxinus_rare_oct<-sp_per_trap_Fraxinus_rare_oct[,-1]
sp_per_trap_Fraxinus_rare_oct[sp_per_trap_Fraxinus_rare_oct>=1]<-1

#### now the real species richness estimation takes place!
##########################

Forest_strata_common_oct= list("near ground" = sp_per_trap_ng_common_oct, "lower canopy" = sp_per_trap_lc_common_oct, "upper canopy" = sp_per_trap_uc_common_oct)

Forest_strata_rare_oct= list("near ground" = sp_per_trap_ng_rare, "lower canopy" = sp_per_trap_lc_rare, "upper canopy" = sp_per_trap_uc_rare_oct)

Forest_treesp_common_oct= list("Tilia" = sp_per_trap_Tilia_common_oct, "Quercus" = sp_per_trap_Quercus_common_oct, "Fraxinus" = sp_per_trap_Fraxinus_common_oct)

Forest_treesp_rare_oct= list("Tilia" = sp_per_trap_Tilia_rare_oct, "Quercus" = sp_per_trap_Quercus_rare_oct, "Fraxinus" = sp_per_trap_Fraxinus_rare_oct)

###### non-asymptotic analyses ######
#install.packages ("iNEXT")
library(iNEXT)
library (ggplot2)
library(ggpubr)

#for strata and common sp
out_strata_common_oct <- iNEXT(Forest_strata_common_oct, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_strata_common_oct
# gives species richness AND diversity estimates!

#for strata and rare sp
out_strata_rare_oct <- iNEXT(Forest_strata_rare_oct, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_strata_rare_oct
# gives species richness AND diversity estimates!

#for tree sp and common sp
out_tree_common_oct <- iNEXT(Forest_treesp_common_oct, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_tree_common_oct
# gives species richness AND diversity estimates!

#for tree sp and rare sp
out_tree_rare_oct <- iNEXT(Forest_treesp_rare_oct, q=0, datatype ="incidence_raw", endpoint=30)
# show detailed output
out_tree_rare_oct
# gives species richness AND diversity estimates!

plot13<-ggiNEXT(out_strata_common_oct, type=1) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot14<-ggiNEXT(out_strata_common_oct, type=2) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot15<-ggiNEXT(out_strata_common_oct, type=3) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

plot16<-ggiNEXT(out_strata_rare_oct, type=1) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot17<-ggiNEXT(out_strata_rare_oct, type=2) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot18<-ggiNEXT(out_strata_rare_oct, type=3) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

plot19<-ggiNEXT(out_tree_common_oct, type=1) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot20<-ggiNEXT(out_tree_common_oct, type=2) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot21<-ggiNEXT(out_tree_common_oct, type=3) + ggtitle("Common species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

plot22<-ggiNEXT(out_tree_rare_oct, type=1) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot23<-ggiNEXT(out_tree_rare_oct, type=2) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))
plot24<-ggiNEXT(out_tree_rare_oct, type=3) + ggtitle("Rare species")+
  theme(plot.title = element_text(hjust = 0.5, size =12), text = element_text(size=10))

Strata_Plot_oct <- ggarrange(plot13, plot14, plot16, plot17, 
                             common.legend = TRUE, legend = "bottom",
                             ncol = 2, nrow = 2)

ggsave("Strata_Plot_oct.png", width = 10, height = 10)

Tree_Plot_oct <- ggarrange(plot19, plot20, plot22, plot23, 
                           common.legend = TRUE, legend = "bottom",
                           ncol = 2, nrow = 2)
ggsave("Tree_Plot_oct.png", width = 10, height = 10)


#octaves

#test significance
Estimate_oct<-c(
  ((out_strata_common_oct[["iNextEst"]][["near ground"]]$qD)[out_tree_common_oct[["iNextEst"]][["near ground"]]$t==30]),
  ((out_strata_common_oct[["iNextEst"]][["lower canopy"]]$qD)[out_tree_common_oct[["iNextEst"]][["lower canopy"]]$t==30]),
  ((out_strata_common_oct[["iNextEst"]][["upper canopy"]]$qD)[out_tree_common_oct[["iNextEst"]][["upper canopy"]]$t==30]), 
  ((out_strata_rare_oct[["iNextEst"]][["near ground"]]$qD)[out_tree_rare_oct[["iNextEst"]][["near ground"]]$t==30]),
  ((out_strata_rare_oct[["iNextEst"]][["lower canopy"]]$qD)[out_tree_rare_oct[["iNextEst"]][["lower canopy"]]$t==30]),
  ((out_strata_rare_oct[["iNextEst"]][["upper canopy"]]$qD)[out_tree_rare_oct[["iNextEst"]][["upper canopy"]]$t==30]),
  ((out_tree_common_oct[["iNextEst"]][["Quercus"]]$qD)[out_tree_common[["iNextEst"]][["Quercus"]]$t==30]),
  ((out_tree_common_oct[["iNextEst"]][["Tilia"]]$qD)[out_tree_common[["iNextEst"]][["Tilia"]]$t==30]),
  ((out_tree_common_oct[["iNextEst"]][["Fraxinus"]]$qD)[out_tree_common[["iNextEst"]][["Fraxinus"]]$t==30]),
  ((out_tree_rare_oct[["iNextEst"]][["Quercus"]]$qD)[out_tree_rare_oct[["iNextEst"]][["Quercus"]]$t==30]),
  ((out_tree_rare_oct[["iNextEst"]][["Tilia"]]$qD)[out_tree_rare_oct[["iNextEst"]][["Tilia"]]$t==30]),
  ((out_tree_rare_oct[["iNextEst"]][["Fraxinus"]]$qD)[out_tree_rare_oct[["iNextEst"]][["Fraxinus"]]$t==30])
)

Alpha_frame_oct<-data.frame(Factor, Unit, Rarity, Estimate_oct)

aov.res.oct <- aov(Estimate_oct~ Rarity + Unit, data=Alpha_frame_oct)
summary(aov.res.oct)

aov.res.oct<-as.data.frame(summary(aov.res.oct)[[1]])
names(aov.res.oct)[2] <- "Sum.Sq"
names(aov.res.oct)[3] <- "Mean.Sq"
names(aov.res.oct)[4] <- "F.value"
names(aov.res.oct)[5] <- "P"

##############################

### Beta diversity of xylobiont beetle communities

#red list
#calculate for rare beetle spp
beetles_rare_strat <- aggregate(beetles_strata$Anzahl[beetles_strata$rare=="y"],
                                by = list(beetles_strata$Falle[beetles_strata$rare=="y"], beetles_strata$Species[beetles_strata$rare=="y"]),
                                FUN = sum) #aggregate per trap and species
beetles_rare_strat$Group.1<- as.numeric(beetles_rare_strat$Group.1)

beetles_rare_treesp <- aggregate(beetles_tree$Anzahl[beetles_tree$rare=="y"],
                                 by = list(beetles_tree$Falle[beetles_tree$rare=="y"], beetles_tree$Species[beetles_tree$rare=="y"]),
                                 FUN = sum) #aggregate per trap and species
beetles_rare_treesp$Group.1<- as.numeric(beetles_rare_treesp$Group.1)

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_rare_strat_1<-dcast(beetles_rare_strat, Group.2~Group.1)
beetles_rare_strat_1[is.na(beetles_rare_strat_1)] <- 0
rownames(beetles_rare_strat_1) <- beetles_rare_strat_1[,1]
beetles_rare_strat_1 <- beetles_rare_strat_1[,-1]

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_rare_treesp_1<-dcast(beetles_rare_treesp, Group.2~Group.1)
beetles_rare_treesp_1[is.na(beetles_rare_treesp_1)] <- 0
rownames(beetles_rare_treesp_1) <- beetles_rare_treesp_1[,1]
beetles_rare_treesp_1 <- beetles_rare_treesp_1[,-1]

###calculate similarity matrix (regional species-overlap indices) using anne chaos spader

#calculate for common beetle spp
beetles_common_strat <- aggregate(beetles_strata$Anzahl[beetles_strata$rare=="n"],
                                  by = list(beetles_strata$Falle[beetles_strata$rare=="n"], beetles_strata$Species[beetles_strata$rare=="n"]),
                                  FUN = sum) #aggregate per trap and species
beetles_common_strat$Group.1<- as.numeric(beetles_common_strat$Group.1)

beetles_common_treesp <- aggregate(beetles_tree$Anzahl[beetles_tree$rare=="n"],
                                   by = list(beetles_tree$Falle[beetles_tree$rare=="n"], beetles_tree$Species[beetles_tree$rare=="n"]),
                                   FUN = sum) #aggregate per trap and species
beetles_common_treesp$Group.1<- as.numeric(beetles_common_treesp$Group.1)

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_common_strat_1<-dcast(beetles_common_strat, Group.2~Group.1)
beetles_common_strat_1[is.na(beetles_common_strat_1)] <- 0
rownames(beetles_common_strat_1) <- beetles_common_strat_1[,1]
beetles_common_strat_1 <- beetles_common_strat_1[,-1]

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_common_treesp_1<-dcast(beetles_common_treesp, Group.2~Group.1)
beetles_common_treesp_1[is.na(beetles_common_treesp_1)] <- 0
rownames(beetles_common_treesp_1) <- beetles_common_treesp_1[,1]
beetles_common_treesp_1 <- beetles_common_treesp_1[,-1]

###calculate similarity matrix (regional species-overlap indices) using anne chaos spader

common_Chao_Strata<-SimilarityMult(beetles_common_strat_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
common_Beta_Strata_Chao<- as.data.frame(common_Chao_Strata[["similarity.matrix"]])
names(common_Beta_Strata_Chao)<-c(1:44)
common_beta_strata_matrix<-as.matrix((common_Beta_Strata_Chao))

#create dissimilarity from similarity 
einser_matrix<-matrix(1, nrow = 44, ncol = 44)

common_Beta_Strata_Chao<-einser_matrix-common_beta_strata_matrix

#save as csv
write.csv(common_Beta_Strata_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2021\\Beta diversity paper\\Paper 2 Nov20\\common_Beta_Strata_Chao.csv", row.names = FALSE)



common_Chao_Tree<-SimilarityMult(beetles_common_treesp_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
common_Beta_Tree_Chao<- as.data.frame(common_Chao_Tree[["similarity.matrix"]])
names(common_Beta_Tree_Chao)<-c(1:30)
common_beta_tree_matrix<-as.matrix((common_Beta_Tree_Chao))

#create dissimilarity from similarity 
einser_matrix<-matrix(1, nrow = 30, ncol = 30)

common_Beta_Tree_Chao<-einser_matrix-common_beta_tree_matrix

#save as csv
write.csv(common_Beta_Tree_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2021\\Beta diversity paper\\Paper 2 Nov20\\common_Beta_Tree_Chao.csv", row.names = FALSE)



# rare and strata
rare_Chao_Strata<-SimilarityMult(beetles_rare_strat_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
rare_Beta_Strata_Chao<- as.data.frame(rare_Chao_Strata[["similarity.matrix"]])
names(rare_Beta_Strata_Chao)<-c(1:44)
rare_beta_strata_matrix<-as.matrix((rare_Beta_Strata_Chao))

#create dissimilarity from similarity 
einser_matrix<-matrix(1, nrow = 44, ncol = 44)

rare_Beta_Strata_Chao<-einser_matrix-rare_beta_strata_matrix

#save as csv
write.csv(rare_Beta_Strata_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2020\\Beta diversity paper\\Paper 2 Nov20\\rare_Beta_Strata_Chao.csv", row.names = FALSE)



#rare and treesp
rare_Chao_Tree<-SimilarityMult(beetles_rare_treesp_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
rare_Beta_Tree_Chao<- as.data.frame(rare_Chao_Tree[["similarity.matrix"]])
names(rare_Beta_Tree_Chao)<-c(1:30)
rare_beta_tree_matrix<-as.matrix((rare_Beta_Tree_Chao))


#create dissimilarity from similarity 
einser_matrix<-matrix(1, nrow = 30, ncol = 30)

rare_Beta_Tree_Chao<-einser_matrix-rare_beta_tree_matrix

#save as csv
write.csv(rare_Beta_Tree_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2020\\Beta diversity paper\\Paper 2 Nov20\\rare_Beta_Tree_Chao.csv", row.names = FALSE)


#octaves

#calculate for rare beetle spp
beetles_octaves_rare_strat <- aggregate(beetles_octaves_strata$Anzahl[beetles_octaves_strata$rare=="y"],
                                        by = list(beetles_octaves_strata$Falle[beetles_octaves_strata$rare=="y"], beetles_octaves_strata$Species[beetles_octaves_strata$rare=="y"]),
                                        FUN = sum) #aggregate per trap and species
beetles_octaves_rare_strat$Group.1<- as.numeric(beetles_octaves_rare_strat$Group.1)

beetles_octaves_rare_treesp <- aggregate(beetles_octaves_tree$Anzahl[beetles_octaves_tree$rare=="y"],
                                         by = list(beetles_octaves_tree$Falle[beetles_octaves_tree$rare=="y"], beetles_octaves_tree$Species[beetles_octaves_tree$rare=="y"]),
                                         FUN = sum) #aggregate per trap and species
beetles_octaves_rare_treesp$Group.1<- as.numeric(beetles_octaves_rare_treesp$Group.1)

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_octaves_rare_strat_1<-dcast(beetles_octaves_rare_strat, Group.2~Group.1)
beetles_octaves_rare_strat_1[is.na(beetles_octaves_rare_strat_1)] <- 0
rownames(beetles_octaves_rare_strat_1) <- beetles_octaves_rare_strat_1[,1]
beetles_octaves_rare_strat_1 <- beetles_octaves_rare_strat_1[,-1]

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_octaves_rare_treesp_1<-dcast(beetles_octaves_rare_treesp, Group.2~Group.1)
beetles_octaves_rare_treesp_1[is.na(beetles_octaves_rare_treesp_1)] <- 0
rownames(beetles_octaves_rare_treesp_1) <- beetles_octaves_rare_treesp_1[,1]
beetles_octaves_rare_treesp_1 <- beetles_octaves_rare_treesp_1[,-1]



#calculate for common beetle spp
beetles_octaves_common_strat <- aggregate(beetles_octaves_strata$Anzahl[beetles_octaves_strata$rare=="n"],
                                          by = list(beetles_octaves_strata$Falle[beetles_octaves_strata$rare=="n"], beetles_octaves_strata$Species[beetles_octaves_strata$rare=="n"]),
                                          FUN = sum) #aggregate per trap and species
beetles_octaves_common_strat$Group.1<- as.numeric(beetles_octaves_common_strat$Group.1)

beetles_octaves_common_treesp <- aggregate(beetles_octaves_tree$Anzahl[beetles_octaves_tree$rare=="n"],
                                           by = list(beetles_octaves_tree$Falle[beetles_octaves_tree$rare=="n"], beetles_octaves_tree$Species[beetles_octaves_tree$rare=="n"]),
                                           FUN = sum) #aggregate per trap and species
beetles_octaves_common_treesp$Group.1<- as.numeric(beetles_octaves_common_treesp$Group.1)

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_octaves_common_strat_1<-dcast(beetles_octaves_common_strat, Group.2~Group.1)
beetles_octaves_common_strat_1[is.na(beetles_octaves_common_strat_1)] <- 0
rownames(beetles_octaves_common_strat_1) <- beetles_octaves_common_strat_1[,1]
beetles_octaves_common_strat_1 <- beetles_octaves_common_strat_1[,-1]

#create abundance matrix equivalent to SpadeRs SimilarityMultData$Abu
beetles_octaves_common_treesp_1<-dcast(beetles_octaves_common_treesp, Group.2~Group.1)
beetles_octaves_common_treesp_1[is.na(beetles_octaves_common_treesp_1)] <- 0
rownames(beetles_octaves_common_treesp_1) <- beetles_octaves_common_treesp_1[,1]
beetles_octaves_common_treesp_1 <- beetles_octaves_common_treesp_1[,-1]

###calculate similarity matrix (regional species-overlap indices) using anne chaos spader

common_octaves_Chao_Strata<-SimilarityMult(beetles_octaves_common_strat_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
common_octaves_Beta_Strata_Chao<- as.data.frame(common_octaves_Chao_Strata[["similarity.matrix"]])
names(common_octaves_Beta_Strata_Chao)<-c(1:44)
common_octaves_beta_strata_matrix<-as.matrix((common_octaves_Beta_Strata_Chao))

#create dissimilarity from similarity ---??? geht das?
einser_matrix<-matrix(1, nrow = 44, ncol = 44)

common_octaves_Beta_Strata_Chao<-einser_matrix-common_octaves_beta_strata_matrix

#save as csv
write.csv(common_octaves_Beta_Strata_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2020\\Beta diversity paper\\Paper 2 Nov20\\common_octaves_Beta_Strata_Chao.csv", row.names = FALSE)


common_octaves_Chao_Tree<-SimilarityMult(beetles_octaves_common_treesp_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
common_octaves_Beta_Tree_Chao<- as.data.frame(common_octaves_Chao_Tree[["similarity.matrix"]])
names(common_octaves_Beta_Tree_Chao)<-c(1:30)
common_octaves_beta_tree_matrix<-as.matrix((common_octaves_Beta_Tree_Chao))

#create dissimilarity from similarity ---??? geht das?
einser_matrix<-matrix(1, nrow = 30, ncol = 30)

common_octaves_Beta_Tree_Chao<-einser_matrix-common_octaves_beta_tree_matrix

#save as csv
write.csv(common_octaves_Beta_Tree_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2020\\Beta diversity paper\\Paper 2 Nov20\\common_octaves_Beta_Tree_Chao.csv", row.names = FALSE)



# rare and strata
rare_octaves_Chao_Strata<-SimilarityMult(beetles_octaves_rare_strat_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
rare_octaves_Beta_Strata_Chao<- as.data.frame(rare_octaves_Chao_Strata[["similarity.matrix"]])
names(rare_octaves_Beta_Strata_Chao)<-c(1:44)
rare_octaves_beta_strata_matrix<-as.matrix((rare_octaves_Beta_Strata_Chao))


#create dissimilarity from similarity 
einser_matrix<-matrix(1, nrow = 44, ncol = 44)
rare_octaves_Beta_Strata_Chao<-einser_matrix-rare_octaves_beta_strata_matrix

#save as csv
write.csv(rare_octaves_Beta_Strata_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2020\\Beta diversity paper\\Paper 2 Nov20\\rare_octaves_Beta_Strata_Chao.csv", row.names = FALSE)



#rare and treesp
rare_octaves_Chao_Tree<-SimilarityMult(beetles_octaves_rare_treesp_1, datatype = c("abundance"), q = 1, nboot = 200, goal = "absolute")
rare_octaves_Beta_Tree_Chao<- as.data.frame(rare_octaves_Chao_Tree[["similarity.matrix"]])
names(rare_octaves_Beta_Tree_Chao)<-c(1:30)
rare_octaves_beta_tree_matrix<-as.matrix((rare_octaves_Beta_Tree_Chao))


#create dissimilarity from similarity 
einser_matrix<-matrix(1, nrow = 30, ncol = 30)
rare_octaves_Beta_Tree_Chao<-einser_matrix-rare_octaves_beta_tree_matrix

#save as csv
write.csv(rare_octaves_Beta_Tree_Chao,"D:\\USB-Laufwerk\\PhD_21.12.2019\\Xylobiont beetles\\2020\\Beta diversity paper\\Paper 2 Nov20\\rare_octaves_Beta_Tree_Chao.csv", row.names = FALSE)


### Replacement versus richness differences

#red list
#common
betaDecompAb<- beta(log(t(beetles_common_strat_1+1)), abund=T)


#datasets with trap numbers as columns
crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")
crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

###do this for all the strata/treespecies!
betaperstratum<-beetles_common_strat
betapertree<-beetles_common_treesp

betaperstratum$Stratum<-crane.env.strat$Stratum[match(betaperstratum$Group.1,crane.env.strat$Fallennummer)]
betapertree$Tree<-crane.env.tree$TreeSp[match(betapertree$Group.1,crane.env.tree$Fallennummer)]

#create abundance matrix equivalent to BCI - for near ground
beetles_ng<-dcast( betaperstratum[which(betaperstratum$Stratum=="Understorey"),], Group.1~Group.2, length)
rownames(  beetles_ng) <-   beetles_ng[,1]
beetles_ng<-   beetles_ng[,-1]
#lower canopy
beetles_lc<-dcast( betaperstratum[which(betaperstratum$Stratum=="Lower"),], Group.1~Group.2, length)
rownames(  beetles_lc) <-   beetles_lc[,1]
beetles_lc<-   beetles_lc[,-1]
#upper canopy
beetles_uc<-dcast( betaperstratum[which(betaperstratum$Stratum=="Upper"),], Group.1~Group.2, length)
rownames(  beetles_uc) <-   beetles_uc[,1]
beetles_uc<-   beetles_uc[,-1]

#create abundance matrix equivalent to BCI - for Tilia
beetles_Til<-dcast( betapertree[which(betapertree$Tree=="Tilia"),], Group.1~Group.2, length)
rownames(  beetles_Til) <-   beetles_Til[,1]
beetles_Til<-   beetles_Til[,-1]

#create abundance matrix equivalent to BCI - for Fraxinus
beetles_Frax<-dcast( betapertree[which(betapertree$Tree=="Fraxinus"),], Group.1~Group.2, length)
rownames(  beetles_Frax) <-   beetles_Frax[,1]
beetles_Frax<-   beetles_Frax[,-1]

#create abundance matrix equivalent to BCI - for Quercus
beetles_Quer<-dcast( betapertree[which(betapertree$Tree=="Quercus"),], Group.1~Group.2, length)
rownames(  beetles_Quer) <-   beetles_Quer[,1]
beetles_Quer<-   beetles_Quer[,-1]

#beta composition calculation per subset

betaDecompNG<- beta(log(beetles_ng+1), abund=T)
betaDecompLC<- beta(log(beetles_lc+1), abund=T)
betaDecompUC<- beta(log(beetles_uc+1), abund=T)

betaDecompTil<- beta(log(beetles_Til+1), abund=T)
betaDecompFrax<- beta(log(beetles_Frax+1), abund=T)
betaDecompQuer<- beta(log(beetles_Quer+1), abund=T)

beta_subset<-c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
beta_repl<-c(mean(betaDecompAb$Brepl), mean(betaDecompNG$Brepl), mean(betaDecompLC$Brepl), mean(betaDecompUC$Brepl), mean(betaDecompTil$Brepl), mean(betaDecompFrax$Brepl), mean(betaDecompQuer$Brepl))
beta_rich<-c(mean(betaDecompAb$Brich), mean(betaDecompNG$Brich), mean(betaDecompLC$Brich), mean(betaDecompUC$Brich), mean(betaDecompTil$Brich), mean(betaDecompFrax$Brich), mean(betaDecompQuer$Brich))

Beta.subsets.data<-data.frame(beta_subset, beta_repl, beta_rich)

#red list
#RARE

rare_betaDecompAb<- beta(log(t(beetles_rare_strat_1+1)), abund=T)

#datasets with trap numbers as columns
crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")
crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

###do this for all the strata/treespecies!
rarebetaperstratum<-beetles_rare_strat
rarebetapertree<-beetles_rare_treesp

rarebetaperstratum$Stratum<-crane.env.strat$Stratum[match(rarebetaperstratum$Group.1,crane.env.strat$Fallennummer)]
rarebetapertree$Tree<-crane.env.tree$TreeSp[match(rarebetapertree$Group.1,crane.env.tree$Fallennummer)]

#create abundance matrix equivalent to BCI - for near ground
rare_beetles_ng<-dcast( rarebetaperstratum[which(rarebetaperstratum$Stratum=="Understorey"),], Group.1~Group.2, length)
rownames(  rare_beetles_ng) <-   rare_beetles_ng[,1]
rare_beetles_ng<-   rare_beetles_ng[,-1]
#lower canopy
rare_beetles_lc<-dcast( rarebetaperstratum[which(rarebetaperstratum$Stratum=="Lower"),], Group.1~Group.2, length)
rownames(  rare_beetles_lc) <-   rare_beetles_lc[,1]
rare_beetles_lc<-   rare_beetles_lc[,-1]
#upper canopy
rare_beetles_uc<-dcast( rarebetaperstratum[which(rarebetaperstratum$Stratum=="Upper"),], Group.1~Group.2, length)
rownames(  rare_beetles_uc) <-   rare_beetles_uc[,1]
rare_beetles_uc<-   rare_beetles_uc[,-1]

#create abundance matrix equivalent to BCI - for Tilia
rare_beetles_Til<-dcast( rarebetapertree[which(rarebetapertree$Tree=="Tilia"),], Group.1~Group.2, length)
rownames(  rare_beetles_Til) <-   rare_beetles_Til[,1]
rare_beetles_Til<-   rare_beetles_Til[,-1]

#create abundance matrix equivalent to BCI - for Fraxinus
rare_beetles_Frax<-dcast( rarebetapertree[which(rarebetapertree$Tree=="Fraxinus"),], Group.1~Group.2, length)
rownames(  rare_beetles_Frax) <-   rare_beetles_Frax[,1]
rare_beetles_Frax<-   rare_beetles_Frax[,-1]

#create abundance matrix equivalent to BCI - for Quercus
rare_beetles_Quer<-dcast( rarebetapertree[which(rarebetapertree$Tree=="Quercus"),], Group.1~Group.2, length)
rownames(  rare_beetles_Quer) <-   rare_beetles_Quer[,1]
rare_beetles_Quer<-   rare_beetles_Quer[,-1]

#beta composition calculation per subset

rare_betaDecompNG<- beta(log(rare_beetles_ng+1), abund=T)
rare_betaDecompLC<- beta(log(rare_beetles_lc+1), abund=T)
rare_betaDecompUC<- beta(log(rare_beetles_uc+1), abund=T)

rare_betaDecompTil<- beta(log(rare_beetles_Til+1), abund=T)
rare_betaDecompFrax<- beta(log(rare_beetles_Frax+1), abund=T)
rare_betaDecompQuer<- beta(log(rare_beetles_Quer+1), abund=T)

rare_beta_subset<-c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
rare_beta_repl<-c(mean(rare_betaDecompAb$Brepl), mean(rare_betaDecompNG$Brepl), mean(rare_betaDecompLC$Brepl), mean(rare_betaDecompUC$Brepl), mean(rare_betaDecompTil$Brepl), mean(rare_betaDecompFrax$Brepl), mean(rare_betaDecompQuer$Brepl))
rare_beta_rich<-c(mean(rare_betaDecompAb$Brich), mean(rare_betaDecompNG$Brich), mean(rare_betaDecompLC$Brich), mean(rare_betaDecompUC$Brich), mean(rare_betaDecompTil$Brich), mean(rare_betaDecompFrax$Brich), mean(rare_betaDecompQuer$Brich))

rare.Beta.subsets.data<-data.frame(rare_beta_subset, rare_beta_repl, rare_beta_rich)
 
#red list
#long dataframe creation
common.Beta.subsets.long<-melt(Beta.subsets.data, id.vars=c("beta_subset"))

# Stacked barplot 
positions <- c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
Brepl_brich_contributions_common<-ggplot(common.Beta.subsets.long, aes(fill=variable, y=value, x=beta_subset, width=0.76)) + 
  geom_bar(position="stack", stat="identity") + scale_x_discrete(limits = positions) +
  labs( x ="Subset", y = "Mean total beta diversity") +
  ggtitle("Common beetle species")+
  scale_fill_discrete(labels = c("Replacement", "Richness differences"))+
  theme(legend.title=element_blank(), axis.text.x = element_text(size=7.5))

#long dataframe creation
rare.Beta.subsets.long<-melt(rare.Beta.subsets.data, id.vars=c("rare_beta_subset"))

# Stacked barplot 
positions <- c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
Brepl_brich_contributions_rare<-ggplot(rare.Beta.subsets.long, aes(fill=variable, y=value, x=rare_beta_subset, width=0.76)) + 
  geom_bar(position="stack", stat="identity") + scale_x_discrete(limits = positions) +
  labs( x ="Subset", y = "Mean total beta diversity") +
  ggtitle("Rare beetle species")+
  scale_fill_discrete(labels = c("Replacement", "Richness differences"))+
  theme(legend.title=element_blank(), axis.text.x = element_text(size=7.5))

rare_common_Brepl_brich_contributions<-grid.arrange(Brepl_brich_contributions_common, Brepl_brich_contributions_rare, nrow = 2)
ggsave("rare_common_Brepl_brich_contributions.png", plot = rare_common_Brepl_brich_contributions, width = 25, height = 17, units = "cm")

#graph goes to Appendix
 
#octaves
#common
oct_betaDecompAb<- beta(log(t(beetles_octaves_common_strat_1+1)), abund=T)


#datasets with trap numbers as columns
crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")
crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

###do this for all the strata/treespecies!
oct_betaperstratum<-beetles_octaves_common_strat
oct_betapertree<-beetles_octaves_common_treesp

oct_betaperstratum$Stratum<-crane.env.strat$Stratum[match(oct_betaperstratum$Group.1,crane.env.strat$Fallennummer)]
oct_betapertree$Tree<-crane.env.tree$TreeSp[match(oct_betapertree$Group.1,crane.env.tree$Fallennummer)]

#create abundance matrix equivalent to BCI - for near ground
oct_beetles_ng<-dcast( oct_betaperstratum[which(oct_betaperstratum$Stratum=="Understorey"),], Group.1~Group.2, length)
rownames(  oct_beetles_ng) <-   oct_beetles_ng[,1]
oct_beetles_ng<-   oct_beetles_ng[,-1]
#lower canopy
oct_beetles_lc<-dcast( oct_betaperstratum[which(oct_betaperstratum$Stratum=="Lower"),], Group.1~Group.2, length)
rownames(  oct_beetles_lc) <-   oct_beetles_lc[,1]
oct_beetles_lc<-   oct_beetles_lc[,-1]
#upper canopy
oct_beetles_uc<-dcast( oct_betaperstratum[which(oct_betaperstratum$Stratum=="Upper"),], Group.1~Group.2, length)
rownames(  oct_beetles_uc) <-   oct_beetles_uc[,1]
oct_beetles_uc<-   oct_beetles_uc[,-1]

#create abundance matrix equivalent to BCI - for Tilia
oct_beetles_Til<-dcast( oct_betapertree[which(oct_betapertree$Tree=="Tilia"),], Group.1~Group.2, length)
rownames(  oct_beetles_Til) <-   oct_beetles_Til[,1]
oct_beetles_Til<-   oct_beetles_Til[,-1]

#create abundance matrix equivalent to BCI - for Fraxinus
oct_beetles_Frax<-dcast( oct_betapertree[which(oct_betapertree$Tree=="Fraxinus"),], Group.1~Group.2, length)
rownames(  oct_beetles_Frax) <-   oct_beetles_Frax[,1]
oct_beetles_Frax<-   oct_beetles_Frax[,-1]

#create abundance matrix equivalent to BCI - for Quercus
oct_beetles_Quer<-dcast( oct_betapertree[which(oct_betapertree$Tree=="Quercus"),], Group.1~Group.2, length)
rownames(  oct_beetles_Quer) <-   oct_beetles_Quer[,1]
oct_beetles_Quer<-   oct_beetles_Quer[,-1]

#beta composition calculation per subset

oct_betaDecompNG<- beta(log(oct_beetles_ng+1), abund=T)
oct_betaDecompLC<- beta(log(oct_beetles_lc+1), abund=T)
oct_betaDecompUC<- beta(log(oct_beetles_uc+1), abund=T)

oct_betaDecompTil<- beta(log(oct_beetles_Til+1), abund=T)
oct_betaDecompFrax<- beta(log(oct_beetles_Frax+1), abund=T)
oct_betaDecompQuer<- beta(log(oct_beetles_Quer+1), abund=T)

beta_subset<-c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
oct_beta_repl<-c(mean(oct_betaDecompAb$Brepl), mean(oct_betaDecompNG$Brepl), mean(oct_betaDecompLC$Brepl), mean(oct_betaDecompUC$Brepl), mean(oct_betaDecompTil$Brepl), mean(oct_betaDecompFrax$Brepl), mean(oct_betaDecompQuer$Brepl))
oct_beta_rich<-c(mean(oct_betaDecompAb$Brich), mean(oct_betaDecompNG$Brich), mean(oct_betaDecompLC$Brich), mean(oct_betaDecompUC$Brich), mean(oct_betaDecompTil$Brich), mean(oct_betaDecompFrax$Brich), mean(oct_betaDecompQuer$Brich))

oct_Beta.subsets.data<-data.frame(beta_subset, oct_beta_repl, oct_beta_rich)


#octaves
#RARE

rare_oct_betaDecompAb<- beta(log(t(beetles_octaves_rare_strat_1+1)), abund=T)

#datasets with trap numbers as columns
crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")
crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

###do this for all the strata/treespecies!
rareoct_betaperstratum<-beetles_octaves_rare_strat
rareoct_betapertree<-beetles_octaves_rare_treesp

rareoct_betaperstratum$Stratum<-crane.env.strat$Stratum[match(rareoct_betaperstratum$Group.1,crane.env.strat$Fallennummer)]
rareoct_betapertree$Tree<-crane.env.tree$TreeSp[match(rareoct_betapertree$Group.1,crane.env.tree$Fallennummer)]

#create abundance matrix equivalent to BCI - for near ground
rare_oct_beetles_ng<-dcast( rareoct_betaperstratum[which(rareoct_betaperstratum$Stratum=="Understorey"),], Group.1~Group.2, length)
rownames(  rare_oct_beetles_ng) <-   rare_oct_beetles_ng[,1]
rare_oct_beetles_ng<-   rare_oct_beetles_ng[,-1]
#lower canopy
rare_oct_beetles_lc<-dcast( rareoct_betaperstratum[which(rareoct_betaperstratum$Stratum=="Lower"),], Group.1~Group.2, length)
rownames(  rare_oct_beetles_lc) <-   rare_oct_beetles_lc[,1]
rare_oct_beetles_lc<-   rare_oct_beetles_lc[,-1]
#upper canopy
rare_oct_beetles_uc<-dcast( rareoct_betaperstratum[which(rareoct_betaperstratum$Stratum=="Upper"),], Group.1~Group.2, length)
rownames(  rare_oct_beetles_uc) <-   rare_oct_beetles_uc[,1]
rare_oct_beetles_uc<-   rare_oct_beetles_uc[,-1]

#create abundance matrix equivalent to BCI - for Tilia
rare_oct_beetles_Til<-dcast( rareoct_betapertree[which(rareoct_betapertree$Tree=="Tilia"),], Group.1~Group.2, length)
rownames(  rare_oct_beetles_Til) <-   rare_oct_beetles_Til[,1]
rare_oct_beetles_Til<-   rare_oct_beetles_Til[,-1]

#create abundance matrix equivalent to BCI - for Fraxinus
rare_oct_beetles_Frax<-dcast( rareoct_betapertree[which(rareoct_betapertree$Tree=="Fraxinus"),], Group.1~Group.2, length)
rownames(  rare_oct_beetles_Frax) <-   rare_oct_beetles_Frax[,1]
rare_oct_beetles_Frax<-   rare_oct_beetles_Frax[,-1]

#create abundance matrix equivalent to BCI - for Quercus
rare_oct_beetles_Quer<-dcast( rareoct_betapertree[which(rareoct_betapertree$Tree=="Quercus"),], Group.1~Group.2, length)
rownames(  rare_oct_beetles_Quer) <-   rare_oct_beetles_Quer[,1]
rare_oct_beetles_Quer<-   rare_oct_beetles_Quer[,-1]

#beta composition calculation per subset

rare_oct_betaDecompNG<- beta(log(rare_oct_beetles_ng+1), abund=T)
rare_oct_betaDecompLC<- beta(log(rare_oct_beetles_lc+1), abund=T)
rare_oct_betaDecompUC<- beta(log(rare_oct_beetles_uc+1), abund=T)

rare_oct_betaDecompTil<- beta(log(rare_oct_beetles_Til+1), abund=T)
rare_oct_betaDecompFrax<- beta(log(rare_oct_beetles_Frax+1), abund=T)
rare_oct_betaDecompQuer<- beta(log(rare_oct_beetles_Quer+1), abund=T)

rare_beta_subset<-c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
rare_oct_beta_repl<-c(mean(rare_oct_betaDecompAb$Brepl), mean(rare_oct_betaDecompNG$Brepl), mean(rare_oct_betaDecompLC$Brepl), mean(rare_oct_betaDecompUC$Brepl), mean(rare_oct_betaDecompTil$Brepl), mean(rare_oct_betaDecompFrax$Brepl), mean(rare_oct_betaDecompQuer$Brepl))
rare_oct_beta_rich<-c(mean(rare_oct_betaDecompAb$Brich), mean(rare_oct_betaDecompNG$Brich), mean(rare_oct_betaDecompLC$Brich), mean(rare_oct_betaDecompUC$Brich), mean(rare_oct_betaDecompTil$Brich), mean(rare_oct_betaDecompFrax$Brich), mean(rare_oct_betaDecompQuer$Brich))

rare.oct_Beta.subsets.data<-data.frame(rare_beta_subset, rare_oct_beta_repl, rare_oct_beta_rich)

#octaves
#long dataframe creation
oct_common.Beta.subsets.long<-melt(oct_Beta.subsets.data, id.vars=c("beta_subset"))

# Stacked barplot 
positions <- c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
oct_Brepl_brich_contributions_common<-ggplot(oct_common.Beta.subsets.long, aes(fill=variable, y=value, x=beta_subset, width=0.76)) + 
  geom_bar(position="stack", stat="identity") + scale_x_discrete(limits = positions) +
  labs( x ="Subset", y = "Mean total beta diversity") +
  ggtitle("Octaves - Common beetle species")+
  scale_fill_discrete(labels = c("Replacement", "Richness differences"))+
  theme(legend.title=element_blank(), axis.text.x = element_text(size=7.5))

#long dataframe creation
oct_rare.Beta.subsets.long<-melt(rare.oct_Beta.subsets.data, id.vars=c("rare_beta_subset"))

# Stacked barplot 
positions <- c("All", "Understorey ", " Lower canopy ", " Upper canopy ", "Tilia", "Fraxinus", "Quercus")
oct_Brepl_brich_contributions_rare<-ggplot(oct_rare.Beta.subsets.long, aes(fill=variable, y=value, x=rare_beta_subset, width=0.76)) + 
  geom_bar(position="stack", stat="identity") + scale_x_discrete(limits = positions) +
  labs( x ="Subset", y = "Mean total beta diversity") +
  ggtitle("Octaves - Rare beetle species")+
  scale_fill_discrete(labels = c("Replacement", "Richness differences"))+
  theme(legend.title=element_blank(), axis.text.x = element_text(size=7.5))

oct_rare_common_Brepl_brich_contributions<-grid.arrange(oct_Brepl_brich_contributions_common, oct_Brepl_brich_contributions_rare, nrow = 2)
ggsave("oct_rare_common_Brepl_brich_contributions.png", plot = oct_rare_common_Brepl_brich_contributions, width = 25, height = 17, units = "cm")

#graph goes to Appendix
 
#############################


### Environmental influences

#### Assessment  of stochastic and deterministic processes of trapwise dissimilarity using Raup-Crick models


#red list
#####employ the model########

# the raup-crick function can be found here: https://github.com/nacmarino/Scripts/blob/master/raup_crick.R

source("raup-crick.R")

#######for the strata
#common
rc_input<-t(beetles_common_strat_1)
common_rc_nora<-raup_crick(rc_input,plot_names_in_col1 = FALSE,reps = 1000)

common_rc_mat<-as.matrix(common_rc_nora)

common_rc_long<-melt(common_rc_mat,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(common_rc_mat),]

crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")

common_rc_long$s1<-crane.env.strat$Stratum[match(common_rc_long$V1,crane.env.strat$Fallennummer)]
common_rc_long$s2<-crane.env.strat$Stratum[match(common_rc_long$V2,crane.env.strat$Fallennummer)]
common_rc_long$sameStrat<-with(common_rc_long, s1==s2)

common_rc_long$g1<-with(common_rc_long,s1=="Understorey")
common_rc_long$g2<-with(common_rc_long,s2=="Understorey")
common_rc_long$ground<-with(common_rc_long,g2==g1)

table(common_rc_long$ground)


#boxplot showing the null model for traps within the three different strata
p <- ggplot(common_rc_long[common_rc_long$sameStrat,], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC", "NG", "UC"))+
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
p1 <- ggplot(common_rc_long, aes(x=sameStrat, y=raupCrick, fill=sameStrat)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))

p3 <- ggplot(common_rc_long[which(common_rc_long$sameStrat=="FALSE"&common_rc_long$s2=="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Upper"), labels=c("LC \n - NG", "UC \n - NG"))+
  theme(text = element_text(size=10))

p4 <- ggplot(common_rc_long[which(common_rc_long$sameStrat=="FALSE"&common_rc_long$s2!="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC \n - UC", "NG \n - C", "UC \n - LC"))+
  theme(text = element_text(size=10))

common_RaupCrick_stratum<-arrangeGrob(p, p1, p4, p3, 
                                      nrow = 2)



# Writing to file
ggsave("common_RaupCrick_stratum.png", plot = common_RaupCrick_stratum, width = 20, height = 17, units = "cm")

#rare
rc_input_rare<-t(beetles_rare_strat_1)
rare_rc_nora<-raup_crick(rc_input_rare,plot_names_in_col1 = FALSE,reps = 1000)

rare_rc_mat<-as.matrix(rare_rc_nora)
rare_rc_long<-melt(rare_rc_mat,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(rare_rc_mat),]

crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")

rare_rc_long$s1<-crane.env.strat$Stratum[match(rare_rc_long$V1,crane.env.strat$Fallennummer)]
rare_rc_long$s2<-crane.env.strat$Stratum[match(rare_rc_long$V2,crane.env.strat$Fallennummer)]
rare_rc_long$sameStrat<-with(rare_rc_long, s1==s2)

rare_rc_long$g1<-with(rare_rc_long,s1=="Understorey")
rare_rc_long$g2<-with(rare_rc_long,s2=="Understorey")
rare_rc_long$ground<-with(rare_rc_long,g2==g1)

table(rare_rc_long$ground)


#boxplot showing the null model for traps within the three different strata
p5 <- ggplot(rare_rc_long[rare_rc_long$sameStrat,], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC", "NG", "UC"))+
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
p6 <- ggplot(rare_rc_long, aes(x=sameStrat, y=raupCrick, fill=sameStrat)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))

p7 <- ggplot(rare_rc_long[which(rare_rc_long$sameStrat=="FALSE"&rare_rc_long$s2=="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Upper"), labels=c("LC \n - NG", "UC \n - NG"))+
  theme(text = element_text(size=10))

p8 <- ggplot(rare_rc_long[which(rare_rc_long$sameStrat=="FALSE"&rare_rc_long$s2!="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC \n - UC", "NG \n - C", "UC \n - LC"))+
  theme(text = element_text(size=10))

rare_RaupCrick_stratum<-arrangeGrob(p5, p6, p8, p7, 
                                    nrow = 2)


# Writing to file
ggsave("rare_RaupCrick_stratum.png", plot = rare_RaupCrick_stratum, width = 20, height = 17, units = "cm")


#######for the tree species

#common
crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

rc_input_tree<-t(beetles_common_treesp_1)

common_rc_nora_tree<-raup_crick(rc_input_tree,plot_names_in_col1 = FALSE,reps = 1000)

common_rc_mat_tree<-as.matrix(common_rc_nora_tree)
common_rc_long_tree<-melt(common_rc_mat_tree,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(common_rc_mat_tree),]

common_rc_long_tree$t1<-crane.env.tree$TreeSp[match(common_rc_long_tree$V1,crane.env.tree$Fallennummer)]
common_rc_long_tree$t2<-crane.env.tree$TreeSp[match(common_rc_long_tree$V2,crane.env.tree$Fallennummer)]
common_rc_long_tree$sameTree<-with(common_rc_long_tree, t1==t2)


#boxplot showing the null model for traps within the three different strata
p9 <- ggplot(common_rc_long_tree[common_rc_long_tree$sameTree,], aes(x=t1, y=raupCrick, fill=t1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
p10 <- ggplot(common_rc_long_tree, aes(x=sameTree, y=raupCrick, fill=sameTree)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))

p11 <- ggplot(common_rc_long_tree[which(common_rc_long_tree$sameTree=="FALSE"&common_rc_long_tree$t1=="Tilia"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Quercus"), labels=c("Tilia \n - Fraxinus", "Tilia \n - Quercus"))+
  theme(text = element_text(size=10))

p12 <- ggplot(common_rc_long_tree[which(common_rc_long_tree$sameTree=="FALSE"&common_rc_long_tree$t1=="Fraxinus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Quercus", "Tilia"), labels=c("Fraxinus \n - Quercus", "Fraxinus \n - Tilia"))+
  theme(text = element_text(size=10))

p13 <- ggplot(common_rc_long_tree[which(common_rc_long_tree$sameTree=="FALSE"&common_rc_long_tree$t1=="Quercus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Tilia"), labels=c("Quercus \n - Fraxinus", "Quercus \n - Tilia"))+
  theme(text = element_text(size=10))


#save
common_RaupCrick_tree <- arrangeGrob(p9, p10, p11, p12, p13 ,nrow = 3) 


# Writing to file
ggsave("common_RaupCrick_tree.png", plot = common_RaupCrick_tree, width = 20, height = 17, units = "cm")


#rare

crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

rc_input_rare_tree<-t(beetles_rare_treesp_1)

rare_rc_nora_tree<-raup_crick(rc_input_rare_tree,plot_names_in_col1 = FALSE,reps = 1000)

rare_rc_mat_tree<-as.matrix(rare_rc_nora_tree)
rare_rc_long_tree<-melt(rare_rc_mat_tree,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(rare_rc_mat_tree),]

rare_rc_long_tree$t1<-crane.env.tree$TreeSp[match(rare_rc_long_tree$V1,crane.env.tree$Fallennummer)]
rare_rc_long_tree$t2<-crane.env.tree$TreeSp[match(rare_rc_long_tree$V2,crane.env.tree$Fallennummer)]
rare_rc_long_tree$sameTree<-with(rare_rc_long_tree, t1==t2)


#boxplot showing the null model for traps within the three different strata
p14 <- ggplot(rare_rc_long_tree[rare_rc_long_tree$sameTree,], aes(x=t1, y=raupCrick, fill=t1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
p15 <- ggplot(rare_rc_long_tree, aes(x=sameTree, y=raupCrick, fill=sameTree)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))

p16 <- ggplot(rare_rc_long_tree[which(rare_rc_long_tree$sameTree=="FALSE"&rare_rc_long_tree$t1=="Tilia"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Quercus"), labels=c("Tilia \n - Fraxinus", "Tilia \n - Quercus"))+
  theme(text = element_text(size=10))

p17 <- ggplot(rare_rc_long_tree[which(rare_rc_long_tree$sameTree=="FALSE"&rare_rc_long_tree$t1=="Fraxinus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Quercus", "Tilia"), labels=c("Fraxinus \n - Quercus", "Fraxinus \n - Tilia"))+
  theme(text = element_text(size=10))

p18 <- ggplot(rare_rc_long_tree[which(rare_rc_long_tree$sameTree=="FALSE"&rare_rc_long_tree$t1=="Quercus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Tilia"), labels=c("Quercus \n - Fraxinus", "Quercus \n - Tilia"))+
  theme(text = element_text(size=10))


#save
rare_RaupCrick_tree <- arrangeGrob(p14, p15, p16, p17, p18 ,nrow = 3) 


# Writing to file
ggsave("rare_RaupCrick_tree.png", plot = rare_RaupCrick_tree, width = 20, height = 17, units = "cm")

# octaves

#####employ the model########

source("raup-crick.R")

#######for the strata
#common
oct_rc_input<-t(beetles_common_strat_1)
oct_oct_common_rc_nora<-raup_crick(oct_rc_input,plot_names_in_col1 = FALSE,reps = 1000)

oct_common_rc_mat<-as.matrix(oct_oct_common_rc_nora)

oct_common_rc_long<-melt(oct_common_rc_mat,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(oct_common_rc_mat),]

crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")

oct_common_rc_long$s1<-crane.env.strat$Stratum[match(oct_common_rc_long$V1,crane.env.strat$Fallennummer)]
oct_common_rc_long$s2<-crane.env.strat$Stratum[match(oct_common_rc_long$V2,crane.env.strat$Fallennummer)]
oct_common_rc_long$sameStrat<-with(oct_common_rc_long, s1==s2)

oct_common_rc_long$g1<-with(oct_common_rc_long,s1=="Understorey")
oct_common_rc_long$g2<-with(oct_common_rc_long,s2=="Understorey")
oct_common_rc_long$ground<-with(oct_common_rc_long,g2==g1)

table(oct_common_rc_long$ground)


#boxplot showing the null model for traps within the three different strata
q <- ggplot(oct_common_rc_long[oct_common_rc_long$sameStrat,], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC", "NG", "UC"))+
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
q1 <- ggplot(oct_common_rc_long, aes(x=sameStrat, y=raupCrick, fill=sameStrat)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))

q3 <- ggplot(oct_common_rc_long[which(oct_common_rc_long$sameStrat=="FALSE"&oct_common_rc_long$s2=="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Upper"), labels=c("LC \n - NG", "UC \n - NG"))+
  theme(text = element_text(size=10))

q4 <- ggplot(oct_common_rc_long[which(oct_common_rc_long$sameStrat=="FALSE"&oct_common_rc_long$s2!="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC \n - UC", "NG \n - C", "UC \n - LC"))+
  theme(text = element_text(size=10))

oct_common_RaupCrick_stratum<-arrangeGrob(q, q1, q4, q3, 
                                          nrow = 2)


# Writing to file
ggsave("oct_common_RaupCrick_stratum.png", plot = oct_common_RaupCrick_stratum, width = 20, height = 17, units = "cm")

#rare
oct_rc_input_rare<-t(beetles_octaves_rare_strat_1)
oct_rare_rc_nora<-raup_crick(oct_rc_input_rare,plot_names_in_col1 = FALSE,reps = 1000)

oct_rare_rc_mat<-as.matrix(oct_rare_rc_nora)
oct_rare_rc_long<-melt(oct_rare_rc_mat,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(oct_rare_rc_mat),]

crane.env.strat<-read.csv("crane.env.2017.csv", sep=";")

oct_rare_rc_long$s1<-crane.env.strat$Stratum[match(oct_rare_rc_long$V1,crane.env.strat$Fallennummer)]
oct_rare_rc_long$s2<-crane.env.strat$Stratum[match(oct_rare_rc_long$V2,crane.env.strat$Fallennummer)]
oct_rare_rc_long$sameStrat<-with(oct_rare_rc_long, s1==s2)

oct_rare_rc_long$g1<-with(oct_rare_rc_long,s1=="Understorey")
oct_rare_rc_long$g2<-with(oct_rare_rc_long,s2=="Understorey")
oct_rare_rc_long$ground<-with(oct_rare_rc_long,g2==g1)

table(oct_rare_rc_long$ground)


#boxplot showing the null model for traps within the three different strata
q5 <- ggplot(oct_rare_rc_long[oct_rare_rc_long$sameStrat,], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC", "NG", "UC"))+
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
q6 <- ggplot(oct_rare_rc_long, aes(x=sameStrat, y=raupCrick, fill=sameStrat)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))


q7 <- ggplot(oct_rare_rc_long[which(oct_rare_rc_long$sameStrat=="FALSE"&oct_rare_rc_long$s2=="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Upper"), labels=c("LC \n - NG", "UC \n - NG"))+
  theme(text = element_text(size=10))

q8 <- ggplot(oct_rare_rc_long[which(oct_rare_rc_long$sameStrat=="FALSE"&oct_rare_rc_long$s2!="Understorey"),], aes(x=s1, y=raupCrick, fill=s1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Lower", "Understorey", "Upper"), labels=c("LC \n - UC", "NG \n - C", "UC \n - LC"))+
  theme(text = element_text(size=10))

oct_rare_RaupCrick_stratum<-arrangeGrob(q5, q6, q8, q7, 
                                        nrow = 2)


# Writing to file
ggsave("oct_rare_RaupCrick_stratum.png", plot = oct_rare_RaupCrick_stratum, width = 20, height = 17, units = "cm")


#######for the tree species

#common
crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

oct_rc_input_tree<-t(beetles_octaves_common_treesp_1)

oct_oct_common_rc_nora_tree<-raup_crick(oct_rc_input_tree,plot_names_in_col1 = FALSE,reps = 1000)

oct_common_rc_mat_tree<-as.matrix(oct_oct_common_rc_nora_tree)
oct_common_rc_long_tree<-melt(oct_common_rc_mat_tree,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(oct_common_rc_mat_tree),]

oct_common_rc_long_tree$t1<-crane.env.tree$TreeSp[match(oct_common_rc_long_tree$V1,crane.env.tree$Fallennummer)]
oct_common_rc_long_tree$t2<-crane.env.tree$TreeSp[match(oct_common_rc_long_tree$V2,crane.env.tree$Fallennummer)]
oct_common_rc_long_tree$sameTree<-with(oct_common_rc_long_tree, t1==t2)


#boxplot showing the null model for traps within the three different strata
q9 <- ggplot(oct_common_rc_long_tree[oct_common_rc_long_tree$sameTree,], aes(x=t1, y=raupCrick, fill=t1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
q10 <- ggplot(oct_common_rc_long_tree, aes(x=sameTree, y=raupCrick, fill=sameTree)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))

q11 <- ggplot(oct_common_rc_long_tree[which(oct_common_rc_long_tree$sameTree=="FALSE"&oct_common_rc_long_tree$t1=="Tilia"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Quercus"), labels=c("Tilia \n - Fraxinus", "Tilia \n - Quercus"))+
  theme(text = element_text(size=10))

q12 <- ggplot(oct_common_rc_long_tree[which(oct_common_rc_long_tree$sameTree=="FALSE"&oct_common_rc_long_tree$t1=="Fraxinus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Quercus", "Tilia"), labels=c("Fraxinus \n - Quercus", "Fraxinus \n - Tilia"))+
  theme(text = element_text(size=10))

q13 <- ggplot(oct_common_rc_long_tree[which(oct_common_rc_long_tree$sameTree=="FALSE"&oct_common_rc_long_tree$t1=="Quercus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Tilia"), labels=c("Quercus \n - Fraxinus", "Quercus \n - Tilia"))+
  theme(text = element_text(size=10))


#save
oct_common_RaupCrick_tree <- arrangeGrob(q9, q10, q11, q12, q13 ,nrow = 3) 


# Writing to file
ggsave("oct_common_RaupCrick_tree.png", plot = oct_common_RaupCrick_tree, width = 20, height = 17, units = "cm")


#rare

crane.env.tree <- read.csv(file="crane.env.treesp.2017.csv", header=TRUE, sep=";")
crane.env.tree$Fallennummer<-as.numeric(crane.env.tree$Fallennummer)

oct_rc_input_rare_tree<-t(beetles_octaves_rare_treesp_1)

oct_rare_rc_nora_tree<-raup_crick(oct_rc_input_rare_tree,plot_names_in_col1 = FALSE,reps = 1000)

oct_rare_rc_mat_tree<-as.matrix(oct_rare_rc_nora_tree)
oct_rare_rc_long_tree<-melt(oct_rare_rc_mat_tree,varnames = c("V1","V2"),value.name = "raupCrick")[lower.tri(oct_rare_rc_mat_tree),]

oct_rare_rc_long_tree$t1<-crane.env.tree$TreeSp[match(oct_rare_rc_long_tree$V1,crane.env.tree$Fallennummer)]
oct_rare_rc_long_tree$t2<-crane.env.tree$TreeSp[match(oct_rare_rc_long_tree$V2,crane.env.tree$Fallennummer)]
oct_rare_rc_long_tree$sameTree<-with(oct_rare_rc_long_tree, t1==t2)


#boxplot showing the null model for traps within the three different strata
q14 <- ggplot(oct_rare_rc_long_tree[oct_rare_rc_long_tree$sameTree,], aes(x=t1, y=raupCrick, fill=t1)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none" ) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  theme(text = element_text(size=10)) 

#boxplot showing the null model for traps within a stratum vs traps between strata
q15 <- ggplot(oct_rare_rc_long_tree, aes(x=sameTree, y=raupCrick, fill=sameTree)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("FALSE", "TRUE"), labels=c("Between", "Within"))+
  theme(text = element_text(size=10))

q16 <- ggplot(oct_rare_rc_long_tree[which(oct_rare_rc_long_tree$sameTree=="FALSE"&oct_rare_rc_long_tree$t1=="Tilia"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Quercus"), labels=c("Tilia \n - Fraxinus", "Tilia \n - Quercus"))+
  theme(text = element_text(size=10))

q17 <- ggplot(oct_rare_rc_long_tree[which(oct_rare_rc_long_tree$sameTree=="FALSE"&oct_rare_rc_long_tree$t1=="Fraxinus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Quercus", "Tilia"), labels=c("Fraxinus \n - Quercus", "Fraxinus \n - Tilia"))+
  theme(text = element_text(size=10))

q18 <- ggplot(oct_rare_rc_long_tree[which(oct_rare_rc_long_tree$sameTree=="FALSE"&oct_rare_rc_long_tree$t1=="Quercus"),], aes(x=t2, y=raupCrick, fill=t2)) + 
  geom_boxplot() +
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),text = element_text(size=10), legend.position = "none") +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
  scale_x_discrete(breaks=c("Fraxinus", "Tilia"), labels=c("Quercus \n - Fraxinus", "Quercus \n - Tilia"))+
  theme(text = element_text(size=10))


#save
oct_rare_RaupCrick_tree <- arrangeGrob(q14, q15, q16, q17, q18 ,nrow = 3) 


# Writing to file
ggsave("oct_rare_RaupCrick_tree.png", plot = oct_rare_RaupCrick_tree, width = 20, height = 17, units = "cm")

######################

#### Dissimilarity between strata and tree species

#red list

#### NMDS plots 

#per strata
strat_common_NMDS<-metaMDS(common_Beta_Strata_Chao)

####per tree sp
#NMDS
treesp_common_NMDS<-metaMDS(common_Beta_Tree_Chao)

####rare 
#per strata
strat_rare_NMDS<-metaMDS(rare_Beta_Strata_Chao)

####per tree sp
#NMDS
treesp_rare_NMDS<-metaMDS(rare_Beta_Tree_Chao)

####plot all tree sp nmds
# Open a png file
png("tree_nmds.png", width = 1350, height = 1000) 

####plot all strata nmds
par(mfrow=c(1,2))

tree_common<-plot(treesp_common_NMDS, disp = "sites", type = "n", main = "a)", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#add the factors as polygons
ordihull(treesp_common_NMDS,groups=crane.env.tree$TreeSp,draw="polygon",label=F)
with(crane.env.tree, ordihull(treesp_common_NMDS, TreeSp, show.groups = "Tilia", draw = "polygon", col="palegreen3"))
with(crane.env.tree, ordihull(treesp_common_NMDS, TreeSp, show.groups = "Fraxinus", draw = "polygon", col="skyblue"))
with(crane.env.tree, ordihull(treesp_common_NMDS, TreeSp, show.groups = "Quercus", draw = "polygon", col="orange"))
legend("topright", legend = c('Tilia', 'Fraxinus', 'Quercus'), pch = 15, 
       pt.cex = 3, cex = 1.7, title = NULL, col = c('palegreen3', 'skyblue', 'orange'))

tree_rare<-plot(treesp_rare_NMDS, disp = "sites", type = "n", main = "b)",cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#add the factors as polygons
ordihull(treesp_rare_NMDS,groups=crane.env.tree$TreeSp,draw="polygon",label=F)
with(crane.env.tree, ordihull(treesp_rare_NMDS, TreeSp, show.groups = "Tilia", draw = "polygon", col="palegreen3"))
with(crane.env.tree, ordihull(treesp_rare_NMDS, TreeSp, show.groups = "Fraxinus", draw = "polygon", col="skyblue"))
with(crane.env.tree, ordihull(treesp_rare_NMDS, TreeSp, show.groups = "Quercus", draw = "polygon", col="orange"))
legend("topright", legend = c('Tilia', 'Fraxinus', 'Quercus'), pch = 15, 
       pt.cex = 3, cex = 1.7, title = NULL, col = c('palegreen3', 'skyblue', 'orange'))

# Close the png file
dev.off()

# Open a png file
png("strata_nmds.png", width = 1350, height = 1000) 

####plot all strata nmds
par(mfrow=c(1,2))

plot(strat_common_NMDS, disp = "sites", type = "n", main = "a)",cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5) 
#add environmental factors
#add the factors as polygons
ordihull(strat_common_NMDS, groups=crane.env.strat$Stratum, draw="polygon",label=F)
with(crane.env.strat, ordihull(strat_common_NMDS, Stratum, show.groups = "Lower", draw = "polygon", col="palegreen3"))
with(crane.env.strat, ordihull(strat_common_NMDS, Stratum, show.groups = "Upper", draw = "polygon", col="skyblue"))
with(crane.env.strat, ordihull(strat_common_NMDS, Stratum, show.groups = "Understorey", draw = "polygon", col="orange"))
legend("topright", legend = c('LC', 'UC', 'NG'), pch = 15, pt.cex = 3, cex = 1.7, 
       title = NULL,
       col = c('palegreen3', 'skyblue', 'orange'))

#rare
plot(strat_rare_NMDS, disp = "sites", type = "n",main = "b)",cex.main = 1.5,  cex.lab = 1.5, cex.axis = 1.5) 
#add environmental factors
#add the factors as polygons
ordihull(strat_rare_NMDS, groups=crane.env.strat$Stratum, draw="polygon",label=F)
with(crane.env.strat, ordihull(strat_rare_NMDS, Stratum, show.groups = "Lower", draw = "polygon", col="palegreen3"))
with(crane.env.strat, ordihull(strat_rare_NMDS, Stratum, show.groups = "Upper", draw = "polygon", col="skyblue"))
with(crane.env.strat, ordihull(strat_rare_NMDS, Stratum, show.groups = "Understorey", draw = "polygon", col="orange"))
legend("topright", legend = c('LC', 'UC', 'NG'), pch = 15, pt.cex = 3, cex = 1.7, 
       title = NULL,
       col = c('palegreen3', 'skyblue', 'orange'))

# Close the png file
dev.off() 

#octaves

#### NMDS plots 

#per strata
oct_strat_common_NMDS<-metaMDS(common_octaves_Beta_Strata_Chao)

####per tree sp
#NMDS
oct_tree_sp_common_NMDS<-metaMDS(common_octaves_Beta_Tree_Chao)

####rare 
#per strata
oct_strat_rare_NMDS<-metaMDS(rare_octaves_Beta_Strata_Chao)

####per tree sp
#NMDS
oct_tree_sp_rare_NMDS<-metaMDS(rare_octaves_Beta_Tree_Chao)

####plot all tree sp nmds
# Open a png file
png("oct_tree_nmds.png", width = 1350, height = 1000) 

####plot all trees nmds
par(mfrow=c(1,2))

oct_tree_common<-plot(oct_tree_sp_common_NMDS, disp = "sites", type = "n", main = "a)",cex.main = 1.5,  cex.lab = 1.5, cex.axis = 1.5)
#add the factors as polygons
ordihull(oct_tree_sp_common_NMDS,groups=crane.env.tree$TreeSp,draw="polygon",label=F)
with(crane.env.tree, ordihull(oct_tree_sp_common_NMDS, TreeSp, show.groups = "Tilia", draw = "polygon", col="palegreen3"))
with(crane.env.tree, ordihull(oct_tree_sp_common_NMDS, TreeSp, show.groups = "Fraxinus", draw = "polygon", col="skyblue"))
with(crane.env.tree, ordihull(oct_tree_sp_common_NMDS, TreeSp, show.groups = "Quercus", draw = "polygon", col="orange"))
legend("topright", legend = c('Tilia', 'Fraxinus', 'Quercus'), pch = 15, 
       pt.cex = 3, cex = 1.7, title = NULL, col = c('palegreen3', 'skyblue', 'orange'))

oct_tree_rare<-plot(oct_tree_sp_rare_NMDS, disp = "sites", type = "n",main = "b)",cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#add the factors as polygons
ordihull(oct_tree_sp_rare_NMDS,groups=crane.env.tree$TreeSp,draw="polygon",label=F)
with(crane.env.tree, ordihull(oct_tree_sp_rare_NMDS, TreeSp, show.groups = "Tilia", draw = "polygon", col="palegreen3"))
with(crane.env.tree, ordihull(oct_tree_sp_rare_NMDS, TreeSp, show.groups = "Fraxinus", draw = "polygon", col="skyblue"))
with(crane.env.tree, ordihull(oct_tree_sp_rare_NMDS, TreeSp, show.groups = "Quercus", draw = "polygon", col="orange"))
legend("topright", legend = c('Tilia', 'Fraxinus', 'Quercus'), pch = 15, 
       pt.cex = 3, cex = 1.7, title = NULL, col = c('palegreen3', 'skyblue', 'orange'))

# Close the png file
dev.off()

# Open a png file
png("oct_strata_nmds.png", width = 1350, height = 1000) 

####plot all strata nmds
par(mfrow=c(1,2))

plot(oct_strat_common_NMDS, disp = "sites", type = "n",main = "a)",cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5) 
#add environmental factors
#add the factors as polygons
ordihull(oct_strat_common_NMDS, groups=crane.env.strat$Stratum, draw="polygon",label=F)
with(crane.env.strat, ordihull(oct_strat_common_NMDS, Stratum, show.groups = "Lower", draw = "polygon", col="palegreen3"))
with(crane.env.strat, ordihull(oct_strat_common_NMDS, Stratum, show.groups = "Upper", draw = "polygon", col="skyblue"))
with(crane.env.strat, ordihull(oct_strat_common_NMDS, Stratum, show.groups = "Understorey", draw = "polygon", col="orange"))
legend("topright", legend = c('LC', 'UC', 'NG'), pch = 15, pt.cex = 3, cex = 1.7, 
       title = NULL,
       col = c('palegreen3', 'skyblue', 'orange'))

#rare
plot(oct_strat_rare_NMDS, disp = "sites", type = "n",main = "b)",cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5) 
#add environmental factors
#add the factors as polygons
ordihull(oct_strat_rare_NMDS, groups=crane.env.strat$Stratum, draw="polygon",label=F)
with(crane.env.strat, ordihull(oct_strat_rare_NMDS, Stratum, show.groups = "Lower", draw = "polygon", col="palegreen3"))
with(crane.env.strat, ordihull(oct_strat_rare_NMDS, Stratum, show.groups = "Upper", draw = "polygon", col="skyblue"))
with(crane.env.strat, ordihull(oct_strat_rare_NMDS, Stratum, show.groups = "Understorey", draw = "polygon", col="orange"))
legend("topright", legend = c('LC', 'UC', 'NG'), pch = 15, pt.cex = 3, cex = 1.7, 
       title = NULL,
       col = c('palegreen3', 'skyblue', 'orange'))

# Close the png file
dev.off() 

###Analysis of variance using distance matrix

adonis2(common_Beta_Strata_Chao ~ Stratum, data = crane.env.strat)
adonis2(rare_Beta_Strata_Chao ~ Stratum, data = crane.env.strat)
adonis2(common_Beta_Tree_Chao ~ TreeSp, data = crane.env.tree)
adonis2(rare_Beta_Tree_Chao ~ TreeSp, data = crane.env.tree)

adonis_results<-read.csv2("adonis_results_23.06.2021.csv", sep=";")
adonis_results$SumOfSqs<-as.numeric(adonis_results$SumOfSqs)
adonis_results$R2<-as.numeric(adonis_results$R2)
adonis_results$F<-as.numeric(adonis_results$F)
adonis_results$Pr..F.<-as.numeric(adonis_results$Pr..F.)
adonis_results$Df<-as.numeric(adonis_results$Df)

adonis_results$SumOfSqs<-round(adonis_results$SumOfSqs, 2)
adonis_results$R2<-round(adonis_results$R2, 2)
adonis_results$F<-round(adonis_results$F, 2)


adonis.res<-flextable(adonis_results)
adonis.res <- theme_box(adonis.res)
adonis.res <- merge_v(adonis.res, "X")
adonis.res<-set_header_labels(adonis.res, X = "Rarity", X.1= "Factor", Df= "Degree of freedom", SumOfSqs = "Sum of squares", MeanSqs = "Mean squares", F= "F statistics", R2 = "R-squared", 'Pr..F.' = "P value")
adonis.res<- align(adonis.res, align = "center")
adonis.res<-set_caption(adonis.res, "Table 1: PERMANOVA quantifying the impact of different environmental sources of variation \n on the dissimilarities between traps
                        ")
adonis.res<-autofit(adonis.res)
adonis.res

### if you have questions about these analyses please contact nora.haack@idiv.de

