setwd("C:/Users/PC/Documents/Beetle Case Study/PPA data")

packages <- c("ggpubr", "knitr", "ggplot2", "broom", "devtools", "reshape2", "vegan", "dendextend", "flextable", "officer","BiodiversityR", "dplyr", "tidyr", "BAT", "gridExtra", "spaa", "tinytex", "SpadeR")

install.packages(setdiff(packages, rownames(installed.packages()))) 
for (pkg in packages) {
  library(pkg, character.only = TRUE)
}

library("ggpubr")
library("ggplot2")
install.packages("scaales")
library("scales")

#load data
dry_13 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_13_dry.csv")
int_13 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_13_intermediate.csv")
moist_13 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_13_moist.csv")

dry_16 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_16_dry.csv")
int_16 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_16_intermediate.csv")
moist_16 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_16_moist.csv")

dry_no_increase <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_dry_init_dry.csv")
dry_increase_half <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_middle_init_dry.csv")
dry_increase_1 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_moist_init_dry.csv")
int_no_increase <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_middle_init_intermediate.csv")
int_increase_half <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_moist_init_intermediate.csv")
int_increase_1 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_wet_init_intermediate.csv")
moist_no_increase <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_moist_init_moist.csv")
moist_increase_half <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_wet_init_moist.csv")
moist_increase_1 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/Simulations/Results_combined_census_v_wet_init_moist.csv")

tree_7 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/tree_species7.csv")
tree_8 <- read.csv("C:/Users/PC/Documents/Beetles Case Study/PPA data/tree_species8.csv")

#join tree name
dry_13 <- merge(x = dry_13, y = tree_7, by.x ="sp", by.y = "Sp")
int_13 <- merge(x = int_13, y = tree_7, by.x ="sp", by.y = "Sp")
moist_13 <- merge(x = moist_13, y = tree_7, by.x ="sp", by.y = "Sp")

dry_16 <- merge(x = dry_16, y = tree_7, by.x ="sp", by.y = "Sp")
int_16 <- merge(x = int_16, y = tree_7, by.x ="sp", by.y = "Sp")
moist_16 <- merge(x = moist_16, y = tree_7, by.x ="sp", by.y = "Sp")

dry_no_increase <- merge(x = dry_no_increase, y = tree_8, by.x ="sp", by.y = "Sp")
dry_increase_half <- merge(x = dry_increase_half, y = tree_8, by.x ="sp", by.y = "Sp")
dry_increase_1 <- merge(x = dry_increase_1, y = tree_8, by.x ="sp", by.y = "Sp")
int_no_increase <- merge(x = int_no_increase, y = tree_8, by.x ="sp", by.y = "Sp")
int_increase_half <- merge(x = int_increase_half, y = tree_8, by.x ="sp", by.y = "Sp")
int_increase_1 <- merge(x = int_increase_1, y = tree_8, by.x ="sp", by.y = "Sp")
moist_no_increase <- merge(x = moist_no_increase, y = tree_8, by.x ="sp", by.y = "Sp") 
moist_increase_half <- merge(x = moist_increase_half, y = tree_8, by.x ="sp", by.y = "Sp")
moist_increase_1 <- merge(x = moist_increase_1, y = tree_8, by.x ="sp", by.y = "Sp")

#summary and test

summary <- dry_16%>%
  group_by(sp, time, cl) %>%
  summarise(n_entry = n_distinct(dbh)) 

test <- dry_13[dry_13$sp == 1 & dry_13$time == 0,] 

test2 <- moist_increase_1 %>%
  distinct(sp)

#Note: in moist, int, dry of 13 and 16, each has 7 tree species. In combined_data, each has 8 tree species.

#calculation basal area
#2013-2020 data dry
dry_13_cal <- dry_13[dry_13$Name != "Carpinus betulus" & dry_13$Name != "Acer platanoides", ] %>%
  select(dbh, Name, time) %>%
  mutate(
    ba = (dbh/200)^2
  ) %>%
  group_by(time, Name) %>%
  summarise(sum_ba = sum(ba)) 

dry_13_plot <- ggplot(dry_13_cal) + 
  geom_line(aes(x = time, y = sum_ba, color = Name), linewidth = 2) + 
  scale_color_manual(values = c("Acer" = "#228B22",
                                "Fraxinus" = "#00FF00",
                                "Quercus robur" = "#FFA500",
                                "Ulmus" = "blue",
                                "Tilia" = "red")) +
  scale_y_continuous(
    limits = c(0, 12),
    labels = number_format(accuracy = 0.1)) +
  xlab("Time (years)") + ylab("Basal area (m²/ha)") +
  ggtitle("Dry (2.3 m DTG)") +
  guides(color = guide_legend(override.aes = list(linewidth = 8))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = "right", 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 18), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 18))

#2013-2020 data int
int_13_cal <- int_13[int_13$Name != "Carpinus betulus" & int_13$Name != "Acer platanoides", ] %>%
  select(dbh, Name, time) %>%
  mutate(
    ba = (dbh/200)^2
  ) %>%
  group_by(time, Name) %>%
  summarise(sum_ba = sum(ba)) 

int_13_plot <- ggplot(int_13_cal) + 
  geom_line(aes(x = time, y = sum_ba, color = Name), linewidth = 2) + 
  scale_color_manual(values = c("Acer" = "#228B22",
                                "Fraxinus" = "#00FF00",
                                "Quercus robur" = "#FFA500",
                                "Ulmus" = "blue",
                                "Tilia" = "red")) +
  scale_y_continuous(
    limits = c(0, 12),
    labels = number_format(accuracy = 0.1)) +
  xlab("Time (years)") + ylab("Basal area (m²/ha)") +
  ggtitle("Intermediate (1.8 m DTG)") +
  guides(color = guide_legend(override.aes = list(linewidth = 8))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = "right", 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 18), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 18))

#2013-2020 data moist
moist_13_cal <- moist_13[moist_13$Name != "Carpinus betulus" & moist_13$Name != "Acer platanoides", ] %>%
  select(dbh, Name, time) %>%
  mutate(
    ba = (dbh/200)^2
  ) %>%
  group_by(time, Name) %>%
  summarise(sum_ba = sum(ba)) 

moist_13_plot <- ggplot(moist_13_cal) + 
  geom_line(aes(x = time, y = sum_ba, color = Name), linewidth = 2) + 
  scale_color_manual(values = c("Acer" = "#228B22",
                                "Fraxinus" = "#00FF00",
                                "Quercus robur" = "#FFA500",
                                "Ulmus" = "blue",
                                "Tilia" = "red")) + 
  scale_y_continuous(
    limits = c(0, 12),
    labels = number_format(accuracy = 0.1)) +
  xlab("Time (years)") + ylab("Basal area (m²/ha)") +
  ggtitle("Moist (1.5 m DTG)") +
  guides(color = guide_legend(override.aes = list(linewidth = 8))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = "right", 
        axis.title = element_text(size = 22), 
        axis.text = element_text(size = 18), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 18))

plot_2013 <- ggarrange(dry_13_plot, int_13_plot, moist_13_plot,
                       common.legend = TRUE, legend = "bottom",
                       ncol = 3, nrow = 1)

ggsave("plot_2013.png", width = 30, height = 10)







  


