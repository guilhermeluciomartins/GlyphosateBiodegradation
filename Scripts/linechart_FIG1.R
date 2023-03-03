### Multiple line chart ### Guilherme (28/09/2022)

### Definir diretorio de trabalho ###
getwd()
setwd("C:/Users/gui_l/OneDrive/Doutorado/Colaborações/Experimento glifosato")
path<-"C:/Users/gui_l/OneDrive/Doutorado/Colaborações/Experimento glifosato"
list.files(path)

# Install packages
install.packages("reshape")
install.packages("reshape2")

# Load packages
library(reshape)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(tidyverse)

### Define the color palette ###
c1 <- c("seagreen4", "seagreen3", "seagreen2",
        "mediumpurple4", "mediumpurple3", "mediumpurple2",
        "orange4", "orange3", "orange2",
        "steelblue4", "steelblue3", "steelblue2",
        "tomato4", "tomato3", "tomato2")

c2 <- c("seagreen3", "mediumpurple3", "orange3","steelblue3", "tomato3")

### Insert data ###
library(readxl)
data <- read_excel("CO2_min.xlsx",
                  sheet = "Sheet1")
head(data)

### Reordenar a apresentacao dos grupos ###
data$Soil = factor(data$Soil, levels = unique(data$Soil))
data$Dilution = factor(data$Dilution, levels = unique(data$Dilution))
data$Code = factor(data$Code, levels = unique(data$Code))


# Reshape to Long format with melt function
df <- melt(data, na.rm = TRUE, value.name = "CCO2", id.vars = c("Sample_ID", "Soil", "Dilution", "Code"))


# Calculates mean, sd, se and IC
my_sum <- df %>%
  group_by(Soil, Dilution, Code, variable) %>%
  summarise( 
    n=n(),
    mean=mean(CCO2),
    sd=sd(CCO2)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))



# Plot line graph

g1 <- ggplot(my_sum, aes(x=variable, y=mean, group = Code)) + 
  geom_line(aes(color = Soil)) +
  geom_point(aes(color = Soil, shape = Dilution), size=3.5) +
  scale_shape_manual(values=c(15, 16, 17, 18), )+
  geom_errorbar(aes(x=variable, ymin=mean-se, ymax=mean+se), width=0.2)+
  theme_bw()+
  scale_color_manual("Soil", values = c2)+
  theme(legend.position="bottom")+
  ggtitle("Mineralization rate of ¹⁴C-Glyphosate")+
  labs(y = "¹⁴C-CO² Evolved (%)", x = "Incubation time (days)")+
  facet_wrap("Soil", ncol = 5, nrow = 1)

g1 + guides(col = FALSE)



#g1 <- ggplot(data, aes(x=Time, y=CCO2, group = Group)) + 
  #geom_line(aes(color = Soil)) +
  #geom_point(aes(color = Soil, shape = Dilution), size=2) +
  #scale_shape_manual(values=c(15, 16, 17, 18))+
  #theme_bw()+
  #scale_color_manual("Soil", values = c2)+
  #theme(legend.position="bottom")+
  #ggtitle("Mineralization rate of ¹⁴C-Glyphosate")+
  #labs(y = "¹⁴C-CO² Evolved (%)", x = "Incubation time (days)")+
  #facet_wrap("Soil", ncol = 5, nrow = 1)

#g1 + guides(col = FALSE)



#  save plot with 600 dpi resolution
dev.print(tiff, "CO2_mineralization3.tiff", compression = "lzw", res=600, height=10, width=25, units="cm")



##############################################################################################

# GA's plot
f1 <- ggplot(df, aes(x = Time, y = variable, linetype = variable, color = Soil)) + 
  geom_line() +
  geom_point(size = 3) + 
  theme(legend.position="none", axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y =element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), #element_rect(fill=NA, linetype=1, colour="black"),
        panel.background = element_blank(), panel.spacing = unit(1, "lines"),
        axis.text.x=element_text(size=16, face="bold", color="black"),
        axis.text.y=element_text(size=14, face="bold", color = "black"),
        axis.title.y=element_text(size=16, face = "bold"),
        axis.title.x = element_text (size = 16, face = "bold"), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        aspect.ratio = 1.5/1) + coord_flip() +
  scale_fill_manual(values=c("#8E6698", "#4477AA", "#AA7744", "#DDDD77", "#117733"))
f1
