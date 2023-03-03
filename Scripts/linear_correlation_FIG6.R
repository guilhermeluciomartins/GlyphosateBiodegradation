### Script para correlação ### (Guilherme Martins - 03/10/2022)


### Installing packages:###
install.packages("Hmisc")
install.packages("corrplot")
install.packages("PerformanceAnalytics")
install.packages("psych")

### Loading packages###
library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(patchwork)
library(psych)
library(corrplot)
library(ggpubr)


### Get and Set working directory
getwd()
setwd("C:/Users/gui_l/OneDrive/Doutorado/Colaborações/Experimento glifosato")


# Definir paleta de cores
c1 <- c("seagreen4", "seagreen3", "seagreen2", "seagreen1",
        "mediumpurple4", "mediumpurple3", "mediumpurple2", "mediumpurple1",
        "orange4", "orange3", "orange2", "orange1",
        "steelblue4", "steelblue3", "steelblue2", "steelblue1",
        "tomato4", "tomato3", "tomato2", "tomato1")

c2 <- c("seagreen3", "mediumpurple3", "orange3","steelblue3", "tomato3")


### Load data for 16S
db <- readxl::read_xlsx("corr_div_CO22.xlsx", sheet = "16S", col_names = TRUE, col_types = NULL, na = "NA")


### Reorganize data
db$Soil <- factor(db$Soil, levels=unique(db$Soil))
db$Dilution <- factor(db$Dilution, levels=unique(db$Dilution))
db$Code <- factor(db$Code, levels=unique(db$Code))


### Plot correlation
g1 <- ggscatter(db, y = "C_CO2", x = "Observed_ASV",
               add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Observed ASV", ylab = "14C-CO² Evolved (%)", title = "Bacterial diversity and glyphosate mineralization",
               add.params = list(color = "red", fill = "lightgray"),
               color= "Soil", palette = c2,
               shape = "Dilution", size=3)+
  theme_bw()+
  scale_shape_manual(values=c(15, 16, 17, 18))

g1

###########################################################################

### Load data for ITS
db <- readxl::read_xlsx("corr_div_CO22.xlsx", sheet = "ITS", col_names = TRUE, col_types = NULL, na = "NA")


### Reorganize data
db$Soil <- factor(db$Soil, levels=unique(db$Soil))
db$Dilution <- factor(db$Dilution, levels=unique(db$Dilution))
db$Code <- factor(db$Code, levels=unique(db$Code))


### Plot correlation
g2 <- ggscatter(db, y = "C_CO2", x = "Observed_ASV",
                add = "reg.line", conf.int = TRUE,
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "Observed ASV", ylab = "14C-CO² Evolved (%)", title = "Fungal diversity and glyphosate mineralization",
                add.params = list(color = "red", fill = "lightgray"),
                color= "Soil", palette = c2,
                shape = "Dilution", size=3)+
  theme_bw()+
  scale_shape_manual(values=c(15, 16, 17, 18))

g2




### Arrange plots on one page ###
ggarrange(g1, g2, 
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 2, nrow = 1)


#  save plot with 800 dpi resolution
dev.print(tiff, "correlation_CO2_div2.tiff", compression = "lzw", res=800, height=15, width=25, units="cm")



################
# Ga plot

xy<-ggplot(db, aes(NaOHi, pqqc2)) + 
  geom_point(size = 4) +
  facet_wrap(db$Genotype, nrow = 1)+
  geom_smooth(method = lm, se =T, color = "black") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "`,`")),
           r.accuracy = 0.01,
           p.accuracy = 0.01, size = 4,
           label.x = 10, label.y = 150) +
  stat_regline_equation(aes(label = ..eq.label..),
                        label.x = 10, label.y = 200, size = 4) +
  theme_bw()+
  theme(legend.position="right",plot.title = element_text(size=11), axis.text.y=element_text(size = 13), axis.text.x=element_text(size = 13)) +
  ggtitle("") +
  xlab("pqqC gene copies") +
  ylab("Inorganic P") +
  scale_fill_manual(values=c("#8E6698", "#4477AA", "#AA7744", "#DDDD77", "#117733"))
xy
