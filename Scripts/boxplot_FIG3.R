### Boxplot for diversity index ###
rm(list=ls(all=TRUE)) #apaga todo o ambiente R

### Definir diretorio de trabalho ###
getwd()
setwd("C:/Users/glmartins/OneDrive/Doutorado/Colaborações/Experimento glifosato")
path<-"C:/Users/glmartins/OneDrive/Doutorado/Colaborações/Experimento glifosato"
list.files(path)

### Carregar pacotes ###
library(multcompView)
library(tidyverse)
library(tidyr)
library(ggthemes)
library(patchwork)
library(ggplot2)
library(ggpubr)

### Define color palette ###
c1 <- c("seagreen4", "seagreen3", "seagreen2",
        "mediumpurple4", "mediumpurple3", "mediumpurple2",
        "orange4", "orange3", "orange2",
        "steelblue4", "steelblue3", "steelblue2",
        "tomato4", "tomato3", "tomato2")

c2 <- c("seagreen3", "mediumpurple3", "orange3","steelblue3", "tomato3")


### Abrir conjunto de dados ###
library(readxl)
div <- read_excel("diversity.xlsx",
                      sheet = "ITS")
div

### Selecionar dados do dataframe com dplyr ###
df_div <- div %>% select(Soil, Dilution, Code, Shannon, Pielou.evenness)
df_div

### Dar nome aos dados a serem analisados ###
Soil=df_div$Soil
Dilution=df_div$Dilution
Code=df_div$Code
Diversity=df_div$Shannon
Evenness=df_div$Pielou.evenness

### Criar dataframe com os dados ###
data=data.frame(Soil, Dilution, Code, Diversity, Evenness)

### Reordenar a apresentacao dos grupos ###
data$Soil = factor(data$Soil, levels = unique(data$Soil))
data$Dilution = factor(data$Dilution, levels = unique(data$Dilution))


# Reshapre data
data_1 = reshape2::melt(data, id.vars = c("Soil", "Dilution", "Code"), variable.name = "diversity")

# Compute the analysis of variance
res.aov <- aov(Evenness ~ Soil * Dilution, data = data)
# Summary of the analysis
summary(res.aov)

# ANOVA for Shannon index
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Soil           4  0.922   0.230   3.641  0.01560 *  
# Dilution       2 17.826   8.913 140.802 5.66e-16 ***
# Soil:Dilution  8  2.697   0.337   5.325  0.00033 ***
# Residuals     30  1.899   0.063                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ANOVA for Evenness index
# Df  Sum Sq Mean Sq F value   Pr(>F)    
# Soil           4 0.03541 0.00885   1.826  0.14986    
# Dilution       2 0.17256 0.08628  17.799 8.01e-06 ***
# Soil:Dilution  8 0.12489 0.01561   3.220  0.00921 ** 
# Residuals     30 0.14543 0.00485                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# we can conclude that there are significant differences between the groups highlighted with “*" in the model summary

print(model.tables(res.aov,"means"),digits=3) 

TukeyHSD(res.aov)

# 1. Homogeneity of variances
plot(res.aov, 1)

# 2. Normality
plot(res.aov, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

# ANOVA test with no assumption of equal variances
oneway.test(Evenness ~ Soil * Dilution, data = data)


### Gerar bloxplots ###
# Fungi
g1 <- ggboxplot(data_1, "Dilution", "value", fill = "Soil")+
  geom_jitter(height=0.1, width=0.15, alpha=0.15, size=1.5)+
  stat_compare_means(
    comparisons = list( c("CS", "D1"), c("D1", "D3"), c("CS", "D3")), 
    label = "p.signif", method = "t.test"
  )+
  theme_bw()+
  labs(x = NULL, y = "Diversity index")+
  scale_fill_manual(values = c2)+
  theme(legend.position="bottom")


fung <- g1 + facet_wrap("diversity", scales = "free_y") + ggtitle("Soil fungal diversity")
fung



########################################################################################################

### 16S Bacterial diversity ###

### Abrir conjunto de dados ###
library(readxl)
div2 <- read_excel("diversity.xlsx",
                  sheet = "16S")
div2

### Selecionar dados do dataframe com dplyr ###
df_div2 <- div2 %>% select(Soil, Dilution, Code, Shannon, Pielou.evenness)
df_div2

### Dar nome aos dados a serem analisados ###
Soil=df_div2$Soil
Dilution=df_div2$Dilution
Code=df_div2$Code
Diversity=df_div2$Shannon
Evenness=df_div2$Pielou.evenness

### Criar dataframe com os dados ###
data2=data.frame(Soil, Dilution, Code, Diversity, Evenness)

### Reordenar a apresentacao dos grupos ###
data2$Soil = factor(data2$Soil, levels = unique(data2$Soil))
data2$Dilution = factor(data2$Dilution, levels = unique(data2$Dilution))

# Compute the analysis of variance
res.aov <- aov(Evenness ~ Soil * Dilution, data = data2)
# Summary of the analysis
summary(res.aov)

# ANOVA for Shannon index
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Soil           4  0.922   0.230   3.641  0.01560 *  
# Dilution       2 17.826   8.913 140.802 5.66e-16 ***
# Soil:Dilution  8  2.697   0.337   5.325  0.00033 ***
# Residuals     30  1.899   0.063                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# ANOVA for Evenness index
# Df  Sum Sq Mean Sq F value   Pr(>F)    
# Soil           4 0.09188 0.02297  33.445 1.16e-10 ***
# Dilution       2 0.15970 0.07985 116.268 7.40e-15 ***
# Soil:Dilution  8 0.04435 0.00554   8.073 9.39e-06 ***
# Residuals     30 0.02060 0.00069                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# we can conclude that there are significant differences between the groups highlighted with “*" in the model summary

print(model.tables(res.aov,"means"),digits=3) 

TukeyHSD(res.aov)

# 1. Homogeneity of variances
plot(res.aov, 1)

# 2. Normality
plot(res.aov, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

# ANOVA test with no assumption of equal variances
oneway.test(Shannon ~ Soil * Dilution, data = data)


# Reshapre data
data_2 = reshape2::melt(data2, id.vars = c("Soil", "Dilution", "Code"), variable.name = "diversity")


### Gerar bloxplot ###
# Bacteria
g2 <- ggboxplot(data_2, "Dilution", "value", fill = "Soil")+
  geom_jitter(height=0.1, width=0.15, alpha=0.15, size=1.5)+
  stat_compare_means(
    comparisons = list( c("CS", "D1"), c("D1", "D3"), c("CS", "D3")), 
    label = "p.signif", method = "t.test"
  )+
  theme_bw()+
  labs(x = NULL, y = "Diversity index")+
  scale_fill_manual(values = c2)+
  theme(legend.position="bottom")
  

bac <- g2 + facet_wrap("diversity", scales = "free_y") + ggtitle("Soil bacterial diversity")
bac


### Arrange plots on one page ###
ggarrange(bac, fung, 
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 1, nrow = 2)

#  save plot with 600 dpi resolution
dev.print(tiff, "alpha_diversity2.tiff", compression = "lzw", res=600, height=12, width=20, units="cm")
