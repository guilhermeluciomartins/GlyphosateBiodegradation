# NMDS - Guilherme Martins

### Install pairwiseAdonis package ###
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Carregar bibliotecas
library (vegan)
library(ape)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(cluster)
library(dplyr)
library(tidyr)
library(ARTool)
library(lsmeans)
library(pairwiseAdonis)
library(ggpubr)

# Definir diretorio dos arquivos
getwd()
setwd("C:/Users/glmartins/OneDrive/Doutorado/Colaborações/Experimento glifosato")

# Carregar dados
bac <- read.table("NMDS_16S.txt", header = TRUE, sep = "\t", dec = ".")
bac$Soil <- factor(bac$Soil, levels=unique(bac$Soil))
bac$Dilution <- factor(bac$Dilution, levels=unique(bac$Dilution))
bac$Code <- factor(bac$Code, levels=unique(bac$Code))


# Definir paleta de cores
c1 <- c("seagreen4", "seagreen3", "seagreen2",
        "mediumpurple4", "mediumpurple3", "mediumpurple2",
        "orange4", "orange3", "orange2",
        "steelblue4", "steelblue3", "steelblue2",
        "tomato4", "tomato3", "tomato2")

c2 <- c("seagreen3", "mediumpurple3", "orange3","steelblue3", "tomato3")

############################################################################
# NMDS para dados de 16S


# Permanova
adonis2(bac[5:2145] ~ Soil*Dilution, data = bac, permutations = 999, method = "bray", p.adjusted = p.adjust(p.value,method="fdr"))
#Df SumOfSqs      R2       F Pr(>F)    
#Soil           4   4.3850 0.29704 14.2944  0.001 ***
#Dilution       2   4.2355 0.28691 27.6141  0.001 ***
#Soil:Dilution  8   3.8414 0.26021  6.2611  0.001 ***
#Residual      30   2.3007 0.15585                   
#Total         44  14.7627 1.00000    

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interação muito siginificativa entre Solos e niveis de Diluição


pairwise.adonis(bac[5:2145], bac$Soil, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#pairs Df SumsOfSqs   F.Model         R2 p.value  p.adjusted sig
#1    NF vs CF  1 1.6147841 6.2227222 0.28001620   0.001 0.002500000   *
#2   NF vs CFN  1 1.7201743 6.4582190 0.28756595   0.001 0.002500000   *
#3    NF vs NT  1 1.6401075 6.7363477 0.29628099   0.001 0.002500000   *
#4   NF vs NTN  1 1.6712298 6.7641750 0.29714123   0.001 0.002500000   *
#5   CF vs CFN  1 0.2214825 0.7891567 0.04700395   0.517 0.517000000    
#6    CF vs NT  1 0.8469736 3.2857183 0.17037054   0.010 0.014285714   .
#7   CF vs NTN  1 0.8786407 3.3616267 0.17362315   0.021 0.026250000   .
#8   CFN vs NT  1 1.0231242 3.8662368 0.19461345   0.004 0.006666667   *
#9  CFN vs NTN  1 1.0922140 4.0719320 0.20286697   0.002 0.004000000   *
#10  NT vs NTN  1 0.2538446 1.0346354 0.06073716   0.257 0.285555556    

# Diferenças significativa entre NF vs CF; NF vs CFN; NF vs NT; NF vs NTN; CFN vs NTN
# Não houve diferenças significativas entre uso de N em nenhum solo e pouca diferença entre CF vs NT


pairwise.adonis(bac[5:2145], bac$Dilution, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1 CT vs D1  1 2.5765013 11.061976 0.28319039   0.001     0.0015   *
#2 CT vs D3  1 2.9922073 11.833548 0.29707492   0.001     0.0015   *
#3 D1 vs D3  1 0.7845552  2.947615 0.09524531   0.007     0.0070   *

# Houve diferenças significativas entre os niveis de diluição


# Ordenacao NMDS
  matrix.nmds<-metaMDS(bac[5:2145], k=2, distance = 'bray', noshare=FALSE, autotransform=TRUE) # k=2 transforma em nmds de duas dimensoes, autotransform = true para dados de sequencias
matrix.nmds$stress #0.1166385
# A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation. 

NMDS1<-matrix.nmds$points[,1]
NMDS2<-matrix.nmds$points[,2]
NMDS=data.frame(NMDS1=NMDS1, NMDS2=NMDS2, Soil=bac$Soil, Dilution=bac$Dilution, Code=bac$Code)

# Plotar a NMDS
NMDS_16S <- ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=Soil, shape=Dilution), size=4,  alpha=0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  stat_ellipse(aes(color = Soil), type = "norm", linetype = 2, alpha = 0.5)+
  scale_colour_manual(values=c2) +
  scale_shape_manual(values=c(15, 16, 17, 18)) +
  theme_bw() +
  ggtitle("Bacterial community structure")+
  annotate("text", x = 5.3, y = 3.5, hjust = 1 , label = "
  stress = 0.116
  
  PERMANOVA: p < 0.001
  Soil: R² = 0.297
  Dilution: R² = 0.286", size = 2.5)+
  theme(legend.position="right",
        legend.text = element_text(size=10),
        axis.text.x= element_text(angle= 0, size= 10),
        axis.text.y= element_text(size= 10),
        axis.title= element_text(size= 10, face= "bold"))

NMDS_16S


############################################################################
# NMDS para dados de ITS


# Carregar dados
fun <- read.table("NMDS_ITS.txt", header = TRUE, sep = "\t", dec = ".")
fun$Soil <- factor(fun$Soil, levels=unique(fun$Soil))
fun$Dilution <- factor(fun$Dilution, levels=unique(fun$Dilution))
fun$Code <- factor(fun$Code, levels=unique(fun$Code))


# Permanova
adonis2(fun[5:700] ~ Soil*Dilution, data = fun, permutations = 999, method = "bray", p.adjusted = p.adjust(p.value,method="fdr"))
#Df SumOfSqs      R2      F Pr(>F)    
#Soil           4   4.4666 0.30748 7.7976  0.001 ***
#Dilution       2   1.7792 0.12248 6.2121  0.001 ***
#Soil:Dilution  8   3.9844 0.27429 3.4779  0.001 ***
#Residual      30   4.2961 0.29575                  
#Total         44  14.5262 1.00000                  

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Interação muito siginificativa entre Solos e niveis de Diluição


pairwise.adonis(fun[5:700], fun$Soil, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1    NF vs CF  1 1.7764885 5.9299310 0.27040354   0.001 0.00250000   *
#2   NF vs CFN  1 2.0855277 7.6645571 0.32388340   0.001 0.00250000   *
#3    NF vs NT  1 2.0874397 7.8708586 0.32972667   0.001 0.00250000   *
#4   NF vs NTN  1 2.0650938 8.1550540 0.33761274   0.001 0.00250000   *
#5   CF vs CFN  1 0.2387090 0.9039101 0.05347343   0.461 0.51222222    
#6    CF vs NT  1 0.6360302 2.4729413 0.13386831   0.040 0.05000000   .
#7   CF vs NTN  1 0.7835476 3.1953707 0.16646569   0.016 0.02666667   .
#8   CFN vs NT  1 0.6465987 2.8147731 0.14960442   0.022 0.03142857   .
#9  CFN vs NTN  1 0.7450875 3.4220138 0.17619253   0.007 0.01400000   .
#10  NT vs NTN  1 0.1018728 0.4831651 0.02931264   0.726 0.72600000   

# Diferenças significativa entre NF vs CF; NF vs CFN; NF vs NT; NF vs NTN; CFN vs NTN
# Não houve diferenças significativas entre uso de N em nenhum solo e pouca diferença entre CF vs NT


pairwise.adonis(fun[5:700], fun$Dilution, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#1 CT vs D1  1 1.1708426 4.206259 0.13060378   0.003      0.006   *
#2 CT vs D3  1 0.9963356 3.174164 0.10182034   0.004      0.006   *
#3 D1 vs D3  1 0.5015672 1.575988 0.05328607   0.118      0.118   

# Houve diferenças significativas entre os niveis de diluição


# Ordenacao NMDS
matrix.nmds<-metaMDS(fun[5:700], k=2, distance = 'bray', noshare=FALSE, autotransform=TRUE) # k=2 transforma em nmds de duas dimensoes, autotransform = true para dados de sequencias
matrix.nmds$stress #0.09273138
# A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation. 

NMDS1<-matrix.nmds$points[,1]
NMDS2<-matrix.nmds$points[,2]
NMDS=data.frame(NMDS1=NMDS1, NMDS2=NMDS2, Soil=fun$Soil, Dilution=fun$Dilution, Code=fun$Code)

# Plotar a NMDS
NMDS_ITS <- ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=Soil, shape=Dilution), size=4,  alpha=0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  stat_ellipse(aes(color = Soil), type = "norm", linetype = 2, alpha = 0.5)+
  scale_colour_manual(values=c2) +
  scale_shape_manual(values=c(15, 16, 17, 18)) +
  theme_bw() +
  ggtitle("Fungal community structure")+
  annotate("text", x = 2.7, y = 2.4, hjust = 1 , label = "
  stress = 0.092
  
  PERMANOVA: p < 0.001
  Soil: R² = 0.307
  Dilution: R² = 0.122", size = 2.5)+
  theme(legend.position="right",
        legend.text = element_text(size=10),
        axis.text.x= element_text(angle= 0, size= 10),
        axis.text.y= element_text(size= 10),
        axis.title= element_text(size= 10, face= "bold"))

NMDS_ITS


### Arrange plots on one page ###
ggarrange(NMDS_16S, NMDS_ITS, 
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 2, nrow = 1)


#  save plot with 600 dpi resolution
dev.print(tiff, "NMDS_16S_ITS.tiff", compression = "lzw", res=600, height=12, width=20, units="cm")
