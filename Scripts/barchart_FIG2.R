### Barchart for enzymatic activity ### (Guilherme MartiCT - 12/12/2022)

rm(list=ls(all=TRUE))
getwd()

setwd("C:/Users/glmartins/OneDrive/Doutorado/Colaborações/Experimento glifosato/Dados")
path <-"C:/Users/glmartins/OneDrive/Doutorado/Colaborações/Experimento glifosato/Dados"
list.files(path)

# Load ggplot2
library(ggplot2)
library(dplyr)
library(ggpubr)
library(agricolae)
library(multcompView)

### Define the color palette ###
c1 <- c("seagreen4", "seagreen3", "seagreen2",
        "mediumpurple4", "mediumpurple3", "mediumpurple2",
        "orange4", "orange3", "orange2",
        "steelblue4", "steelblue3", "steelblue2",
        "tomato4", "tomato3", "tomato2")

c2 <- c("seagreen3", "mediumpurple3", "orange3","steelblue3", "tomato3")


### Open excel files ###
library(readxl)
data <- read_excel("enzymes.xlsx")


### Reorganize data ###
data$Soil = factor(data$Soil, levels = unique(data$Soil))
data$Dilution = factor(data$Dilution, levels = unique(data$Dilution))


############ BETA-GLUCOSIDASE ############ 


### Calculate ANOVA and Tukey-HSD in factorial treatments ###
# Comparison of treatment combination
anova <- aov(B.gli ~ Dilution * Soil, data = data)

# Extract ANOVA table
print(anova)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

# table with factors and 3rd quantile
dt <- data %>% 
  group_by(Soil, Dilution) %>%
  summarise(w=mean(B.gli), sd = sd(B.gli) / sqrt(n())) %>%
  arrange(desc(w)) %>% 
  ungroup() %>% 
  mutate(Dilution = factor(Dilution,
                            levels = c("CS",
                                       "D1",
                                       "D3",
                                       "SS"),
                            ordered = TRUE))  %>%
  mutate(Soil = factor(Soil,
                       levels = c("NF",
                                  "CT",
                                  "CTN",
                                  "NT",
                                  "NTN"),
                       ordered = TRUE))

# extracting the compact letter display and adding to the Tk table
cld2 <- data.frame(letters = cld$`Dilution:Soil`$Letters)
dt$cld <- cld2$letters

print(dt)
# A tibble: 20 × 5
#Soil  Dilution     w    sd cld  
#<chr> <ord>    <dbl> <dbl> <chr>
#1 NF    CT       189.  5.20  a    
#2 NT    CT       161.  6.92  b    
#3 NTN   CT       159.  7.71  b    
#4 CFN   D3        69.8 1.80  c    
#5 CFN   SS        67.7 4.41  c    
#6 CFN   CT        64.7 2.21  c    
#7 CF    CT        62.8 1.76  c    
#8 CFN   D1        55.7 2.21  cd   
#9 CF    D3        38.1 4.78  de   
#10 NTN   SS        36.2 5.87  e    
#11 CF    D1        36.0 0.408 ef   
#12 CF    SS        33.3 3.37  efg  
#13 NTN   D3        33.1 0.992 efg  
#14 NT    D3        28.5 1.25  efg  
#15 NF    SS        26.3 2.93  efg  
#16 NF    D1        21.5 1.19  efg  
#17 NT    SS        21.0 1.31  efg  
#18 NT    D1        17.0 0.713 fg   
#19 NTN   D1        15.9 0.413 g    
#20 NF    D3        14.6 0.928 g   


# Plot Bar Chart with Standard Deviation
g1 <- ggplot(dt, aes(x = Dilution, y = w)) +
  geom_bar(stat="identity", aes(fill=Soil), show.legend = FALSE) +
  geom_errorbar(aes(x=Dilution, ymin=w-sd, ymax=w+sd), width=0.2) +
  ggtitle("β-Glucosidase activity")+
  geom_text(aes(label = cld, x = Dilution, y = w+sd), vjust = -0.5) +
  theme_bw()+
  scale_fill_manual(values = c2)+
  xlab(NULL)+
  ylab("µg nitrophenol g-¹ soil h-¹")+
  facet_wrap("Soil", ncol = 5, nrow = 1)

g1


############ ACID PHOSPHATASE ############ 

### Calculate ANOVA and Tukey-HSD in factorial treatments ###
# Comparison of treatment combinations
anova <- aov(Acid.phos ~ Dilution * Soil, data = data)

# Extract ANOVA table
print(anova)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

# table with factors and 3rd quantile
dt <- data %>% 
  group_by(Soil, Dilution) %>%
  summarise(w=mean(Acid.phos), sd = sd(Acid.phos) / sqrt(n())) %>%
  arrange(desc(w)) %>% 
  ungroup() %>% 
  mutate(Dilution = factor(Dilution,
                           levels = c("CS",
                                      "D1",
                                      "D3",
                                      "SS"),
                           ordered = TRUE))  %>%
  mutate(Soil = factor(Soil,
                       levels = c("NF",
                                  "CT",
                                  "CTN",
                                  "NT",
                                  "NTN"),
                       ordered = TRUE))

# extracting the compact letter display and adding to the Tk table
cld2 <- data.frame(letters = cld$`Dilution:Soil`$Letters)
dt$cld <- cld2$letters

print(dt)
# A tibble: 20 × 5
#Soil  Dilution     w    sd cld  
#<ord> <ord>    <dbl> <dbl> <chr>
#1 NT    CT       798.  65.7  a    
#2 NTN   CT       751.  64.4  a    
#3 NF    CT       702.  17.1  a    
#4 CFN   CT       337.  20.5  b    
#5 CF    CT       294.  21.3  b    
#6 NF    D1       253.  19.4  bc   
#7 NF    D3       158.  19.2  cd   
#8 NTN   D3       139.  12.4  cd   
#9 NT    D1       139.  11.4  cd   
#10 NTN   D1       130.   9.35 cd   
#11 NT    D3       126.   8.52 cd   
#12 NF    SS       125.  11.5  cd   
#13 CFN   D1       125.  12.3  cd   
#14 CFN   D3       116.   6.18 d    
#15 NTN   SS       113.   4.48 d    
#16 CFN   SS       110.   4.65 d    
#17 NT    SS        98.0 22.2  d    
#18 CF    D1        86.7 10.0  d    
#19 CF    D3        83.0  4.47 d    
#20 CF    SS        43.8  2.73 d 

# Plot Bar Chart with Standard Deviation
g2 <- ggplot(dt, aes(x = Dilution, y = w)) +
  geom_bar(stat="identity", aes(fill=Soil), show.legend = FALSE) +
  geom_errorbar(aes(x=Dilution, ymin=w-sd, ymax=w+sd), width=0.2) +
  ggtitle("Acid phosphatase activity")+
  geom_text(aes(label = cld, x = Dilution, y = w+sd), vjust = -0.5) +
  theme_bw()+
  scale_fill_manual(values = c2)+
  xlab(NULL)+
  ylab("µg nitrophenol g-¹ soil h-¹")+
  facet_wrap("Soil", ncol = 5, nrow = 1)

g2


### Arrange plots on one page ###
ggarrange(g1, g2, 
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "none",
          ncol = 1, nrow = 2)


#  save plot with 600 dpi resolution
dev.print(tiff, "enzyme_activity2.tiff", compression = "lzw", res=600, height=23, width=25, units="cm")
