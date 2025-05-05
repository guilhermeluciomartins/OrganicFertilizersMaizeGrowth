# Stacked barchart for carbon analysis - (Guilherme - 27/01/2023)

# 1. Prepare R environment -----------------------------------------------------

# Load packages
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(dplyr)
library(agricolae)

# get and set directory
getwd()
setwd("C:/Users/GuilhermeL/OneDrive/Mestrado/Artigo 2")
path <-"C:/Users/GuilhermeL/OneDrive/Mestrado/Artigo 2"
list.files(path)

# Color palette
c1 <- c("gray70", "gray30")

# Open excel files
library(readxl)
data <- read_excel("manuscript-data2.xlsx",
                  sheet = "Soil residual")
data

# Select data
data <- data %>% select(-SampleID, -pH, -C.total, -N.total, -P, -K, -Ca, -Mg, -S, -Cu, -Zn)







# 2. Statistical analysis ------------------------------------------------------

# Fit the three-way ANOVA model
model <- aov(C.total ~ Soil * Group * Type, data = data)

# Summarize the ANOVA results
summary(model)

# C.total
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
# Soil             1  384.8   384.8 401.935  < 2e-16 ***
# Group            4   37.0     9.2   9.654 3.08e-05 ***
# Type             1    0.6     0.6   0.624    0.435    
# Soil:Group       4    7.2     1.8   1.878    0.138    
# Soil:Type        1    0.1     0.1   0.099    0.755    
# Group:Type       2    0.7     0.3   0.342    0.713    
# Soil:Group:Type  2    0.4     0.2   0.212    0.810    
# Residuals       32   30.6     1.0                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#                 Df Sum Sq Mean Sq F value   Pr(>F)    
# Soil             1 15.227  15.227  58.864 9.59e-09 ***
# Group            4  9.432   2.358   9.116 4.94e-05 ***
# Type             1  0.411   0.411   1.589    0.217    
# Soil:Group       4  0.388   0.097   0.375    0.825    
# Soil:Type        1  0.697   0.697   2.693    0.111    
# Group:Type       2  0.665   0.332   1.285    0.291    
# Soil:Group:Type  2  0.565   0.283   1.092    0.348    
# Residuals       32  8.278   0.259                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Solos diferiram para C.POM, por isso os dados foram analisados separados. Houve diferença significativa entre groupos de resíduos organicos, mas nao houve entre os tipos de residuos organicos.


#                 Df Sum Sq Mean Sq  F value   Pr(>F)    
# Soil             1  553.1   553.1 1098.411  < 2e-16 ***
# Group            4   17.8     4.4    8.819 6.45e-05 ***
# Type             1    0.0     0.0    0.035   0.8536    
# Soil:Group       4    7.8     1.9    3.855   0.0115 *  
# Soil:Type        1    1.3     1.3    2.595   0.1170    
# Group:Type       2    1.7     0.8    1.677   0.2030    
# Soil:Group:Type  2    0.1     0.1    0.101   0.9045    
# Residuals       32   16.1     0.5 

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Houve diferenca significativa entre os solos para C.MAOM, alem disso houve diferenca entre os grupos de residuos organicos, mas nao entre a forma compostada e fresca.


# Two way ANOVA for each soil type to test the effects of organic residues types

# Filter data
data1 <- data[1:24,]  # CLAY SOIL
data2 <- data[25:48,] # SANDY CLAY LOAM SOIL

# Fit the model for clay soil
model <- aov(C.MAOM ~ Group * Type, data = data2)

# Summarize the ANOVA results
summary(model)


# ANOVA for clay soil 
# Df Sum Sq Mean Sq F value Pr(>F)  
# Group        4  3.627  0.9066   4.604 0.0115 *
# Type         1  0.019  0.0188   0.095 0.7616  
# Group:Type   2  1.226  0.6132   3.113 0.0721 .
# Residuals   16  3.151  0.1969                 

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# C.POM: Houve diferença significativa entre os grupos de residuos organicos, mas não para os tipos. Houve interação entre grupos e tipos de residuos organicos.
# C.MAOM: Não houve diferenças significativas.


# ANOVA for the sandy clay loam soil (C.MAOM)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Group        4 20.750   5.188  19.320 5.71e-06 ***
# Type         1  0.511   0.511   1.903    0.187    
# Group:Type   2  0.519   0.259   0.966    0.402    
# Residuals   16  4.296   0.269                     

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# houve diferenca apenas entre as fontes de residuos organicos (groups), mas nao entre as formas dos residuos (typeS)


# Visualize the results
par(mfrow=c(1,3))
boxplot(C.POM ~ Soil, data = data, main = "Soil Type")
boxplot(C.POM ~ Group, data = data, main = "Organic residues groups")
boxplot(C.POM ~ Type, data = data, main = "Residue type")

# 1. Homogeneity of variances
plot(model, 1)

# 2. Normality
plot(model, 2)

# 3. Spread-Location
plot(model, 3)

# Extract the residuals
aov_residuals <- residuals(object = model)

# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals)



# Fit the model for sandy clay loam soil
model <- aov(C.MAOM ~ Group * Type, data = data2)

# Summarize the ANOVA results
summary(model)

#         Df Sum Sq Mean Sq F value  Pr(>F)   
# Group        4  6.193  1.5483   4.832 0.00951 **
# Type         1  1.089  1.0889   3.399 0.08386 . 
# Group:Type   2  0.004  0.0018   0.006 0.99429   
# Residuals   16  5.126  0.3204                   
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# C.POM: Houve diferença significativa entre os Groups e Types de residuos organicos, mas não houve interação entre esses dois.
# C.MAOM: Houve diferenças significativas para os Groups de residuos organicos, mas nao para os tipos.



# Tukey-HSD post-hoc tests
# Perform an ANOVA
fit <- aov(C.MAOM ~ Group * Type, data = data2)
# Perform Tukey's HSD test
Tukey_HSD <- TukeyHSD(fit, conf.level = 0.95)
# Print the results of the Tukey's HSD test
print(Tukey_HSD)


outHSD<- HSD.test(fit,c("Group", "Type"),
                  main = "Dry mass ~ Soil + Group + Type + Group:Type", console=TRUE)


# Clay Soil statistics 

#                  C.MAOM groups
# CM:Fresh       14.60213      a
# PL:Fresh       14.19765      a
# FC:Composted   14.06228      a
# PL:Composted   13.93008      a
# FC:Fresh       13.92947      a
# CM:Composted   13.46181      a
# Fert:Control   13.07642      a
# Unfert:Control 12.93381      a

#                   C.POM groups
# FC:Fresh       2.222328      a
# PL:Fresh       2.188213      a
# FC:Composted   2.178091      a
# CM:Composted   2.039142      a
# PL:Composted   1.674790      a
# CM:Fresh       1.287810      a
# Unfert:Control 1.188692      a
# Fert:Control   1.110401      a



# Sandy clay loam soil statistics

#                 C.MAOM groups
# CM:Fresh       8.575588      a
# CM:Composted   8.434906      a
# FC:Composted   6.966390      b
# PL:Composted   6.836646      b
# Fert:Control   6.721183      b
# FC:Fresh       6.349193      b
# PL:Fresh       6.302157      b
# Unfert:Control 5.694522      b

#                   C.POM groups
# PL:Fresh       3.579042      a
# FC:Fresh       3.544121     ab
# FC:Composted   3.074538     ab
# PL:Composted   3.046817     ab
# CM:Fresh       2.977077     ab
# CM:Composted   2.503139     ab
# Fert:Control   2.210971     ab
# Unfert:Control 1.965326      b


# 3. Preparing data for plots --------------------------------------------------

# Reshape data (clay soil)
data_1 = reshape2::melt(data1, id.vars = c("Group", "Type", "Soil", "Unified"), variable.name = "Carbon")

# Calculates mean, sd, se and IC
my_sum_clay <- data_1 %>%
  group_by(Unified, Carbon) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate(se=sd/sqrt(n))  %>%
  mutate(ic=se * qt((1-0.05)/2 + .5, n-1))  

# Reorder data
my_sum_clay$Unified = factor(my_sum_clay$Unified, levels = c("CT_unfert", "CT_fert",
                                                   "FC_fresh", "FC_comp",
                                                   "PL_fresh", "PL_comp",
                                                   "CM_fresh", "CM_comp"))

my_sum_clay$Carbon = factor(my_sum_clay$Carbon, levels = c( "C.POM", "C.MAOM"))



# Reshape data
data_2 = reshape2::melt(data2, id.vars = c("Group", "Type", "Soil", "Unified"), variable.name = "Carbon")

# Calculates mean, sd, se and IC
my_sum_sandy <- data_2 %>%
  group_by(Unified, Carbon) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate(se=sd/sqrt(n))  %>%
  mutate(ic=se * qt((1-0.05)/2 + .5, n-1))  

# Reorder data
my_sum_sandy$Unified = factor(my_sum$Unified, levels = c("CT_unfert", "CT_fert",
                                                   "FC_fresh", "FC_comp",
                                                   "PL_fresh", "PL_comp",
                                                   "CM_fresh", "CM_comp"))

my_sum_sandy$Carbon = factor(my_sum$Carbon, levels = c( "C.POM", "C.MAOM"))




# 4. Plot ----------------------------------------------------------------------

# Plot graphics
g1 <-  ggplot(my_sum_clay, aes(x = Unified, y = mean))+
  geom_bar(position = "stack", stat = "identity", color="black", fill = c1)+
  geom_errorbar(aes(x=Unified, ymin=mean-ic, ymax=mean+ic), width=0.4)+
  #geom_text(data = value_max, aes(x=Unified, y = 0.1 + max_value, label = sig.letters$groups), vjust=-1)+
  theme_hc()+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(),
        legend.position = "right",
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10))+
  xlab(NULL)+
  ylab("Dry mass weight (g plant-¹)")

g1


# Plot graphics (clay soil)
cplot1 <- ggplot(my_sum_clay, aes(fill = Carbon, x = Unified, y = mean, group = desc(Unified)))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = c1)+
  geom_errorbar(aes(x=Unified, ymin=mean-se, ymax=mean+se), width=0.4)+
  theme_hc()+
  theme(axis.title.y = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10))+
  scale_y_reverse(limits = c(17, 0))+
  xlab(NULL)+
  ylab("Fractions of soil organic carbon (g kg-¹)")
  
cplot1



# Plot graphics (Sandy clay loam soil)
cplot <- ggplot(my_sum_sandy, aes(fill = Carbon, x = Unified, y = mean, group = desc(Unified)))+
          geom_bar(position = "stack", stat = "identity")+
          scale_fill_manual(values = c1)+
          geom_errorbar(aes(x=Unified, ymin=mean-se, ymax=mean+se), width=0.4)+
          theme_hc()+
          theme(axis.title.y = element_text(size = 12),
          legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 10))+
          scale_y_reverse(limits = c(17, 0))+
          xlab(NULL)+
          ylab("Fractions of soil organic carbon (g kg-¹)")

cplot


# Merge the plots
figure_top <- ggarrange(g1_clay, g1_sand,  
                    #labels = c("A", "B"),
                    common.legend = T,
                    legend = "right",
                    ncol = 2, nrow = 1)

figure_top

# Saved as SVG (h = 1100 X w = 450)

figure_bottom <- ggarrange(cplot1, cplot,  
                        #labels = c("A", "B"),
                        common.legend = T,
                        legend = "right",
                        ncol = 2, nrow = 1)

figure_bottom

# Saved as SVG (h = 1100 X w = 450)


ggarrange(figure_top, figure_bottom,  
          #labels = c("A", "B"),
          common.legend = F,
          legend = "right",
          ncol = 1, nrow = 2)

