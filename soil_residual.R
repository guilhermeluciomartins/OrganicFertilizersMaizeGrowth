# Organic residues residual effects

# Load packages
library(agricolae)
library(tidyverse)
library(ggpubr)
library(ggthemes)

# Get and set directory
getwd()
setwd("C:/Users/GuilhermeL/OneDrive/Mestrado/Artigo 2")

# Color palette
color2 <- c("#e2e2e2", "#c5c5c5", "#a8a8a8", "#8b8b8b", "#6e6e6e", "#515151")


# Open excel files
library(readxl)
soil <- read_excel("manuscript-data2.xlsx",
                    sheet = "Soil residual")
soil


# Drop columns
soil = subset(soil, select = -c(SampleID, pH, C.POM, C.MAOM, C.total, Cu, Zn))

# Separar dados por solo
soil_clay <- soil[1:24,]

soil_sand <- soil[25:48,]


### Calculate the three-way ANOVA and Tukey-HSD in factorial treatments ###
anova <- aov(Ca ~ Soil * Group * Type, data = soil)

# Extract ANOVA table
print(anova)
summary(anova)

# ANOVA for both soils
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
# Soil             1 0.8039  0.8039  49.261 5.90e-08 ***
# Group            4 1.8416  0.4604  28.212 4.34e-10 ***
# Type             1 0.0295  0.0295   1.811   0.1879    
# Soil:Group       4 0.8374  0.2093  12.828 2.44e-06 ***
# Soil:Type        1 0.0535  0.0535   3.276   0.0797 .  
# Group:Type       2 0.0577  0.0289   1.769   0.1868    
# Soil:Group:Type  2 0.0031  0.0016   0.096   0.9091    
# Residuals       32 0.5222  0.0163                     
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### Calculate the two-way ANOVA and Tukey-HSD in factorial treatments ###
anova <- aov(S ~ Group * Type, data = soil_clay)

# Perform Tukey's HSD test
Tukey_HSD <- TukeyHSD(anova, conf.level = 0.95)
# Print the results of the Tukey's HSD test
print(Tukey_HSD)

outHSD<- HSD.test(anova, c("Group", "Type"),
                  main = "Dry mass ~ Group + Type + Group:Type", console=TRUE)

# Post-hoc Tukey test for clay soil

#                  N.total groups
# CM:Fresh       1.2102100      a
# PL:Fresh       1.1346192      a
# PL:Composted   1.0646833      a
# FC:Fresh       1.0567121      a
# FC:Composted   0.9784472      a
# CM:Composted   0.9552815      a
# Fert:Control   0.8042424      a
# Unfert:Control 0.1011321      b

#                    P groups
# CM:Composted   2.290      a
# CM:Fresh       2.040      a
# PL:Composted   1.280     ab
# PL:Fresh       1.076     ab
# FC:Fresh       0.821     ab
# FC:Composted   0.235      b
# Fert:Control   0.097      b
# Unfert:Control 0.084      b

#                        K groups
# CM:Fresh       0.1671333      a
# PL:Fresh       0.1500111      a
# FC:Fresh       0.1479667      a
# Unfert:Control 0.1472000      a
# PL:Composted   0.1349333      a
# Fert:Control   0.1344222      a
# FC:Composted   0.1341667      a
# CM:Composted   0.1157667      a

#                       Ca groups
# CM:Fresh       1.3609695      a
# PL:Composted   1.1556509     ab
# FC:Fresh       0.8381402     bc
# PL:Fresh       0.6240221     cd
# CM:Composted   0.5924910     cd
# FC:Composted   0.4751661     cd
# Unfert:Control 0.3453754      d
# Fert:Control   0.2742471      d

#                       Mg groups
# PL:Composted   0.2287587      a
# PL:Fresh       0.2250927      a
# Unfert:Control 0.2096954      a
# FC:Fresh       0.2082290      a
# CM:Fresh       0.2016302      a
# CM:Composted   0.1986974      a
# FC:Composted   0.1605710      a
# Fert:Control   0.1099801      a


# Post-hoc Tukey test for sandy clay loam soil

#                  N.total groups
# CM:Fresh       0.8332375      a
# CM:Composted   0.7505481     ab
# FC:Composted   0.7104449    abc
# PL:Composted   0.6637918    abc
# PL:Fresh       0.6314540    abc
# FC:Fresh       0.6007641     bc
# Fert:Control   0.5240540      c
# Unfert:Control 0.5203970      c


#                       P groups
# CM:Fresh       6.9761728      a
# PL:Composted   3.7045679      b
# CM:Composted   3.6798765      b
# FC:Fresh       0.6644444      c
# PL:Fresh       0.6613580      c
# FC:Composted   0.5193827      c
# Fert:Control   0.1261728      c
# Unfert:Control 0.1237037      c


#                         K groups
# CM:Fresh       0.11088235      a
# PL:Fresh       0.09201961      a
# Unfert:Control 0.09023529      a
# PL:Composted   0.08207843      a
# FC:Composted   0.08054902      a
# FC:Fresh       0.07850980      a
# CM:Composted   0.07774510      a
# Fert:Control   0.07366667      a


#                       Ca groups
# PL:Composted   1.0089947      a
# CM:Composted   0.6775516     ab
# PL:Fresh       0.6474871      b
# CM:Fresh       0.5719592      b
# FC:Fresh       0.5565603      b
# FC:Composted   0.5176964     bc
# Unfert:Control 0.3263100     bc
# Fert:Control   0.1767207      c


#                        Mg groups
# CM:Fresh       0.15837138      a
# CM:Composted   0.12831015     ab
# PL:Composted   0.11584573     ab
# Unfert:Control 0.11291293     ab
# PL:Fresh       0.11144653     ab
# FC:Composted   0.09091690      b
# FC:Fresh       0.09018370      b
# Fert:Control   0.03446044      c


#                          S groups
# PL:Fresh       0.010797814      a
# FC:Composted   0.008038251      a
# PL:Composted   0.006890710      a
# CM:Fresh       0.006344262      a
# FC:Fresh       0.006098361      a
# CM:Composted   0.005633880      a
# Fert:Control   0.004240437      a
# Unfert:Control 0.003857923      a

#                          S groups
# FC:Fresh       0.012163934      a
# CM:Composted   0.011344262      a
# CM:Fresh       0.009950820      a
# Unfert:Control 0.009240437      a
# PL:Fresh       0.007273224      a
# PL:Composted   0.007136612      a
# FC:Composted   0.006508197      a
# Fert:Control   0.004076503      a




# Reshape data
data2 = reshape2::melt(soil_sand, id.vars = c("Group", "Type", "Soil", "Unified"), variable.name = "Nutrients")

# Calculates mean, sd, se and IC
my_sum <- data2 %>%
  group_by(Unified, Nutrients) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))  

# Reorder data
my_sum$Unified = factor(my_sum$Unified, levels = c("CT_unfert", "CT_fert", "FC_fresh", "FC_comp", "PL_fresh", "PL_comp", "CM_fresh", "CM_comp"))
my_sum$Nutrients = factor(my_sum$Nutrients, levels = c("S", "Ca", "Mg", "K", "P", "N.total"))

# Plot graphics
soil_final <-  ggplot(my_sum, aes(fill = Nutrients, x = Unified, y = mean, group = desc(Unified)))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = color2)+
  theme_hc()+
  #geom_errorbar(aes(x=Unified, ymin=mean-se, ymax=mean+se), width=0.4)+
  guides(fill = guide_legend(title = "Nutrients"))+
  theme(axis.title.y = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10))+
  xlab(NULL)+
  scale_y_reverse()+
  ylab("Soil nutrients residual effects (g kg-¹)")

soil_final




# Arrange plots on one page ###
figure5 <- ggarrange(nutrients, g1, soil_final, cplot, 
                    labels = c("A", "B", "C", "D"),
                    common.legend = F,
                    legend = "right",
                    ncol = 2, nrow = 2)

figure5

# save plot with 900 dpi resolution
dev.print(tiff, "figure5.tiff", compression = "lzw", res=900, height=25, width=38, units="cm")


# save as SGV file width 800 e height 600 
