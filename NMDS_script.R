# NMDS - Guilherme Martins

### Install pairwiseAdonis package ###
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Carregar bibliotecas
library(vegan)
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

# Definir diretorio dos arquivos
getwd()
setwd("C:/Users/gui_l/OneDrive/Mestrado/Artigo 2")

# Carregar dados
data <- read.table("MDS_table2.txt", header = TRUE, sep = "\t", dec = ".")
data$Group <- factor(data$Group, levels=unique(data$Group))
data$Type <- factor(data$Type, levels=unique(data$Type))
data$Soil <- factor(data$Soil, levels=unique(data$Soil))

# Definir paleta de cores
c1 <- c("gray30", "tomato2", "tomato4", "seagreen2", "seagreen4", "steelblue2", "steelblue4", "gray70")
c2 <- c("gray30", "tomato3", "seagreen3", "steelblue3", "gray70")
c3 <- c("gray70","gray40", "#59A14F", "#EDC948", "#9C755F")


# Select rows for the clay and sandy soils
data_clay <- data[c(1:24), ]
data_sandy <- data[c(25:48), ]



############################################################################
# NMDS para dados de 16S


# Permanova for the clay soil

adonis2(data_clay[5:8132] ~ Group*Type, data = data_clay, permutations = 999, method = "bray", p.adjusted = p.adjust(p.value,method="fdr"))
#           Df SumOfSqs      R2      F Pr(>F)    
# Group       3   1.0183 0.19247 2.0709  0.002 ** 
# Type        2   0.8896 0.16814 2.7138  0.001 ***
# Group:Type  2   0.7605 0.14373 2.3198  0.003 ** 
# Residual   16   2.6225 0.49566                  
# Total      23   5.2909 1.00000                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Houve diferenças significativas entre Groups, Types e Soils. Alem disso, houve interação entre Groups e Types de residuos organicos e entre Soil:Groups:Types.


pairwise.adonis(data_clay[5:8132], data_clay$Group, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#          pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 Control vs FC  1 0.2022611 1.282480 0.1136700   0.167     0.2004    
# 2 Control vs PL  1 0.3849957 2.052862 0.1703215   0.114     0.1815    
# 3 Control vs CM  1 0.4306321 1.932532 0.1619549   0.004     0.0240   .
# 4      FC vs PL  1 0.3301610 1.615091 0.1390511   0.121     0.1815    
# 5      FC vs CM  1 0.3731647 1.556703 0.1347013   0.045     0.1350    
# 6      PL vs CM  1 0.3154234 1.170208 0.1047615   0.241     0.2410 

# Houve uma pequena diferença significativas entre o Controle e o CM


pairwise.adonis(data_clay[5:8132], data_clay$Type, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#                pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1      Unfert vs Fert  1 0.3050313 3.057075 0.4331929   0.100      0.119    
# 2     Unfert vs Fresh  1 0.4828549 1.889353 0.1589114   0.040      0.080    
# 3 Unfert vs Composted  1 0.1857783 1.256111 0.1115937   0.119      0.119    
# 4       Fert vs Fresh  1 0.4430098 1.715017 0.1463948   0.072      0.108    
# 5   Fert vs Composted  1 0.3164971 2.100939 0.1736178   0.014      0.042   .
# 6  Fresh vs Composted  1 0.5845745 2.553422 0.1376254   0.006      0.036   .

# Houve diferenças entre Fresh vs Compost! Tambem houve diferença entre Fert e Composted wastes

# ------------------------------------------------------------------------------

# Permanova for the sandy clay loam soil
adonis2(data_sandy[5:8132] ~ Group*Type, data = data_sandy, permutations = 999, method = "bray", p.adjusted = p.adjust(p.value,method="fdr"))
#            Df SumOfSqs      R2      F Pr(>F)    
# Group       3   0.7737 0.18117 1.8100  0.010 ** 
# Type        2   0.7294 0.17080 2.5596  0.001 ***
# Group:Type  2   0.4875 0.11417 1.7108  0.026 *  
# Residual   16   2.2797 0.53386                  
# Total      23   4.2703 1.00000                  

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Houve diferenca significativa entre os Grupos de residuos organicos, entre os tipos e houve interacao entre os dois fatores.


pairwise.adonis(data_sandy[5:8132], data_sandy$Group, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#           pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 Control vs FC  1 0.2534512 1.6253708 0.13981239   0.075     0.1500    
# 2 Control vs PL  1 0.3504074 2.2888681 0.18625540   0.019     0.0990    
# 3 Control vs CM  1 0.2661113 1.8423121 0.15557030   0.033     0.0990    
# 4      FC vs PL  1 0.2673360 1.3026667 0.11525304   0.191     0.2865    
# 5      FC vs CM  1 0.2326521 1.1835326 0.10582816   0.283     0.3396    
# 6      PL vs CM  1 0.1773930 0.9156624 0.08388519   0.391     0.3910   

# Nao houve diferenca significativa entre os grupos para o solo arenoso


pairwise.adonis(data_sandy[5:8132], data_sandy$Type, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr", perm = 999)
#                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1      Unfert vs Fert  1 0.2406765 3.458677 0.4637119   0.100      0.120    
# 2     Unfert vs Fresh  1 0.4674643 2.199092 0.1802668   0.038      0.076    
# 3 Unfert vs Composted  1 0.3317466 3.265370 0.2461575   0.005      0.015   .
# 4       Fert vs Fresh  1 0.2905665 1.319015 0.1165309   0.178      0.178    
# 5   Fert vs Composted  1 0.1597716 1.461580 0.1275200   0.086      0.120    
# 6  Fresh vs Composted  1 0.4887159 2.659214 0.1425148   0.003      0.015   .

# Houve pequenas diferencas significativas entre o Unfert e o Compost; e entre o Fresh e o Compost.

# ------------------------------------------------------------------------------

# Ordenacao NMDS
matrix.nmds<-metaMDS(data_clay[5:8132], k=2, distance = 'bray', noshare=FALSE, autotransform=TRUE) # k=2 transforma em nmds de duas dimensoes, autotransform = true para dados de sequencias
matrix.nmds$stress #0.061 (great!)
# A rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation. 

matrix.nmds<-metaMDS(data_sandy[5:8132], k=2, distance = 'bray', noshare=FALSE, autotransform=TRUE) # k=2 transforma em nmds de duas dimensoes, autotransform = true para dados de sequencias
matrix.nmds$stress #0.107 (great!)




NMDS1<-matrix.nmds$points[,1]
NMDS2<-matrix.nmds$points[,2]
NMDS=data.frame(NMDS1=NMDS1, NMDS2=NMDS2, Group=data$Group, Type=data$Type, Soil=data$Soil)

# Plotar a NMDS
ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=Group, shape=Type), size=4,  alpha=0.8) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  #stat_ellipse(aes(color = Group), type = "norm", linetype = 2, alpha = 0.5)+
  scale_colour_manual(values=c3) +
  scale_shape_manual(values=c(15, 16, 17, 18)) +
  theme_bw() +
  ggtitle("Beta diversity")+
  annotate("text", x = 0.8, y = 1.8, hjust = 1 , label = "
  stress = 0.047
  
  PERMANOVA: p < 0.001
  Group: F = 2.90 | R2 = 0.08
  Type: F = 4.178 | R2 = 0.09
  Soil: F = 33.22 | R2 = 0.33
  G:T:S: F = 1.69 | R2 = 0.03", size = 2.5)+
  theme(legend.position="right",
        legend.text = element_text(size=10),
        axis.text.x= element_text(angle= 0, size= 10),
        axis.text.y= element_text(size= 10),
        axis.title= element_text(size= 10, face= "bold"))

### Arrange plots on one page ###
ggarrange(NMDS_16S, NMDS_ITS, 
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 2, nrow = 1)


#  save plot with 600 dpi resolution
dev.print(tiff, "NMDS_16S_ITS.tiff", compression = "lzw", res=600, height=12, width=20, units="cm")
