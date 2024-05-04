#Leewis et al 2024
#Long-term legacy of phytoremediation on plant succession and soil microbial communities in petroleum-contaminated sub-Arctic soils
#Coding

library(vegan)
library(ggplot2)
  theme_set(theme_bw())
library(phyloseq); packageVersion("phyloseq")
library(dplyr)
library(car)  #companion to applied regression, necessary for much 2way anovas
library(rcompanion) #functions to support extension programs
library(ggpubr) #for plotting data to see if normal (etc)
library(RColorBrewer) 
library(Hmisc)
library(corrplot)

setwd("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/2024_org_analysis")

#colours based on trt (or trt_2). Darker = CO:
co.col.list <- c("#1B9E77", "#66A61E", "#D95F02","#E6AB02", "#A6761D", "#7570B3", "darkviolet") 
de.col.list <- c("#66C2A5", "#A6D854", "#FC8D62", "#FFD92F", "#E5C494", "#8DA0CB", "mediumpurple1") 
#co then de
  all.col.list<- c("#1B9E77", "#66A61E", "#D95F02","#E6AB02", "#A6761D", "#7570B3", "darkviolet", "#66C2A5", "#A6D854", "#FC8D62", "#FFD92F", "#E5C494", "#8DA0CB", "mediumpurple1") 

#trt SHAPE list cntrl1, cntrl2, fert, plnt1, plnt1.fert, plnt2, plnt2.fert
  #scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))
  
#### Soil Chemistry measured on all sub-plots ####
all.dat <- read.table("metadata_2024.txt", row.names=1, header=T)
  de.dat <- subset(all.dat, Contam =="Diesel")
  co.dat <- subset(all.dat, Contam =="CrudeOil")

tph.dat <- read.table("tph_stats.txt", row.names=1, header=T)
  de.tph <- subset(tph.dat, Contam =="Diesel")
  co.tph <- subset(tph.dat, Contam =="CrudeOil")
##Normality Testing - Examine Individual Elements
  #normality test (normal if p>0.05)
  shapiro.test(de.dat$Moo.tot) 
  #test for homogeneity of variance in >2 groups (normal if p>0.05)
  leveneTest(Moo.tot ~ trt*Contam, data = de.dat) #normal
 #Visual Methods to assess normality: 
  #qqPlots (does data follow reference line?)
  ggqqplot(de.dat$Moo.tot)
  #density plot (i.e. is it bell shaped?)
  ggdensity(co.dat$Moo.tot) #uses library(ggpubr)

#Normal data: 
  #ANOVA
  res.aov1 <- aov(Moo.tot ~ Trtmt, data = de.dat)
  summary(res.aov1)
  #post-hoc test
  TukeyHSD(res.aov1)#, which = "Plant")
  
  #NON-normal data: Kruskal-Wallace Test
  kruskal.test(Moo.tot ~ Trtmt, data = co.dat)
  kruskal.test(tphpct ~ Trtmt, data = de.dat)
  
#Mann-Whitney U/Wilcox Rank Sum for post-hoc pair-wise comparisons
pairwise.wilcox.test(co.dat$Moo.tot, co.dat$Trtmt,
                     p.adjust.method = "fdr") 
pairwise.wilcox.test(de.dat$tphpct, de.dat$Trtmt,
                     p.adjust.method = "fdr") 

##### Figure 2 -  Plant % Cover Plot ########
  #A# bare.co, grass.co, tree.co, forb.co
  #B# bare.de, grass.de, tree.de, forb.de
bare.co  <- ggplot(co.dat, aes(x= Trtmt, bare_pct, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y="Percent Cover (est. %)", x=element_blank(), 
       subtitle = "(a) Bare Ground, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)
grass.co  <- ggplot(co.dat, aes(x= Trtmt, grass_pct, fill=Contam))+
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(b) Grass, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)
tree.co  <- ggplot(co.dat, aes(x= Trtmt, tree_pct, fill=Contam))+
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(c) Tree, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)
forb.co  <- ggplot(co.dat, aes(x= Trtmt, forb_pct, fill=Contam))+
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(d) Forb, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)
forb.co
  # DE
bare.de  <- ggplot(de.dat, aes(x= Trtmt, bare_pct, fill=Contam))+
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray80"))+ #de:  
  labs(y="Percent Cover (est. %)", x=element_blank(), 
       subtitle = "(e) Bare Ground, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)
grass.de  <- ggplot(de.dat, aes(x= Trtmt, grass_pct, fill=Contam))+
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray80"))+
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(f) Grass, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)
tree.de  <- ggplot(de.dat, aes(x= Trtmt, tree_pct, fill=Contam))+
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray80"))+
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(g) Tree, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)
forb.de  <- ggplot(de.dat, aes(x= Trtmt, forb_pct, fill=Contam))+
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray80"))+
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(h) Forb, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,110)

figure.pct <- ggarrange(bare.co, grass.co, tree.co, forb.co,
                        bare.de, grass.de, tree.de, forb.de,
                        ncol = 4, nrow = 2)
  figure.pct
  
  
# Figure MPN and PLFA box plots 
    # Only total microbial biomass, bacterial biomass is basically the same, not necessary to show twice with the stacked bar charts made in excel
MPN_plot  <- ggplot(all.dat, aes(x= Trtmt, MPN.Log, fill=Contam))+ 
    geom_boxplot(varwidth=T, outlier.size=-1) +
    geom_jitter(size = 3) +
    theme_bw(base_size=16) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle=45, vjust=1, hjust=1),
          legend.position="none")+
    scale_fill_manual(values = c("gray32", "gray80"))+
    labs(y="DDM (log MPN/g soil)", x=element_blank(),
         subtitle = "(a) DDM") +
    facet_wrap(~Contam, scales="free_x")+
    scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))#+
MooPLFA_plot  <- ggplot(all.dat, aes(x= Trtmt, Moo.tot, fill=Contam))+ 
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) +
  theme_bw(base_size=16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position="none")+
  scale_fill_manual(values = c("gray32", "gray80"))+
  labs(y="PLFA (nmol/g soil)", x=element_blank(),
       subtitle = "(b) Total Microbial Biomass") +
  facet_wrap(~Contam, scales="free_x")+
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))#+

figure.num <- ggarrange(MPN_plot, MooPLFA_plot,
                        ncol = 1, nrow = 2)
figure.num

#### COMBINED MPN VS TREES & TPH Correlation Plot FOR FIGURE 3? #### 
de.dat
co.dat

### Correlation Stats ##### #Be sure to test normality with Shapiro Test, 
  #if not normal use Kendall Rank Correlation Test)
#shapiro test (normal if p>0.05)
shapiro.test(co.dat$TPH_2011)
shapiro.test(de.dat$TPH_2011)

#Individual Correlations #method=c("pearson", "kendall", "spearman")
  #bare_pct	grass_pct	tree_pct	forb_pct
cor.test(co.dat$Moo.tot, co.dat$tot.pctf, method="kendall")
cor.test(de.dat$Moo.tot, de.dat$tot.pctf, method="kendall")

#graphing
co.MPN.TPH.corr.plot <- ggplot(co.dat, aes(TPH_2011, MPN.Log, col=Trtmt))+ 
  geom_jitter(aes(shape=Trtmt), size = 3) + #col=trt, 
  geom_smooth(method="glm", se=F, color = "black")+
  scale_colour_manual(values= co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12))+
  guides(col = FALSE)+
  labs(x = "TPH (ppm)", y = "DDM (log MPN/g soil)", 
       shape ="Treatment", title="(a) Crude Oil")

de.MPN.TPH.corr.plot <- ggplot(de.dat, aes(TPH_2011, MPN.Log, col=Trtmt))+ 
  geom_jitter(aes(shape=Trtmt), size = 3) + #col=trt, 
  geom_smooth(method="glm", se=F, color = "black")+
  scale_colour_manual(values= co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12))+
  guides(col = FALSE)+
  labs(x = "TPH (ppm)", y = "DDM (log MPN/g soil)", 
       color ="Treatment", title="(b) Diesel")

co.corplot <- ggplot(co.dat, aes(tree_pct, MPN.Log, col=Trtmt))+ 
  geom_jitter(aes(shape=Trtmt), size = 3) + #col=trt, 
  geom_smooth(method="glm", se=F, color = "black")+
  scale_colour_manual(values= co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12))+
  guides(col = FALSE)+
  labs(x = "Tree Cover (%)", y = "DDM (log MPN/g soil)", 
       shape="Treatment", title="(c) Crude Oil") #title = "Crude Oil DDM vs % trees"
co.corplot

de.corplot <- ggplot(de.dat, aes(tree_pct, MPN.Log, col=Trtmt))+ 
  geom_jitter(aes(shape=Trtmt), size = 3) + #col=trt, 
  geom_smooth(method="glm", se=F, color = "black")+
  scale_colour_manual(values= co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12))+
  guides(col = FALSE)+
  labs(x = "Tree Cover (%)", y = "DDM (log MPN/g soil)", 
       shape="Treatment", title="(d) Diesel") #, title = "Diesel DDM vs % trees"
de.corplot

figureMPN <- ggarrange(co.MPN.TPH.corr.plot, de.MPN.TPH.corr.plot, co.corplot, de.corplot,
                       #labels = c("A)", "B)", "C)", "D)"),
                       nrow = 2, ncol=2,
                       common.legend = TRUE, legend = "right")
figureMPN # just save as PDF 6 x 8

#### PLFA vs Plant Cover
co.plant.corplot <- ggplot(co.dat, aes(tot.pctf, Moo.tot, col=Trtmt))+ 
  geom_jitter(aes(shape=Trtmt), size = 3) + #col=trt, 
  geom_smooth(method="glm", se=F, color = "black")+
  scale_colour_manual(values= co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12))+
  guides(col = FALSE)+
  labs(x = "Total Plant Cover (est. %)", y = "Microbial Biomass (PLFA nmol/g soil)", 
       shape="Treatment", title="(a) Crude Oil") #, title = "Diesel DDM vs % trees"
co.plant.corplot

de.plant.corplot <- ggplot(de.dat, aes(tot.pctf, Moo.tot, col=Trtmt))+ 
  geom_jitter(aes(shape=Trtmt), size = 3) + #col=trt, 
  geom_smooth(method="glm", se=F, color = "black")+
  scale_colour_manual(values= co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=12))+
  guides(col = FALSE)+
  labs(x = "Total Plant Cover (est. %)", y = "Microbial Biomass (PLFA nmol/g soil)", 
       shape="Treatment", title="(b) Diesel") #, title = "Diesel DDM vs % trees"
de.plant.corplot

figurepctcover <- ggarrange(co.plant.corplot, de.plant.corplot,
                       #labels = c("A)", "B)", "C)", "D)"),
                       nrow = 2, ncol=1,
                       common.legend = TRUE, legend = "right")
figurepctcover #save as 6x8 portrait

#######
##### Figure Supplemental -  PLFAs by treatment [not sure this will be included in final paper] ########
#G.pos	G.neg	Actino	Fungi	Protoz	Uncl
#CO
G.pos.co  <- ggplot(co.dat, aes(x= Trtmt, G.pos, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y="PLFA Relative Abundance (%)", x=element_blank(), 
       subtitle = "(a) G+, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,5)
G.neg.co  <- ggplot(co.dat, aes(x= Trtmt, G.neg, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(b) G-, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,8)
Actino.co  <- ggplot(co.dat, aes(x= Trtmt, Actino, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(c) Actino, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,2.5)
Fungi.co  <- ggplot(co.dat, aes(x= Trtmt, Fungi, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(d) Fungi, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,14)
Protoz.co  <- ggplot(co.dat, aes(x= Trtmt, Protoz, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(e) Proz, Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,2)
Uncl.co  <- ggplot(co.dat, aes(x= Trtmt, Uncl, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(e) Uncl., Crude Oil" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,0.5)

##DE
G.pos.de  <- ggplot(de.dat, aes(x= Trtmt, G.pos, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y="PLFA Relative Abundance (%)", x=element_blank(), 
       subtitle = "(a) G+, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,5)
G.neg.de  <- ggplot(de.dat, aes(x= Trtmt, G.neg, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(b) G-, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,8)
Actino.de  <- ggplot(de.dat, aes(x= Trtmt, Actino, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(c) Actino, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,2.5)
Fungi.de  <- ggplot(de.dat, aes(x= Trtmt, Fungi, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(d) Fungi, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,14)
Protoz.de  <- ggplot(de.dat, aes(x= Trtmt, Protoz, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(e) Proz, Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,2)
Uncl.de  <- ggplot(de.dat, aes(x= Trtmt, Uncl, fill=Contam))+ #Pub.Name
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) + #aes(shape=Trtmt), 
  theme_bw(base_size=11) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(size = 12),
        legend.position="none")+
  scale_fill_manual(values = c("gray32"))+ #de:  "gray80"
  labs(y=element_blank(), x=element_blank(), 
       subtitle = "(e) Uncl., Diesel" ) +
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))+
  ylim(0,0.5)


figurePLFAindivid <- ggarrange(G.pos.co,	G.neg.co,	Actino.co,	Fungi.co,	Protoz.co, Uncl.co,
                               G.pos.de,	G.neg.de,	Actino.de,	Fungi.de,	Protoz.de, Uncl.de,
                            #labels = c("A)", "B)", "C)", "D)"),
                            nrow = 2, ncol=6,
                            legend = "none")
figurePLFAindivid ## Decided not to include in paper, no strong relationship. 



#######################################################################
#### Microbiology Measured only on A,B,C therefore subset of sites ####
#######################################################################

#LOAD sequence data to phyloseq##
otu_mat<- read.table( "/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/2021_16S_phyloseq/otu_table.txt",header=T)
tax_mat<- read.table( "/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/2021_16S_phyloseq/taxonomy.txt",header=T)
samples_df <- read.table( "/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/2021_16S_phyloseq/metadata.txt",header=T)

#define the first column as OTUs or Sample IDs
row.names(otu_mat) <- otu_mat$OTU
otu_mat <- otu_mat %>% select (-OTU) 

row.names(tax_mat) <- tax_mat$OTU
tax_mat <- tax_mat %>% select (-OTU) 

row.names(samples_df) <- samples_df$sampleID
samples_df <- samples_df %>% select (-sampleID) 

#Transform into matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

FLR <- phyloseq(OTU, TAX, samples)
FLR

#Save PhyloSeq Object as RDS
saveRDS(FLR, file="FLR.RDS")

#subset the TWO contaminants
CO_flr <- subset_samples(FLR, Contam =="CrudeOil")
DE_flr <- subset_samples(FLR, Contam =="Diesel")

### AlphaDiversity (note: data rarified in mothur)
plot_richness(FLR)
plot_richness(FLR, measures=c("Chao1", "Shannon"))

plot_richness(CO_flr, x="Pub.Name.contam", measures=c("Chao1"), color="Pub.Name") +
  facet_wrap(~Contam, scales="free_x")

p <- plot_richness(FLR, x="Trtmt", measures=c("Shannon"), color="trt.contam") + 
  facet_wrap(~Contam, scales="free_x") +
  geom_point(size=4)+
  scale_colour_manual(values=trt.col.list)+
  labs(x="Original Treatment", y="Shannon Diveristy (H')")+
  facet_grid(rows = vars(Contam))
p


OTU.plot  <- ggplot(samples_df, aes(x= Trtmt, Observed, fill=Contam))+ 
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) +
  theme_bw(base_size=16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position="none")+
  scale_fill_manual(values = c("gray32", "gray80"))+
  labs(y="No. OTUs observed", x=element_blank(),
       subtitle = "(a) Num. OTU") +
  facet_wrap(~Contam, scales="free_x")+
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))#+
Shann.plot  <- ggplot(samples_df, aes(x= Trtmt, Shannon, fill=Contam))+ 
  geom_boxplot(varwidth=T, outlier.size=-1) +
  geom_jitter(size = 3) +
  theme_bw(base_size=16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position="none")+
  scale_fill_manual(values = c("gray32", "gray80"))+
  labs(y="Shannon Diversity Index", x=element_blank(),
       subtitle = "(b) Shannon ") +
  facet_wrap(~Contam, scales="free_x")+
  scale_x_discrete(limits = c("c1", "c2", "p1", "p2", "f", "p1f", "p2f"))#+

figure.alpha <- ggarrange(OTU.plot, Shann.plot,
                        ncol = 1, nrow = 2)
figure.alpha #save


de.env.al <- subset(samples_df, Contam =="Diesel")
co.env.al <- subset(samples_df, Contam =="CrudeOil")

##Normality Testing - Examine Individual Elements
#normality test (normal if p>0.05)
shapiro.test(co.env.al$Observed) 
#test for homogeneity of variance in >2 groups (normal if p>0.05)
leveneTest(Moo.tot ~ trt*Contam, data = de.dat) #normal
#Visual Methods to assess normality: 
#qqPlots (does data follow reference line?)
ggqqplot(de.dat$Moo.tot)
#density plot (i.e. is it bell shaped?)
ggdensity(co.dat$Moo.tot) #uses library(ggpubr)

#Normal data: 
#ANOVA
res.aov1 <- aov(Shannon ~ Trtmt, data = co.env.al)
summary(res.aov1)
#post-hoc test
TukeyHSD(res.aov1)#, which = "P



#NON-normal data: Kruskal-Wallace Test
kruskal.test(Observed ~ Contam, data = samples_df)
kruskal.test(Observed ~ Trtmt, data = co.env.al)
kruskal.test(Observed ~ Trtmt, data = de.env.al)

#Mann-Whitney U/Wilcox Rank Sum for post-hoc pair-wise comparisons
pairwise.wilcox.test(co.env.al$Moo.tot, co.env.al$Trtmt,
                     p.adjust.method = "fdr") 
pairwise.wilcox.test(de.env.al$tphpct, de.env.al$Trtmt,
                     p.adjust.method = "fdr") 



#Write alphaDiv Table
richness_tab = estimate_richness(FLR)
write.csv(richness_tab, "/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/2021_16S_phyloseq/AlphaDiv_tab.csv")
#end alpha diversity ###############################################

###############################################
################### heatmap ###################
###############################################
p6 <- ggplot(filt_major_taxa_proportions_tab_for_plot.g2, 
             aes(char, Major_Taxa, fill= Proportion)) + 
  geom_tile(color = "white",lwd = 0.5,linetype = 1) + #white lines between each box
  facet_wrap(~site, scales="free_x")+
  labs(fill="Mean\nRelative\nAbund. (%)") +      #<\n insert line break
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,   #align to middle
  )+ 
  scale_fill_gradient(low = "lightblue", high = "red") +
  scale_x_discrete(limits = c("c1", "c2", "p1", 
                              "p2", "f", "p1f", 
                              "p2f"))
p6 # save as portrait 8x11
#end heatmap ###

###############################################
################### Ordinations ###############
###############################################

#### Ordinations ########
setwd("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2")
#PLFA data is in nmol/g#
all.PLFA <- read.table("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/PLFA_individ.txt", row.names=1, header=T)
dePLFA <- read.table("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/PLFA_de_individ.txt", row.names=1, header=T)
coPLFA <- read.table("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/PLFA_co_individ.txt", row.names=1, header=T)

all.16S <- read.table("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/OTU.txt", row.names=1, header=T) #updated august 2021
de16S <- read.table("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/OTU_De.txt", row.names=1, header=T)
co16S <- read.table("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/OTU_CO.txt", row.names=1, header=T)


all.16S.env <- read.table("/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/env.txt", row.names=1, header=T)
de.env <- subset(all.16S.env, Contam =="Diesel")
co.env <- subset(all.16S.env, Contam =="CrudeOil")

#all.veg <- read.table("allveg.txt", row.names=1, header=T)
#de.veg <- subset(all.veg, Contam =="Diesel")x
#co.veg <- subset(all.veg, Contam =="CrudeOil")

#####November 2021: remove highly correlated variables?
#rm the highly correlated variables, leaving only one of them: 
de.env <-subset(de.env, select = -c(Contam,	trt_2,	fert, P, K, C, Silt, Clay, Tree.num, TPH95, TPH96) )
co.env <-subset(co.env, select = -c(Contam,	trt_2,	fert,P, K, C, Silt, Clay, Tree.num, TPH95, TPH96) )

## CO and DE together

NMDS.co.de.PLFA <- metaMDS(all.PLFA, distance="bray", k=2)
stressplot(NMDS.co.de.PLFA)
data.scores = as.data.frame(scores(NMDS.co.de.PLFA))
data.scores$trt = all.16S.env$trt
data.scores$Contam = all.16S.env$Contam
co.de.PLFA = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores,  aes(shape = Contam, colour=trt), size = 4) + 
  labs(shape = "Treatment", colour = "Soil Type", title = "Crude Oil & Diesel, PLFAs") +
  scale_colour_manual(values=all.col.list)+
  #scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
co.de.PLFA

NMDS.co.de.16S <- metaMDS(all.16S, distance="bray", k=2)
stressplot(NMDS.co.de.16S)
data.scores = as.data.frame(scores(NMDS.co.de.16S))
data.scores$trt = all.16S.env$trt
data.scores$Contam = all.16S.env$Contam
co.de.16S = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores,  aes(shape = Contam, colour=trt), size = 4) + 
  labs(shape = "Treatment", colour = "Soil Type", title = "Crude Oil & Diesel, 16S rRNA") +
  scale_colour_manual(values=all.col.list)+
  #scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
co.de.16S

figure.all.ord <- ggarrange(co.de.PLFA, co.de.16S,
                    labels = c("(a)", "(b)"),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "right")
figure.all.ord #save PDF as portrait 8x11

#1 of 4########### DE Ordinations 16S ################
#NDMS of 16S DE data 
NMDS.de.16S <- metaMDS(de16S, distance="bray", k=2)

#check the stress/fit: large scatter indicates bad fit, but this seems ok.
stressplot(NMDS.de.16S)
#calculate stress: non-metric fit R2=1-S2: 
#(sqrt(0.988)-1)/-1 = 0.006018
#basic NMDS plot
plot.NMDS.de.16S <- plot(NMDS.de.16S, display="sites", type="t")
#add ENV vectors
ef <- envfit(NMDS.de.16S, de.env, permu = 999)
ef
plot(ef, p.max = 0.05)

#write to a file, add [, append=TRUE] if want to add to the end of another file
#sink('De_16S_vectors.txt') 
#ef
#sink()

#Extract the significant vectors to use with ggplot. Example based on <https://jkzorz.github.io/2020/04/04/NMDS-extras.html>
#for NMDS output use the following code to extract the sample coordinates in ordination space
data.scores = as.data.frame(scores(NMDS.de.16S))
data.scores$trt = de.env$trt
#extract ALL environmental vectors
en_coord_cont = as.data.frame(scores(ef, "vectors")) * ordiArrowMul(ef)
en_coord_cat = as.data.frame(scores(ef, "factors")) * ordiArrowMul(ef)

#extract SIGNIFICANT environmental vectors only
env.scores.MDS <- as.data.frame(scores(ef, display = "vectors"))* ordiArrowMul(ef) #extracts relevant scores from envifit
env.scores.MDS <- cbind(env.scores.MDS, env.variables = rownames(env.scores.MDS)) #and then gives them their names
env.scores.MDS <- cbind(env.scores.MDS, pval = ef$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.MDS, pval<=0.05) #subset data to show variables significant at 0.05
head(env.scores.MDS)

#Now plot the Ordination with GGplot 
de16S = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores,  aes(shape = trt, colour=trt), size = 4) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(shape = "Treatment", title = "Diesel, 16S rRNA") +
  scale_colour_manual(values=co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))
de16S

#add environmental factors plot only those significant
de16S <- de16S +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scrs, size =1, alpha = 0.5, colour = "grey8") +
  geom_text(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2), colour = "grey8", 
            fontface = "bold", label = row.names(sig.env.scrs))+
  guides(color=FALSE)
de16S  
##2 of 4########### CO Ordinations 16S ################
NMDS.co.16S <- metaMDS(co16S, distance="bray", k=2)
stressplot(NMDS.co.16S)
#(sqrt(0.995)-1)/-1 = 0.002503
plot.NMDS.co.16S <- plot(NMDS.co.16S, display="sites", type="t")
ef <- envfit(NMDS.co.16S, co.env, permu = 999)
ef
plot(ef, p.max = 0.05)
#sink('CO_16S_vectors.txt') 
#ef
#sink()

data.scores = as.data.frame(scores(NMDS.co.16S))
data.scores$trt = co.env$trt
en_coord_cont = as.data.frame(scores(ef, "vectors")) * ordiArrowMul(ef)
en_coord_cat = as.data.frame(scores(ef, "factors")) * ordiArrowMul(ef)
env.scores.MDS <- as.data.frame(scores(ef, display = "vectors"))* ordiArrowMul(ef) #extracts relevant scores from envifit
env.scores.MDS <- cbind(env.scores.MDS, env.variables = rownames(env.scores.MDS)) #and then gives them their names
env.scores.MDS <- cbind(env.scores.MDS, pval = ef$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.MDS, pval<=0.05) #subset data to show variables significant at 0.05
head(env.scores.MDS)

co16S = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores,  aes(shape = trt, colour=trt), size = 4) + 
  labs(shape = "Treatment", title = "Crude Oil, 16S rRNA") +
  scale_colour_manual(values=co.col.list)+
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scrs, size =1, alpha = 0.5, colour = "grey8") +
  geom_text(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2), colour = "grey8", 
            fontface = "bold", label = row.names(sig.env.scrs))+
  guides(color=FALSE)
co16S
##3 of 4########### CO Ordinations PLFA ################
NMDS.co.PLFA <- metaMDS(coPLFA, distance="bray", k=2)
stressplot(NMDS.co.PLFA)
#(sqrt(0.999)-1)/-1 = 0.00050
plot.NMDS.co.PLFA <- plot(NMDS.co.PLFA, display="sites", type="t")
ef <- envfit(NMDS.co.PLFA, co.env, permu = 999)
ef
plot(ef, p.max = 0.05)
#sink('CO_PLFA_vectors.txt') 
#ef
#sink()
data.scores = as.data.frame(scores(NMDS.co.PLFA))
data.scores$trt = co.env$trt
en_coord_cont = as.data.frame(scores(ef, "vectors")) * ordiArrowMul(ef)
en_coord_cat = as.data.frame(scores(ef, "factors")) * ordiArrowMul(ef)
env.scores.MDS <- as.data.frame(scores(ef, display = "vectors"))* ordiArrowMul(ef) #extracts relevant scores from envifit
env.scores.MDS <- cbind(env.scores.MDS, env.variables = rownames(env.scores.MDS)) #and then gives them their names
env.scores.MDS <- cbind(env.scores.MDS, pval = ef$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.MDS, pval<=0.05) #subset data to show variables significant at 0.05
head(env.scores.MDS)

coPLFA <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores,  aes(shape = trt, colour=trt), size = 4) + 
  labs(shape = "Treatment", title = "Crude Oil, PLFAs") +
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  scale_colour_manual(values=co.col.list)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scrs, size =1, alpha = 0.5, colour = "grey8") +
  geom_text(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2), colour = "grey8", 
            fontface = "bold", label = row.names(sig.env.scrs))+
  guides(color=FALSE)

# 4 of 4 ########## DE Ordinations PLFA ################
NMDS.de.PLFA <- metaMDS(dePLFA, distance="bray", k=2)
stressplot(NMDS.de.PLFA)
#(sqrt(0.996)-1)/-1 = 0.002002
plot.NMDS.de.PLFA <- plot(NMDS.de.PLFA, display="sites", type="t")
ef <- envfit(NMDS.de.PLFA, co.env, permu = 999)
ef
plot(ef, p.max = 0.05)
#sink('DE_PLFA_vectors.txt') 
#ef
#sink()
data.scores = as.data.frame(scores(NMDS.de.PLFA))
data.scores$trt = co.env$trt
en_coord_cont = as.data.frame(scores(ef, "vectors")) * ordiArrowMul(ef)
en_coord_cat = as.data.frame(scores(ef, "factors")) * ordiArrowMul(ef)
env.scores.MDS <- as.data.frame(scores(ef, display = "vectors"))* ordiArrowMul(ef) #extracts relevant scores from envifit
env.scores.MDS <- cbind(env.scores.MDS, env.variables = rownames(env.scores.MDS)) #and then gives them their names
env.scores.MDS <- cbind(env.scores.MDS, pval = ef$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores.MDS, pval<=0.05) #subset data to show variables significant at 0.05
head(env.scores.MDS)

dePLFA <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores,  aes(shape = trt, colour=trt), size = 4) + 
  labs(shape = "Treatment", title = "Diesel, PLFAs") +
  scale_shape_manual(values= c(12,13,18, 2, 17, 0, 15))+
  scale_colour_manual(values=co.col.list)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = sig.env.scrs, size =1, alpha = 0.5, colour = "grey8") +
  geom_text(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2), colour = "grey8", 
            fontface = "bold", label = row.names(sig.env.scrs))+
  guides(color=FALSE)

############ Arrange into one panel ###

figure <- ggarrange(coPLFA, dePLFA,co16S, de16S,
                    labels = c("(a)", "(b)", "(c)", "(d)"),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right")
figure

ggexport(figure, filename = "figure3_202212.pdf")



###############################################
####### Muiltivariate Stats ###################
###############################################

dePLFA.dist <- vegdist(dePLFA, method="bray")
coPLFA.dist <- vegdist(coPLFA, method="bray")
allPLFA.dist <- vegdist(all.PLFA, method="bray")

de16S.dist <- vegdist(de16S, method="bray")
co16S.dist <- vegdist(co16S, method="bray")
all16S.dist <- vegdist(all.16S, method="bray")

## ANOSIM to see if groupings plot differently ##
de.ano <-anosim(de16S.dist, de.env$trt, permutations = 999)
de.ano
plot(de.ano)
co.ano <-anosim(co16S.dist, co.env$trt, permutations = 999)
co.ano
plot(co.ano)

## MRPP# non parametric procedure for testing the hypothesis of no difference between two or more groups of entities. 
#dune.mrpp <- with(de.env, mrpp(dePLFA, Trtmt))

de16S.mrpp <-mrpp(de16S.dist, de.env$trt)
de16S.mrpp

co16S.mrpp <-mrpp(co16S.dist, co.env$trt)
co16S.mrpp

allveg <- all.veg[,10:ncol(all.veg)]
all.veg.dist <- vegdist(allveg, method="bray", na.rm=T)
allveg.mrpp <-mrpp(all.veg.dist, all.veg$Contam) 

deveg.mrpp <-mrpp(de.veg.dist, de.env$trt)
deveg.mrpp

coveg.mrpp <-mrpp(co.veg.dist, co.env$trt)
coveg.mrpp



# if want to test the individual comparisons, then it's 7x7=49 comparisons. 
#need to subset first, then run comparison
cntfert <- subset(de.env, trt == 'plnt_1' | trt == 'fert') #subset 
#depthSM.mrpp <- mrpp(dePLFA.dist, cntfert$trt)
#depthSM.mrpp<- mrpp(dePLFA.dist, de.env[de.env$trt == "plnt_1" | de.env$trt == "fert"])  
#depthSM.mrpp


## Mantel test. Compare two distance/similarity matrices that were obtained independently of each other. (e.g. are these two distance matrices related? e.g. ICP-OES vs community structure)
# Compare Veg vs. Microbial Community
#set up veg data from the metadata (i.e. only numeric values which start at column 8) then turn into a distance matrix
devegdata <- de.veg[,9:ncol(de.veg)]
de.veg.dist <- vegdist(devegdata, method="bray", na.rm=T)
covegdata <- co.veg[,9:ncol(co.veg)]
co.veg.dist <- vegdist(covegdata, method="bray", na.rm=T)

mantel(de16S.dist, de.veg.dist, method = "pearson", permutations = 999)    
mantel(all.dist, veg.dist, method = "pearson", permutations = 999)    

mantel(co.veg.dist, de.veg.dist, method = "pearson", permutations = 999)    

### permanova
##perm.blocked <- adonis(all.dist ~ Ter_Site/Y_cm, data=samples_df, permutations=999)
trt.coveg <- adonis(co.veg.dist ~ trt, data=co.env, permutations=999)
fert.coveg <- adonis(co.veg.dist ~ fert, data=co.env, permutations=999)
trt.deveg <- adonis(de.veg.dist ~ trt, data=de.env, permutations=999)
fert.deveg <- adonis(de.veg.dist ~ fert, data=de.env, permutations=999)
pH.coveg <- adonis(co.veg.dist ~ pH*NO3, data=co.env, permutations=999)
NO3.coveg <- adonis(co.veg.dist ~ NO3, data=co.env, permutations=999)

trt.de16S <- adonis(de16S.dist ~ trt, data=de.env, permutations=999)
fert.de16S <- adonis(de16S.dist ~ fert , data=de.env, permutations=999)
#tree.de16S <- adonis(de16S.dist ~ Tree.num , data=de.env, permutations=999)
#veg.de16S <- adonis(de16S.dist ~ Veg.num , data=de.env, permutations=999)
pH.de16S <- adonis(de16S.dist ~ pH , data=de.env, permutations=999)
NO3.de16S <- adonis(de16S.dist ~ NO3 , data=de.env, permutations=999)
TPH.de16S <- adonis(de16S.dist ~ TPH , data=de.env, permutations=999)
P.de16S <- adonis(de16S.dist ~ P , data=de.env, permutations=999)
K.de16S <- adonis(de16S.dist ~ K , data=de.env, permutations=999)
CEC.de16S <- adonis(de16S.dist ~ CEC , data=de.env, permutations=999)
C.de16S <- adonis(de16S.dist ~ C , data=de.env, permutations=999)
Sand.de16S <- adonis(de16S.dist ~ Sand , data=de.env, permutations=999)
Silt.de16S <- adonis(de16S.dist ~ Silt , data=de.env, permutations=999)
Clay.de16S <- adonis(de16S.dist ~ Clay , data=de.env, permutations=999)
treepct.dePLFA <- adonis(dePLFA.dist ~ tree.pct , data=de.env, permutations=999)
grasspct.dePLFA <- adonis(dePLFA.dist ~ grass.pct , data=de.env, permutations=999)
forbpct.dePLFA <- adonis(dePLFA.dist ~ forb.pct , data=de.env, permutations=999)
barepct.dePLFA <- adonis(dePLFA.dist ~ bare.pct , data=de.env, permutations=999)

trt.co16S <- adonis(co16S.dist ~ trt, data=co.env, permutations=999)
fert.co16S <- adonis(co16S.dist ~ fert , data=co.env, permutations=999)
tree.co16S <- adonis(co16S.dist ~ Tree.num , data=co.env, permutations=999)
veg.co16S <- adonis(co16S.dist ~ Veg.num , data=co.env, permutations=999)
pH.co16S <- adonis(co16S.dist ~ pH , data=co.env, permutations=999)
NO3.co16S <- adonis(co16S.dist ~ NO3 , data=co.env, permutations=999)
TPH.co16S <- adonis(co16S.dist ~ TPH , data=co.env, permutations=999)
P.co16S <- adonis(co16S.dist ~ P , data=co.env, permutations=999)
K.co16S <- adonis(co16S.dist ~ K , data=co.env, permutations=999)
CEC.co16S <- adonis(co16S.dist ~ CEC , data=co.env, permutations=999)
C.co16S <- adonis(co16S.dist ~ C , data=co.env, permutations=999)
Sand.co16S <- adonis(co16S.dist ~ Sand , data=co.env, permutations=999)
Silt.co16S <- adonis(co16S.dist ~ Silt , data=co.env, permutations=999)
Clay.co16S <- adonis(co16S.dist ~ Clay , data=co.env, permutations=999)
treepct.coPLFA <- adonis(coPLFA.dist ~ tree.pct , data=co.env, permutations=999)
grasspct.coPLFA <- adonis(coPLFA.dist ~ grass.pct , data=co.env, permutations=999)
forbpct.coPLFA <- adonis(coPLFA.dist ~ forb.pct , data=co.env, permutations=999)
barepct.coPLFA <- adonis(coPLFA.dist ~ bare.pct , data=co.env, permutations=999)

perm.dePLFA <- adonis(dePLFA.dist ~ C * pH * NO3 , data=de.env, permutations=999)

#p.adjust
pvalues <- read.table("permanovaPadj.txt",header=T)
p.adjust(as.matrix(pvalues),method="BH")


###### End Vegan Multivariate #################################################################



#######################################################################
#### DeSeq Coding, Adapted from Papik et al 2023  ####
#######################################################################
library(tidyverse)
library(phyloseq)
library(here)
library(ggplot2)

## load a phyloseq object (not rarefied) 
FLR <- readRDS(here("FLR.RDS"))
packageVersion("DESeq2")
sample_data(FLR)

#subset the TWO contaminants
CO_flr <- subset_samples(FLR, Contam =="CrudeOil")
DE_flr <- subset_samples(FLR, Contam =="Diesel")

# build a model based on which variables are the most important in structuring the communities (or of your interest?)
diagdds = phyloseq_to_deseq2(CO_flr, ~deseqcomp)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

# Deseq2 normalization
diagdds = DESeq(diagdds, test="Wald", fitType="local")

# Generating a result table
res = results(diagdds, cooksCutoff = FALSE)
# applying a lfc shrink function as instructed in deseq2 manual so log2fold change values are realistic (e.g., not greater than 10)
# defining which factors within a variable you want to compare
res <- lfcShrink(diagdds, contrast = c("deseqcomp", "pf", "f"), type="normal") #res <- lfcShrink(diagdds, contrast = c("Variable1", "Level1", "Level2"), type="normal")
  #pvf; pvpf; pfvf
res.joro.azul <- res
res.joro.azul

# Define significance threshold and filter results based on that
alpha = 0.05
res.joro.azul = res.joro.azul[order(res.joro.azul$padj, na.last=NA), ]
sigtab = res.joro.azul[(res.joro.azul$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))

# Make a graph and export a result table
#negative values represent taxa enriched by level2, positive values represent taxa enriched by level1
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

#write.table(sigtabgen, "co_phyloseq_pvf.csv")
#write.table(sigtabgen, "co_phyloseq_pvpf.csv")
#write.table(sigtabgen, "co_phyloseq_pfvf.csv")

#load data for top15 written as txt

deseq_compile <- read.table( "/Users/mleewis/OneDrive/Manuscripts/2017_FarmersLoop2/2024_org_analysis/deseq_top15.txt", row.names=1, header=T)
co_pvf <- subset(deseq_compile, Comp.Contam =="PvF_co")
co_pvpf<- subset(deseq_compile, Comp.Contam =="PvPF_co")
co_pfvf<- subset(deseq_compile, Comp.Contam =="PFvF_co")
de_pvf<- subset(deseq_compile, Comp.Contam =="PvF_de")
de_pvpf<- subset(deseq_compile, Comp.Contam =="PvPF_de")
de_pfvf<- subset(deseq_compile, Comp.Contam =="PFvF_de")

# co_pvf co_pvpf co_pfvf
# de_pvf de_pvpf de_pfvf

#graphing
  #negative values represent taxa enriched by level2, 
  #positive values represent taxa enriched by level1
    #therefore write caption as level2 vs level 1 (level 2 are neg), data object labeled L1vL2
co_pvf.g <- 
  ggplot(co_pvf, aes(x=Genus, y=log2FoldChange, fill=Contam)) + #, fill=Phylum
  coord_flip() +
  geom_point(size=6) + #, aes(shape=Phylum)
#  scale_fill_manual(values = c("gray32"))+ #, "gray80"
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme_bw(base_size=16) +
  theme(strip.text.x = element_text(size = 20, colour = "black"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + 
  geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x=element_blank(), y="log2FoldChange",
       subtitle = "(a) F vs. P, Crude Oil")

co_pvpf.g<- 
  ggplot(co_pvpf, aes(x=Genus, y=log2FoldChange)) + 
  coord_flip() +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme_bw(base_size=16) +
  theme(strip.text.x = element_text(size = 20, colour = "black"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + 
  geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x=element_blank(), y="log2FoldChange",
       subtitle = "(b) PF vs. P, Crude Oil")

co_pfvf.g<- 
  ggplot(co_pfvf, aes(x=Genus, y=log2FoldChange)) + 
  coord_flip() +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme_bw(base_size=16) +
  theme(strip.text.x = element_text(size = 20, colour = "black"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + 
  geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x=element_blank(), y="log2FoldChange",
       subtitle = "(c) F vs. PF, Crude Oil")
de_pvf.g<- 
  ggplot(de_pvf, aes(x=Genus, y=log2FoldChange)) + 
  coord_flip() +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme_bw(base_size=16) +
  theme(strip.text.x = element_text(size = 20, colour = "black"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + 
  geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x=element_blank(), y="log2FoldChange",
       subtitle = "(d) F vs. P, Diesel")

de_pvpf.g<- 
  ggplot(de_pvpf, aes(x=Genus, y=log2FoldChange)) + 
  coord_flip() +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme_bw(base_size=16) +
  theme(strip.text.x = element_text(size = 20, colour = "black"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + 
  geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x=element_blank(), y="log2FoldChange",
       subtitle = "(e) PF vs. P, Diesel")
de_pfvf.g <- 
  ggplot(de_pfvf, aes(x=Genus, y=log2FoldChange)) + 
  coord_flip() +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme_bw(base_size=16) +
  theme(strip.text.x = element_text(size = 20, colour = "black"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + 
  geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x=element_blank(), y="log2FoldChange",
       subtitle = "(f) F vs. PF, Diesel")


figure <- ggarrange(co_pvf.g, co_pvpf.g, co_pfvf.g,
                    de_pvf.g,de_pvpf.g,de_pfvf.g,
                    ncol = 3, nrow = 2, legend = "none")
                    #common.legend = TRUE, legend = "right")

#want Phylum colour, but in order wanted.Need to work on below still. 
test<- ggplot(deseq_compile, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  coord_flip() +
  geom_point(size=6) + #, aes(shape=Phylum)
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme_bw(base_size=16) +
  theme(strip.text.x = element_text(size = 20, colour = "black"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + 
  geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))+
  labs(x=element_blank(), y="log2FoldChange")+
  #facet_wrap(~Comp.Contam, scales="free_y")
  facet_grid(.~Comp.Contam)
  facet_grid(~factor(Comp.Contam, 
                     levels=c("co_pvf","co_pvpf", "co_pfvf", "de_pvf", "de_pvpf", "de_pfvf")))

