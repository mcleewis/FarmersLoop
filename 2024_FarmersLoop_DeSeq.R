#######################################################################
#### DeSeq Coding, Adapted from Papik et al 2023  ####
#######################################################################
library(tidyverse)
library(phyloseq)
library(here)
library(ggplot2)

## load a phyloseq object(not rarefied) 
#FLR <- readRDS(here("FLR.RDS"))
packageVersion("DESeq2")
sample_data(FLR)

#subset the TWO contaminants
CO_flr <- subset_samples(FLR, Contam =="CrudeOil")
DE_flr <- subset_samples(FLR, Contam =="Diesel")

# build a model based on which variables are the most important in structuring the communities (or of your interest?)
diagdds = phyloseq_to_deseq2(FLR, ~deseqcomp)
diagdds = phyloseq_to_deseq2(FLR)

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
 #######comparisons for Leewis 2024 - pvf; pvpf; pfvf ######
# co_pvf 
    # co_pvpf co_pfvf
# de_pvf de_pvpf de_pfvf

res <- lfcShrink(diagdds, contrast = c("deseqcomp", "p", "f"), type="normal") #res <- lfcShrink(diagdds, contrast = c("Variable1", "Level1", "Level2"), type="normal")
  
res.joro.azul <- res
res.joro.azul

# Define significance threshold and filter results based on that
alpha = 0.01
res.joro.azul = res.joro.azul[order(res.joro.azul$padj, na.last=NA), ]
sigtab = res.joro.azul[(res.joro.azul$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(ps)[rownames(sigtab), ], "matrix"))

# Make a graph and export a result table
#negative values represent taxa enriched by level2, 
  #positive values represent taxa enriched by level1
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))


# co_pvf co_pvpf co_pfvf
# de_pvf de_pvpf de_pfvf
co_pvf <- 
  ggplot(sigtabgen, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  coord_flip() +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  theme(strip.text.x = element_text(size = 20, colour = "black"), legend.title=element_text(size=14), legend.text=element_text(size=14), axis.text=element_text(size=14)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) + geom_hline(yintercept=0, color = "grey", size=1) +
  theme(axis.text.y = element_text(face = "italic"))

#write.table(sigtabgen, "co_phyloseq_pvf.csv")
#write.table(sigtabgen, "co_phyloseq_pvpf.csv")
#write.table(sigtabgen, "co_phyloseq_pfvf.csv")

#write.table(sigtabgen, "de_phyloseq_pvf.csv")
#write.table(sigtabgen, "de_phyloseq_pvpf.csv")
#write.table(sigtabgen, "de_phyloseq_pfvf.csv")

# co_pvf co_pvpf co_pfvf
# de_pvf de_pvpf de_pfvf


figure <- ggarrange(coPLFA, dePLFA,co16S, de16S,
                    labels = c("A", "B", "C", "D"),
                    ncol = 3, nrow = 2,
                    common.legend = TRUE, legend = "bottom")
figure


