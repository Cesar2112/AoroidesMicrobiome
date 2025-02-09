######## R packages needed ########
# "phyloseq", "knitr", "xtable", "phangorn", "decipher", "microbiome", "microbiomeutilities", "eulerr", "ggplot2", 
# "vegan","ampvis2", "lme4", "lmerTest"

## Phyloseq:
# ASV_Table
ASV_16S1_Table <- read.delim("ASV_Table_16S1_tophyseq.txt", header = TRUE, sep = "\t", row.names=1) 
ASV_16S1_Table<-as.matrix(ASV_16S1_Table)
ASV = otu_table(ASV_16S1_Table, taxa_are_rows = TRUE)
# Taxonomy_Table
ASV_Tax_16S1 <- read.delim("Taxa_Table_161S1_tophyseq.txt", header = TRUE, sep = "\t", row.names = 1)
ASV_Tax_16S1<-as.matrix(ASV_Tax_16S1)
TAX = tax_table(ASV_Tax_16S1)
# Metadata
Mapfile=read.delim("MappingFileC.txt" , header=TRUE, sep="\t", row.names=NULL)
rownames(sampledata)=sampledata$SampleID
sampledata = sample_data(Mapfile)
# Phylogenetic tree
seqs <- getSequences(TAX)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
phy_tree = phy_tree(fitGTR$tree)

physeq_ori = phyloseq(ASV, TAX, sampledata, phy_tree)
"
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 446 taxa and 30 samples ]
sample_data() Sample Data:       [ 78 samples by 6 sample variables ]
tax_table()   Taxonomy Table:    [ 446 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 446 tips and 444 internal nodes ]"

# Prune taxa not present in any sample (if they exist)
physeq_ori_NoZeroTaxa <- prune_taxa(taxa_sums(physeq_ori) > 0, physeq_ori)
physeq16s1<-physeq_ori_NoZeroTaxa
"phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 365 taxa and 30 samples ]"

# Taxa filtering by prevalecence
prevdf_AO = apply(X = otu_table(physeq16s1),
                  MARGIN = ifelse(taxa_are_rows(physeq16s1), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
prevdf_AO = data.frame(Prevalence = prevdf_AO,
                       TotalAbundance = taxa_sums(physeq16s1),
                       tax_table(physeq16s1))
plyr::ddply(prevdf_AO, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                     sum(df1$Prevalence))}) -> dfprev_AO
kable(dfprev_AO)
write.table(dfprev_AO, "PhylumReadsPrev_161S1_AO_FN.tsv", sep = "\t", quote=F, col.names=NA)
filterPhyla = c("Bacteroidetes", "Campilobacterota", "Chlamydiae", "Deinococcus-Thermus", "Rhodothermaeota") 
(physeq16FIN = subset_taxa(physeq16s1, !Phylum %in% filterPhyla))
# physeq16FIN [ 351 taxa and 30 samples ]

# Change Taxa names par ASV#
taxa_names(Physeq_MS)
taxa_names(Physeq_MS) <- paste0("ASV", seq(ntaxa(physeq16FIN)))

# Overview of taxonomic assignement
TaxOverviwAO <-table(tax_table(physeq16FIN)[, 2])
write.table(TaxOverviwAO, "TaxOverview_AO_FN.tsv", sep = "\t", quote=F, col.names=NA)

# Number of features for each taxonomic level
PhylaTable_AO<-table(tax_table(physeq16FIN)[, "Phylum"], exclude = NULL)
write.table(PhylaTable_AO, "PhylaTable16s1AO_FN.tsv", sep = "\t", quote=F, col.names=NA)

# Relative abundance
physeq16FIN_AbRel = transform_sample_counts(physeq16FIN, function(x) x / sum(x) )

#***** Alpha Diversity*****
set.seed(12345)
#Physeq object used for aplpha diversity has 
plot_richness(physeq16s1Ori, title ="Agelas oroides Alpha Diversity", 
              measures=c("observed","Shannon", "PD", "Chao1")) + 
              geom_point(size=2)+ geom_boxplot()
Div_tabAO <- estimate_richness(physeq16s1Ori, measures=c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
alpha_div_AO <- estimate_pd(physeq16s1Ori)
Div_tabAO <- cbind(sample_data(physeq16s1Ori), alpha_div_AO, Div_tabAO) 
write.table(Div_tabAO, "Div_tab_AO_MS_FN.tsv", sep = "\t", quote=F, col.names=NA)

#Normal distribution test
shapiro.test(Div_tabAO$PD)
"data:  Div_tabAO$PD
W = 0.60427, p-value = 0.002162"
shapiro.test(Div_tabAO$Chao1)
"data:  Div_tabAO$Chao1
W = 0.9813, p-value = 0.49"    # Normal distributed
shapiro.test(Div_tabAO$Shannon)
"data:  Div_tabAO$Shannon
W = 0.81979, p-value = 0.0001538"
shapiro.test(Div_tabAO$Observed)
"data:  Div_tabAO$Observed     # Normal distributed
W = 0.96821, p-value = 0.4915"
shapiro.test(Div_tabAO$Simpson)
"data:  Div_tabAO$Simpson
W = 0.56354, p-value = 0.00000002695"
                                            
Linear fixed models analysis:
Random effects: Agelas oroides specimens 
Fixed effects: Sponge body parts "SpoTissue" and Enriched fractions "Fractions"
model = lmer(Chao1~SpoTissue + (1|spp), data=Div_tabAO) 
"Linear mixed model fit by REML ['lmerMod']
REML criterion at convergence: 181.1
Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-2.48059 -0.54063 -0.08892  0.78904  1.30637 
Random effects:
 Groups   Name        Variance Std.Dev.
 spp      (Intercept) 23.39    4.836   
 Residual             27.45    5.240   
Number of obs: 30, groups:  spp, 6
Fixed effects:
             Estimate Std. Error t value
(Intercept)   37.0000     2.4872  14.876
SpoTissueGnl   6.6667     2.6199   2.545
SpoTissueInt  -0.3333     2.1391  -0.156
Correlation of Fixed Effects:
            (Intr) SpTssG
SpoTissuGnl -0.351       
SpoTissuInt -0.430  0.408"
#P-value with lmerTest package
summary(lmer(Chao1 ~ SpoTissue + (1|spp), data = Div_tabAOb))
"Fixed effects:
		             Estimate Std. Error      df t value    Pr(>|t|)    
		(Intercept)   37.0000     2.4872  8.1087  14.876 0.000000359 ***
		SpoTissueGnl   6.6667     2.6199 22.0000   2.545      0.0185 *  
		SpoTissueInt  -0.3333     2.1391 22.0000  -0.156      0.8776    
		---
		Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

		Correlation of Fixed Effects:
		            (Intr) SpTssG
		SpoTissuGnl -0.351       
		SpoTissuInt -0.430  0.408"
					    
summary(lmer(Chao1 ~ Fraction + (1|spp), data = Div_tabAOb))
"Fixed effects:
                 Estimate Std. Error     df t value    Pr(>|t|)    
(Intercept)        35.750      2.475  7.959  14.447 0.000000542 ***
FractionCellular    2.167      2.090 22.000   1.037     0.31111    
FractionGeneral     7.917      2.560 22.000   3.093     0.00531 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) FrctnC
FractnClllr -0.422       
FractinGnrl -0.345  0.408"
anova(Model01,Model02)
Models:
Model01: Chao1 ~ SpoTissue + (1 | spp)
Model02: Chao1 ~ SpoTissue + Fraction + (1 | spp)
        npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
Model01    5 201.29 208.30 -95.645   191.29                     
Model02    6 202.14 210.55 -95.072   190.14 1.1461  1     0.2844                   

#***** Beta Diversity*****
# Body part example
Physeq16s1.metaAO <- meta(Div_tabAO)
# Create a list with the factors to compare
pares.spps <- combn(seq_along(SpoT_Levels), 2, simplify = TRUE, FUN = function(i)SpoT_Levels[i])
# Calcul divergences 
referenceBetaDiv3 <- apply(abundances(physeq16FIN), 1, median)
div_AOint <- divergence(subset_samples(physeq16FIN, SpoTissue == "Int"), referenceBetaDiv3, method = "bray")
data.frame(div_AOint) -> df_div_AOint
mutate(df_div_AOint, SpoTissue = "Int") -> df_div_AOint
colnames(df_div_AOint3) <- c("divergence", "SpoTissue")
# Graphic
plot <- ggviolin(title= "Agelas oroides microbiome divergences Physeq_MS", data = div.boxplot3,
                 x = "SpoTissue", y = "divergence", fill = "SpoTissue",
                 add = "boxplot", palette = c("#fdbf6f","#b2df8a","#a6cee3"))
plot + stat_compare_means(comparisons = pares.SpoT)
plot + stat_compare_means()#Kruskal Wallis p = 0.15

# Ordination methods
# Example Sponge body part
physeq.mds.unifrac <- ordinate(physeq16FIN, method = "MDS", distance = "unifrac")
physeq.mds.unifrac
evals <- physeq.mds.unifrac$values$Eigenvalues
pUnif <- plot_ordination(physeq16FIN, physeq.mds.unifrac, color = "SpoTissue", title = "Unifrac distance ordination") +
  labs(col = "A oriodes Sponge Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

physeq.mds.Weigtunifrac <- ordinate(physeq16FIN, method = "MDS", distance = "wunifrac")
evals <- physeq.mds.Weigtunifrac$values$Eigenvalues
pWeighUnif <- plot_ordination(physeq16FIN, physeq.mds.Weigtunifrac, color = "SpoTissue", title = "Weighted-Unifrac distance ordination") +
  labs(col = "A oriodes Sponge Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

physeq.mds.dpcoa <- ordinate(physeq16FIN, method = "MDS", distance = "dpcoa")
evals <- physeq.mds.dpcoa$values$Eigenvalues
pdpcoa_b <- plot_ordination(physeq16FIN, physeq.mds.dpcoa, color = "SpoTissue", title = "DPCoA ordination") +
  labs(col = "Sponge Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

#Permanova test with random models were conducted with Primer7

# Barplots
BPPhylaAgelas_MS <- plot_bar(physeq16FIN_AbRel, fill="Phylum") +
  scale_fill_manual(values = getPalette(colorCount)) +
  guides(fill=guide_legend(ncol=2))+ 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_blank())
BPPhylaAgelas_MS +  facet_wrap (~ SpoTissue, scales = "free")

# Venn diagrams
#create a list with factors
SpoTissueListVenn <- unique(as.character(meta(physeq16FIN_AbRel)$SpoTissue))
print(SpoTissueListVenn)#[1] "Ext" "Gnl" "Int"
list_core <- c() #  empty object to store information

for (n in SpoTissueListVenn){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(physeq16FIN_AbRel, SpoTissue == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  print(list_core)
}

---------------------------------------------------------------------------------------------
