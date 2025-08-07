.libPaths("C:/Users/00093406/Dropbox/Software/Rpackages")
library(ggplot2)
library(phyloseq)
library(plyr)
library(nlme)
library(multcomp)
library(ape)
library(picante)
library(tidyverse)
library(rcompanion)
library(multcompView)
library(lsmeans)
library(seqinr)
library(BIOM.utils)
library(DESeq2)
load()
load()
load()
setwd("C:/Users/00093406/Dropbox/Collabs/NESP streams/Nurdi data/microbiome/bioinformatics")

## export tax table to excel ##
ra<-transform_sample_counts(G, function(x){x / sum(x)})
otus <- t(ra@otu_table)
otus=as.data.frame(otus)
tax <- ra@tax_table
tax=as.data.frame(tax)
full_table=cbind(tax,otus)
write.csv(full_table,file = "before_filter_taxa_RA.csv")

#ADD CORRECTED MINERAL N DATA TO PHYLOSEQ OBJECT! ####
# map<-import_qiime_sample_data("frass_all_data.txt") #mapping file with added variables
# sam.new <- data.frame(map[,18:19]) #select variables to add
# #sam.new <- sam.new[-1, ] #remove blank row
# sam.new <- sample_data(sam.new) #turn into sample data
# head(sam.new) #check
# # Merge with original phyloseq object
# G1 <- merge_phyloseq(G, sam.new)
# rare2 <- merge_phyloseq(rare, sam.new)
# glom2 <- merge_phyloseq(glom, sam.new)
# head(sample_data(G1)) #check it worked
# head(sample_data(rare2))
# head(sample_data(glom2))
# G=G1 #re-name to save
# rare=rare2
# glom=glom2
# #re-save files
# save(G,file="frass.phyloseq")#filtered dataset. Not rarefied or glommed
# save(rare,file="frass_rare.phyloseq")#filtered dataset. rarefied
# save(glom,file="frass_glom.phyloseq")

taxa=readRDS("C:/Users/00093406/Dropbox/Collabs/NESP streams/Nurdi data/microbiome/bioinformatics/nurdi_tax.rds")
seqtab.nochim=readRDS("C:/Users/00093406/Dropbox/Collabs/NESP streams/Nurdi data/microbiome/bioinformatics/nurdi_seqtab.nochim.rds")#
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))
saveRDS(ps, file="C:/Users/21211201/Dropbox/PhD/Collabs/Sasha Jenkins/BSF/Frass/DADA2/BSF_CG6MJ.rds")
map<-import_qiime_sample_data("sample_names_nurdi.txt")
map<-sample_data(map)#mapping file
class(map) #should return "sample_data"
dim(map) #return number of samples and number of variables i.e. the data dimensions
G <- merge_phyloseq(ps,map)
G
G1=G

before_filter=as.data.frame(sample_sums(G))
write.csv(before_filter,file="reads_before_filter.csv")
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 40981 taxa and 233 samples ]
# sample_data() Sample Data:       [ 233 samples by 5 sample variables ]
# tax_table()   Taxonomy Table:    [ 40981 taxa by 6 taxonomic ranks ]

#remove class chloroplast 
G<-subset_taxa(G, !is.na(Order) & !Class %in% c("Chloroplast"))
G<-subset_taxa(G, !is.na(Class) & !Order %in% c("Chloroplast"))
G<-subset_taxa(G, !is.na(Family) & !Family %in% c("Mitochondria"))
G

#blank
blank<-subset_samples(G,owner=="blank")
blank10 = names(sort(taxa_sums(blank), TRUE)[1:10])
blank10= prune_taxa(blank10, blank)
blank10df = psmelt(blank10) 
write.csv(blank10df,file="blank10.csv")#write out the blanks for supplementary
G<-subset_samples(G,!owner=="blank") #remove=blank
sample_names(G)

#remove singletons and NAs
any(taxa_sums(G) < 1)
G = prune_taxa(taxa_sums(G) > 0, G)#gets rid of them
any(taxa_sums(G) < 1) #check it worked
any(is.na(otu_table(G))) #FALSE
any(otu_table(G) < 0) #FALSE
any(taxa_sums(G) < 2)
G = prune_taxa(taxa_sums(G) > 2, G)#gets rid of them
any(taxa_sums(G) < 2) #check it worked
any(is.na(otu_table(G))) #FALSE
any(otu_table(G) < 0) #FALSE
G

#Prevalance filtering
table(tax_table(G)[, "Phylum"], exclude = NULL)
prevdf = apply(X = otu_table(G),
               MARGIN = ifelse(taxa_are_rows(G), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(G),
                    tax_table(G))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#phyla that does not have a mean prevalence of 2 or greater (column 1) are suspect of contaminats. 
#Check by entering suspicious phyla names into the following code:
# Phylum        1     2
# 1        Acidobacteria 2.861374  2415
# 2       Actinobacteria 5.510045 13439
# 3      Armatimonadetes 2.084507   148
# 4        Bacteroidetes 2.856663 11898
# 5          Caldiserica 1.750000     7 #suspect - fine
# 6      Calditrichaeota 2.392857   201
# 7           Chlamydiae 1.440000    72 #suspect - fine
# 8          Chloroflexi 3.551756  3843
# 9        Cloacimonetes 2.628571    92
# 10       Cyanobacteria 4.279412   873
# 11     Deferribacteres 3.750000    15
# 12 Deinococcus-Thermus 1.607143    45 #suspect - fine
# 13        Dependentiae 1.471204   281 #suspect - fine
# 14       Elusimicrobia 1.000000     7 #suspect - remove
# 15   Entotheonellaeota 1.400000     7 #suspect - remove
# 16  Epsilonbacteraeota 5.345205  1951
# 17          Euglenozoa 1.000000     1 #suspect - remove
# 18       Euryarchaeota 2.506667   188
# 19       Fibrobacteres 3.646018   412
# 20          Firmicutes 3.051753  7312
# 21        Fusobacteria 2.333333    98
# 22    Gemmatimonadetes 2.736559  1018
# 23    Halanaerobiaeota 1.800000     9 #suspect - remove
# 24     Hydrogenedentes 1.225806    76 #suspect - fine
# 25  Kiritimatiellaeota 4.000000    12
# 26     Latescibacteria 2.052632   156
# 27       Lentisphaerae 1.000000     2 #suspect - remove
# 28      Modulibacteria 3.000000   309
# 29         Nitrospirae 2.538462    99
# 30     Patescibacteria 2.462963   133
# 31      Planctomycetes 1.577114   317 #suspect - fine
# 32      Proteobacteria 3.438402 43819
# 33        Rokubacteria 2.350000    47
# 34        Spirochaetes 2.145969   985
# 35       Synergistetes 1.000000     3 #suspect - remove
# 36         Tenericutes 2.468085   116
# 37      Thaumarchaeota 2.800000   126
# 38         Thermotogae 1.000000     1 #suspect - remove
# 39     Verrucomicrobia 2.206777  2540

subsample<-subset_taxa(G, Phylum == "Dependentiae")
plot_bar(subsample)

#have removed these sus phyla
G<-subset_taxa(G, !Phylum == "Elusimicrobia")
G<-subset_taxa(G, !Phylum == "Entotheonellaeota")
G<-subset_taxa(G, !Phylum=="Euglenozoa")
G<-subset_taxa(G, !Phylum=="Halanaerobiaeota")
G<-subset_taxa(G, !Phylum=="Lentisphaerae")
G<-subset_taxa(G, !Phylum=="Synergistetes")
G<-subset_taxa(G, !Phylum=="Thermotogae")

#check phyla were removed
table(tax_table(G)[, "Phylum"], exclude = NULL) 

#Subset phyla and plot with a line at 1% prevalence #check for suspect phyla and remove i.e. below the 1% line
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(G, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(G),color=Phylum))+  
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

prevalenceThreshold = 0.01 * nsamples(G)
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
G = prune_taxa(keepTaxa, G) 
G

#Export file of filtered taxa for supp material
finalFilter<-as.data.frame(sample_sums(G))
write.csv(finalFilter, file="reads_after_filter.csv")

jen=subset_samples(G,owner=="jen")
any(taxa_sums(jen) < 2)
jen = prune_taxa(taxa_sums(jen) > 2, jen)#gets rid of them
jen.glom<-tax_glom(jen,"Genus", NArm = TRUE)
save(jen,file="nurdi.phyloseq")
save(jen.glom,file="nurdi_glom.phyloseq")

maria=subset_samples(G,owner=="maria")
any(taxa_sums(maria) < 2)
maria = prune_taxa(taxa_sums(maria) > 2, maria)#gets rid of them
maria.glom<-tax_glom(maria,"Genus", NArm = TRUE)
save(maria,file="maria.phyloseq")
save(maria.glom,file="maria_glom.phyloseq")

luka=subset_samples(G,owner=="luka")
any(taxa_sums(luka) < 2)
luka = prune_taxa(taxa_sums(luka) > 2, luka)#gets rid of them
luka.glom<-tax_glom(luka,"Genus", NArm = TRUE)
save(luka,file="luka.phyloseq")
save(luka.glom,file="luka_glom.phyloseq")
luka.chloro=subset_samples(G1,owner=="luka")
save(luka.chloro,file="luka_before_filter.phyloseq")
#luka before filter RA
ra<-transform_sample_counts(luka.chloro, function(x){x / sum(x)})
otus <- t(ra@otu_table)
otus=as.data.frame(otus)
tax <- ra@tax_table
tax=as.data.frame(tax)
full_table=cbind(tax,otus)
write.csv(full_table,file = "luka_before_filter_taxa_RA.csv")

#Top 10 Genus
glomfam<-jen.glom
top<-glomfam
top10 = names(sort(taxa_sums(top), TRUE)[1:10]) 
top10= prune_taxa(top10, top)
top10df = psmelt(top10) #'psmelt' function turns a physolseq object into a data frame! handy!
#topdf<-psmelt(top)
write.csv(top10df,file="top10_genus.csv")#write out the blanks for supplementary
