#Set working directory 
#In the Pipelines subfolderuploaded all of the files to run this script and to generate all of the plots for the paper

setwd("~/Desktop/CHP_manuscript")
path <- "/Users/Desktop/CHP_manuscript"
list.files(path)

#Load required packages 
#Required: Need to load all of these packages for script to work
library(devtools)
library(phyloseq)
library(ggplot2)
library(biomformat)
library(ape)
library(scales)
library(grid)
library(gridExtra)
library(plyr)
library(vegan)
library(lfda)
library(decontam)
library(ggpubr)
library(tidyverse)
library(microbiome)
library(qiime2R)
library(tibble)
library(dplyr)
library(reshape2)
library(dendextend)
library(cowplot)
library(readxl)
library(survival)
library(survminer)
library(tidyr)
library(MASS)
library(rcompanion)
library(ggbiplot)
library(multcompView)
library(FSA) 
library(DescTools)
library(lattice)

############################################################################################################################################################################################################################
Cell counts  

df = read.csv("Cell-counts-barcharts.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)

mean_df = tapply(df$Macrophages, df$Diagnosis, mean) 
mean_df

SD_df = tapply(df$Macrophages, df$Diagnosis, sd) 
SD_df

N_df = tapply(df$Macrophages, df$Diagnosis, length) 
N_df

df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))

#Macrophages
p <- ggbarplot(df, x = "Diagnosis", y = "Macrophages", size = 0.5, title = "Macrophages", ylim = c(0,100), add = c("mean_se", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of total cells in BAL", x = "Diagnosis", color = "Diagnosis")  

#Lymphocytes
p <- ggbarplot(df, x = "Diagnosis", y = "Lymphocytes", size = 0.5, title = "Lymphocytes", ylim = c(0,100), add = c("mean_se", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) + labs(y="% of total cells in BAL", x = "Diagnosis", color = "Diagnosis")  

#Neutrophils 
p <- ggbarplot(df, x = "Diagnosis", y = "Neutrophils", size = 0.5, title = "Neutrophils", ylim = c(0,100), add = c("mean_se", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of total cells in BAL", x = "Diagnosis", color = "Diagnosis") 

##############################################################################################################################################################################################################################
df = read.csv("Survival_analysis_burden.csv", header = TRUE, na.strings=c("","NA"))
CHP_IPF <- as.data.frame(df)
CHP_IPF$Lymphocytes <- factor(CHP_IPF$Lymphocytes, levels = c("low", "high", "IPF"))
surv_object <- Surv(time = CHP_IPF$Survival, event = CHP_IPF$DEATH_FLAG)
fit1 <- survfit(surv_object ~ Lymphocytes, data = CHP_IPF)
summary(fit1)
fit.coxph <- coxph(surv_object ~ Lymphocytes + Age + Sex + Smoking_history + `FVC_%_predicted`, data = CHP_IPF)
summary(fit.coxph)
ggforest(fit.coxph, data = CHP_IPF)
ggsurvplot(fit1, data = CHP_IPF, pval = TRUE, risk.table = TRUE, break.x.by=365, xlab="Time (Days)", legend.labs=c("low", "high", "IPF"), palette = c("grey", "dodgerblue4", "cadetblue"), xlim = c(0, 1460))

####################################################################################################################################################################################################################################################################################
#16S qPCR bacterial burden Controls (n = 28) vs CHP (n = 110) vs IPF (n =45)

df = read.csv("qPCR-Ctrl-CHP-IPF(CHPmanuscript).csv", header = TRUE, row.names = 1, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=CopiesLog, color=Diagnosis)) + 
  geom_boxplot()
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) + labs(y="Log 16S rRNA gene copies/ml", x = "Diagnosis") + scale_y_continuous(trans='log10') + geom_jitter(position=position_jitter(0.1), size = 0.8) 

my_comparisons <- list( c("Controls", "CHP"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)    

my_comparisons <- list( c("Controls", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)   

my_comparisons <- list( c("CHP", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)   
####################################################################################################################################################################################################################################################################################
#Survival curve CHP vs IPF & Supplementary Table E1
df = read.csv("Survival_analysis_burden.csv", header = TRUE, na.strings=c("","NA"))
CHP_IPF <- as.data.frame(df)
CHP_IPF$disease <- factor(CHP_IPF$disease, levels = c("CHP", "IPF"))
surv_object <- Surv(time = CHP_IPF$Survival, event = CHP_IPF$DEATH_FLAG)
fit1 <- survfit(surv_object ~ disease, data = CHP_IPF)
summary(fit1)
fit.coxph <- coxph(surv_object ~ disease + Age + Sex + Smoking_history + `FVC_%_predicted` + honeycombing, data = CHP_IPF)
summary(fit.coxph)
ggforest(fit.coxph, data = CHP_IPF)
ggsurvplot(fit1, data = CHP_IPF, pval = TRUE, risk.table = TRUE, break.x.by=365, xlab="Time (Days)", legend.labs=c("CHP", "IPF"), palette = c("dodgerblue4", "cadetblue"), xlim = c(0, 1460))

###############################################################################################################################################################################################################################
#Decontamination script to generate all files for Microbiome plots 
#The following 3 files have been generated using the QIIME2-CHP-Script)

physeq <- qza_to_phyloseq("unfiltered-table.qza","unfiltered-rooted-tree.qza","taxonomy-unfiltered-final.qza")
mapping = read.csv("CHP_mapping_tapestation.csv", header = TRUE, row.names = 1)
sample_data(physeq)<-mapping

#Looking for contamination
ps<-physeq
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.freq <- isContaminant(ps, method="combined", conc="quant_reading",neg="is.neg", threshold=0.1)
table(contamdf.freq$contaminant)

#Just 30 out of the 3095 ASVs were classified as contaminants
head(which(contamdf.freq$contaminant))

# In this plot the dashed black line shows the model of a noncontaminant sequence feature for which frequency is expected to be independent of the input DNA concentration. The red line shows the model of a contaminant sequence feature, for which frequency is expected to be inversely proportional to input DNA concentration, as contaminating DNA will make up a larger fraction of the total DNA in samples with very little total DNA.
plot_frequency(ps, taxa_names(ps)[sample(which(contamdf.freq$contaminant),16)], conc="quant_reading") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

# Those all look like contaminants based on the associaiton with the red line!
# Next can look at Prevalence measures for these samples
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#FALSE  TRUE 
#3075  50

head(which(contamdf.prev$contaminant))

#Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#Samples seem to split pretty well though should have had more negative contols!
#Now I can do some manual curation to see what they are and how many reads they have
#First get interpreteable names
tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.freq$contaminant),]

#This leaves only the contaminants!
ps.contam <- prune_taxa(contamdf.freq$contaminant, ps)
summarize_phyloseq(ps.contam)

#Look at distribution
plot_bar(ps.contam)

#Need to identify the ASV above 1000 reads and see what they are
filter <- phyloseq::genefilter_sample(ps.contam, filterfun_sample(function(x) x >= 1000))
ps.contam.1k <- prune_taxa(filter, ps.contam)
otu_table<-as.data.frame(ps.contam@otu_table)

#check if tax is now a dataframe
tax_df <- as.data.frame(tax)
is.data.frame(tax_df)

contam<-left_join(rownames_to_column(otu_table), (rownames_to_column(tax_df)))
write.table(contam,file="contamination.txt", col.names=NA, row.names=T,sep="\t")

#If happy to remove all these then use following command
ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps)
ps.noncontam
summarize_phyloseq(ps.noncontam)

CHP<-subset_samples(ps.noncontam, Description=="BAL")
sum(taxa_sums((CHP)) == 0)

summarize_phyloseq(CHP)

#1] Min. number of reads = 7553 
#2] Max. number of reads = 1881019 
#3] Total number of reads = 10073056 
#4] Average number of reads = 55044.0218579235 
#5] Median number of reads = 45188 
#7] Sparsity = 0.964928449729424 
#6] Any OTU sum to 1 or less? YES 
#8] Number of singletons = 69 
#9] Percent of OTUs that are singletons 2.22940226171244 
#10] Number of sample variables are: 13 

CHP

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 3095 taxa and 183 samples ]
#sample_data() Sample Data:       [ 183 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 3095 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 3095 tips and 3061 internal nodes ]

summarize_phyloseq(CHP)

filterPhyla = c("p__Actinobacteria", "p__Bacteroidetes", "p__Firmicutes", "p__Fusobacteria", "p__Proteobacteria") # Filter entries with unidentified Phylum.
ps1 = subset_taxa(CHP, !(!Phylum %in% filterPhyla))
ps1

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2776 taxa and 183 samples ]
#sample_data() Sample Data:       [ 183 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 2776 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 2776 tips and 2745 internal nodes ]

summarize_phyloseq(ps1)

#1] Min. number of reads = 6899 
#2] Max. number of reads = 1870082 
#3] Total number of reads = 9911801 
#4] Average number of reads = 54162.8469945355 
#5] Median number of reads = 44863 
#7] Sparsity = 0.964049778743642 
#6] Any OTU sum to 1 or less? YES 
#8] Number of singletons = 57 
#9] Percent of OTUs that are singletons 2.05331412103746 
#10] Number of sample variables are: 13 

#Create final filtered object!
CHP_decontaminated<-ps1

#Remove contaminants manually (i.e. Sphingomonas) and after you have removed them create a final phyloseq object (i.e. CHP_decontaminated_final)
CHP_decontaminated@tax_table
badTaxa = c("bdb9e5f6c41c701ec4b3bc754f1fde53", "2b2ef93955a0e7fc9940b16759cf6c2a", "80a5cb9001fa1356ec097549337a0a64", "2974d9feee2aa60a8680500782e838c0","9f8b16f780589d38ccbb4a92b6cc258b", "ec62f5ab04d9c4ebc9aa3850f7e776c5", "066a26869b150e119af1ba42c9560827", "e78d9715d827eba81db3c87e840b3a16", "873532098886a3068e333839cff270e3", "6f05bbb7b495e6ffda5b2dacaf47c75c", "f0d3ad42c3754a9af2c18ab6f102bb28", "34b32933f7ccf504f3d591952229610c", "cc6669d31bc57687068868bc38e8e168", "654878c72a7532eac049ef9a09758f62", "74c3b8fe5db4ddf4eca09d816aedd8e2", "b1ebdc7dc82cf510768b062efab89fc9", "2e9f64ce9ec727ee1fde65c75aa71db7", "f65857c80ee4b6c9d502aaaed0ae786f", "78005887af64a050cbdfbfe834f60e1b", "14eeed0c81e495bdd16555cb4bc5d4c7", "871fcc3c84b0af8daad8418bebe9534c", "dd4e2672ac27cfe7a3d6b4fbf50c7ca5", "0b12835fa498ceaf800913bca8764435", "b9faf4c2c3f8a00a325acccd1a46ab93", "d524e847be1a6f4918baa05f5214f0a7", "09b3ffde2779e3c85da11372ef30cd1d", "98cd57089b10be0f057e407e11f06f71", "1006d42c9984f242af419f88f8f30b61", "147d6ecf30b0296147481cf4e4b58c28", "27c8a86b35998a0bc6cc4fc407a715df", "f3389b554fa53e0c16a261d364470eae", "2897f7c0c302e175164c6412f262d5b4")
goodTaxa <- setdiff(taxa_names(CHP_decontaminated), badTaxa)
CHP_decontaminated_final <- prune_taxa(goodTaxa, CHP_decontaminated)

#View feature table after removing manually contaminants as well as using the decontam package to identify and remove contaminants and negative control samples 
CHP_decontaminated_final

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2744 taxa and 183 samples ]
#sample_data() Sample Data:       [ 183 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 2744 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 2744 tips and 2717 internal nodes ]

summarize_phyloseq(CHP_decontaminated_final)

#1] Min. number of reads = 1760 
#2] Max. number of reads = 1637073 
#3] Total number of reads = 8581923 
#4] Average number of reads = 46895.7540983607 
#5] Median number of reads = 38651 
#7] Sparsity = 0.96797184916121 
#6] Any OTU sum to 1 or less? YES 
#8] Number of singletons = 57 
#9] Percent of OTUs that are singletons 2.07725947521866 
#10] Number of sample variables are: 13 

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(CHP_decontaminated_final)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet = prune_taxa(names(keepTaxa), CHP_decontaminated_final)

prunedSet

# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}

CHP_relative = transformSampleCounts(prunedSet, normalizeSample)

CHP_relative

#This is the FINAL filtered phyloseq object (183 samples and 495 taxa)

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 495 taxa and 183 samples ]
#sample_data() Sample Data:       [ 183 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 495 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 495 tips and 490 internal nodes ]

summarize_phyloseq(CHP_relative)

# Now will summarise results into 3 txt files (phylum, class, genus)
# Summarise at Phylum
CHP_relative_phylum<- aggregate_taxa(CHP_relative, 'Phylum')
otu_table<-as.data.frame(CHP_relative_phylum@otu_table)
tax<-as.data.frame(tax)
CHP_phylum<-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(CHP_phylum,file="CHP_phylum.txt", col.names=NA, row.names=T,sep="\t")

# Summarise at Class
CHP_relative_class<- aggregate_taxa(CHP_relative, 'Class')
otu_table<-as.data.frame(CHP_relative_class@otu_table)
tax<-as.data.frame(tax)
CHP_class<-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(CHP_class,file="CHP_class.txt", col.names=NA, row.names=T,sep="\t")

# Summarise at Genus
CHP_relative_genus<- aggregate_taxa(CHP_relative, 'Genus')
otu_table<-as.data.frame(CHP_relative_genus@otu_table)
tax<-as.data.frame(tax)
CHP_genus<-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(CHP_genus,file="CHP_genus.txt", col.names=NA, row.names=T,sep="\t")

###############################################################################################################################################################################################################################
#Taxonomy barcharts and the (A) phylum, (B) class, (C) genus levels 

CHP_phylum = read.csv("Phylum_average.csv", header = TRUE, na.strings=c("","NA"))
df <- melt(CHP_phylum, id.vars = "Diagnosis", measure.vars = c("Fusobacteria", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Firmicutes"), factorsAsStrings = TRUE, na.rm = TRUE)
colnames(df) <- c("Diagnosis", "Phylum", "Abundance")
p <- ggplot(df, aes(df$Diagnosis, df$Abundance)) + 
  geom_bar(aes(y = df$Abundance, x = df$Diagnosis, fill = df$Phylum),
           stat="identity", position = position_stack())
df$Diagnosis <- ordered(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))
levels(df$Diagnosis)
p + scale_fill_manual(values = c("darkgreen", "darkgoldenrod1", "lightskyblue", "darkblue", "firebrick")) + theme_classic() + labs(x = "Diagnosis", y = "Relative abundance (%)") + 
  guides(fill = guide_legend(title = "Phylum"))

CHP_class = read.csv("Class_average.csv", header = TRUE, na.strings=c("","NA"))
df <- melt(CHP_class, id.vars = "Diagnosis", measure.vars = c("Flavobacteriia", "Epsilonproteobacteria", "Coriobacteriia", "Alphaproteobacteria", "Fusobacteriia", "Betaproteobacteria", "Gammaproteobacteria", "Clostridia", "Bacilli", "Bacteroidia"), factorsAsStrings = TRUE, na.rm = TRUE)
colnames(df) <- c("Diagnosis", "Class", "Abundance")
p <- ggplot(df, aes(df$Diagnosis, df$Abundance)) + 
  geom_bar(aes(y = df$Abundance, x = df$Diagnosis, fill = df$Class),
           stat="identity", position = position_stack())
df$Diagnosis <- ordered(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))
levels(df$Diagnosis)
p + scale_fill_manual(values = c("darkseagreen", "forestgreen", "darkgreen", "gold", "sandybrown", "coral", "lightskyblue", "royalblue", "darkblue", "firebrick")) + theme_classic() + labs(x = "Diagnosis", y = "Relative abundance (%)") + 
  guides(fill = guide_legend(title = "Class"))

CHP_genus = read.csv("Genus_average.csv", header = TRUE, na.strings=c("","NA"))

df <- melt(CHP_genus, id.vars = "Diagnosis", measure.vars = c("Others", "Staphylococcus", "Actinobacillus",	"Fusobacterium",	"Rothia",	"Actinomyces","Haemophilus",	"Neisseria",	"Veillonella",	"Streptococcus",	"Prevotella"), factorsAsStrings = TRUE, na.rm = TRUE)
colnames(df) <- c("Diagnosis", "Genus", "Abundance")
p <- ggplot(df, aes(df$Diagnosis, df$Abundance)) + 
  geom_bar(aes(y = df$Abundance, x = df$Diagnosis, fill = df$Genus),
           stat="identity", position = position_stack())
df$Diagnosis <- ordered(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))
levels(df$Diagnosis)
p + scale_fill_manual(values = c("gainsboro", "darkseagreen", "forestgreen", "darkgreen", "gold", "sandybrown", "coral", "lightskyblue", "royalblue", "darkblue", "firebrick")) + theme_classic() + labs(x = "Diagnosis", y = "Relative abundance (%)") + 
  guides(fill = guide_legend(title = "Genus")) 

###################################################################################################################################################################################################################################################################################
#OPTIONAL: To Run Beta-diversity need to make sure that I filter out taxa that are not seen more than 2 times and at least in 1% of samples

#Remove taxa not seen more than 2 times in at least 1% of the samples. This protects against an OTU with small mean & trivially large C.V.
filteredset = filter_taxa(CHP_decontaminated_final, function(x) sum(x > 2) > (0.1*length(x)), TRUE)

filteredset

#otu_table()   OTU Table:         [ 179 taxa and 183 samples ]
#sample_data() Sample Data:       [ 183 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 179 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 179 tips and 178 internal nodes ]

summarize_phyloseq(filteredset)

#1] Min. number of reads = 609 
#2] Max. number of reads = 1511200 
#3] Total number of reads = 7746786 
#4] Average number of reads = 42332.1639344262 
#5] Median number of reads = 36133 
#7] Sparsity = 0.673260677107183 
#6] Any OTU sum to 1 or less? NO 
#8] Number of singletons = 0 
#9] Percent of OTUs that are singletons 0 
#10] Number of sample variables are: 13 

filteredset@tax_table

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(filteredset)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet_beta = prune_taxa(names(keepTaxa), filteredset)

prunedSet_beta

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 178 taxa and 183 samples ]
#sample_data() Sample Data:       [ 183 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 178 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 178 tips and 177 internal nodes ]

# To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}

CHP_relative_beta = transformSampleCounts(prunedSet_beta, normalizeSample)

CHP_relative_beta

#This is the FINAL filtered phyloseq object for beta-diversity analysis (PCA biplot) (183 samples and 178 taxa)

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 178 taxa and 183 samples ]
#sample_data() Sample Data:       [ 183 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 178 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 178 tips and 177 internal nodes ]

# Now will summarise results into 2 txt files (phylum & genus)
# Summarise at Phylum
CHP_relative_beta <- aggregate_taxa(CHP_relative_beta, 'Phylum')
otu_table<-as.data.frame(CHP_relative_beta@otu_table)
tax<-as.data.frame(tax)
CHP_relative_beta_phylum<-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(CHP_relative_beta_phylum,file="CHP-relative-beta-diversity-phylum.txt", col.names=NA, row.names=T,sep="\t")

# Summarise at Genus
CHP_relative_beta<- aggregate_taxa(CHP_relative_beta, 'Genus')
otu_table<-as.data.frame(CHP_relative_beta@otu_table)
tax<-as.data.frame(tax)
CHP_relative_beta_genus <-left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))
write.table(CHP_relative_beta_genus,file="CHP-relative-beta-diversity-genus.txt", col.names=NA, row.names=T,sep="\t")

#PCA
df = read.csv("PCA-phylum-beta.csv", header = TRUE, row.names = 1, na.strings=c("","NA")) #contains CHP_relative_beta dataset + additional filtering in clustvis (removed 7 IPF subjects that didn't make the cut)
df <- as.data.frame(df)

df$Diagnosis <- factor(df$Diagnosis, levels = c("CHP", "IPF"))

head(df)

CHP_IPF_PCA <- prcomp(df[1:2], scale. = TRUE)

ggbiplot(CHP_IPF_PCA, obs.scale = 0, var.scale = 0,
         groups = df$Diagnosis, ellipse = TRUE, circle = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') + theme_classic() + scale_color_manual(values=c("dodgerblue4", "cadetblue")) 

#Alpha-diversity: richness and evenness analysis
#Files used are in CHP_manuscript --> Pipelines --> CHPFinal-1 --> Unfiltered_198 --> Diversity-unfiltered-198-final
#Used unfiltered table for alpha diversity (ran on 198 samples)

#Chao1 comparison CHP vs IPF (Chao1 index)
df = read.csv("chao1", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("CHP", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=chao1, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("CHP", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + stat_compare_means(comparisons = my_comparisons) + scale_color_manual(values=c("dodgerblue4", "cadetblue")) + labs(y="Chao1 index") + geom_jitter(position=position_jitter(0.1), size = 1) 

#Shannon Index comparison CHP vs IPF (Shannon index)
df = read.csv("Shannon.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("CHP", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=shannon, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("CHP", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + stat_compare_means(comparisons = my_comparisons) + scale_color_manual(values=c("dodgerblue4", "cadetblue")) + labs(y="Shannon index") + geom_jitter(position=position_jitter(0.1), size = 1) 

#Observed OTUs comparison CHP vs IPF
df = read.csv("Observed.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("CHP", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=observed_otus, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("CHP", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + stat_compare_means(comparisons = my_comparisons) + scale_color_manual(values=c("dodgerblue4", "cadetblue")) + labs(y="Observed OTUs") + geom_jitter(position=position_jitter(0.1), size = 1) 

####################################################################################################################################################################################################################################################################################
#Negative control taxonomy plot and bacterial burden

Phylum_negative_controls <- read_csv("Phylum_negative_controls.csv")
df <- melt(Phylum_negative_controls, id.vars = "Diagnosis", measure.vars = c("Armatimonadetes", "Spirochaetes", "Acidobacteria", "Cyanobacteria", "Fusobacteria", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Firmicutes"), factorsAsStrings = TRUE, na.rm = TRUE)
colnames(df) <- c("Diagnosis", "Phylum", "Abundance")
p <- ggplot(df, aes(df$Diagnosis, df$Abundance)) + 
  geom_bar(aes(y = df$Abundance, x = df$Diagnosis, fill = df$Phylum),
           stat="identity", position = position_stack())
df$Diagnosis <- ordered(df$Diagnosis, levels = c("All BAL samples", "Neg.controls"))
levels(df$Diagnosis)
p + scale_fill_manual(values = c("lightsteelblue", "yellow", "darkorchid4", "mediumseagreen", "darkgreen", "darkgoldenrod1", "lightskyblue", "darkblue", "firebrick")) + theme_classic() + labs(x = "Diagnosis", y = "Relative abundance (%)") + 
  guides(fill = guide_legend(title = "Phylum"))

df = read.csv("qPCR-negative-controls.csv", header = TRUE, row.names = 1, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF", "Neg. Controls"))
p <- ggplot(df, aes(x=Diagnosis, y=Copies, color=Diagnosis)) + 
  geom_boxplot()
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue", "snow4")) + labs(y="Log 16S rRNA gene copies/ml", x = "Diagnosis") + scale_y_continuous(trans='log10') + geom_jitter(position=position_jitter(0.1), size = 0.8) 

####################################################################################################################################################################################################################################################################################
#Kaplan Meier survival curve and Cox proportional hazards model of disease progression in IPF

df = read.csv("Survival_analysis_burden.csv", header = TRUE, na.strings=c("","NA"))

survival<-df %>%
  filter(disease == "IPF") %>%
  rename(burden=Burden) %>%
  #filter(burden >= 10000)  %>%
  #filter(!is.na(fibrosis_extent)) %>%
  mutate(burdenlog=log(burden)) %>%
  mutate(tertile = ntile(burdenlog, 3))
#survival$composite<-as.numeric(survival$composite)
#survival$patmeasage<-as.numeric(survival$patmeasage)
#survival$age<-as.numeric(survival$age)
survival$DEATH_FLAG<-as.numeric(survival$DEATH_FLAG)
survival$tertile <- factor(survival$tertile, levels = c("1", "2", "3"))
#survival$smoking_history <- factor(survival$smoking_history, levels = c("Never", "Ex", "Current"))
#Univariant analysis
#making formulas
univ_formulas <- sapply(c("tertile","burdenlog","Sex","Age", "Smoking_history"),function(x)as.formula(paste('Surv(Survival,DEATH_FLAG)~',x)))
#making a list of models
univ_models <- lapply(univ_formulas, function(x){coxph(x,data=survival)})
#extract data (here I've gone for HR and confint and P Value)
univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)), pValue = summary(x)$coefficients[, 5]))})
#view univ_results for IPF (write down HR, CI, P cvalue)
univ_results
#Multivariant analysis
surv_object <- Surv(time = survival$Survival, event = survival$DEATH_FLAG)
fit.coxph <- coxph(surv_object ~ tertile + Sex + Age + Smoking_history + `FVC_%_predicted`, data = survival)
summary(fit.coxph)
ggforest(fit.coxph, data = survival)
#how does the model perform?
cox.zph(fit.coxph)
#plot(cox.zph((fit.coxph)))
#tertile plot
fit1 <- survfit(surv_object ~ tertile, data =survival)
ggsurvplot(fit1, data = survival, pval = TRUE, risk.table = FALSE, break.x.by=365, legend.labs=c("low","middle", "high"), palette = c("lightseagreen", "chocolate1", "darkorchid4"), xlim = c(0, 1460))

#Univariate model for CHP (surival analysis results)

df = read.csv("Survival_analysis_burden.csv", header = TRUE, na.strings=c("","NA"))

survival<-df %>%
  filter(disease == "CHP") %>%
  rename(burden=Burden) %>%
  #filter(burden >= 10000)  %>%
  #filter(!is.na(fibrosis_extent)) %>%
  mutate(burdenlog=log(burden)) %>%
  mutate(tertile = ntile(burdenlog, 3))
#survival$composite<-as.numeric(survival$composite)
#survival$patmeasage<-as.numeric(survival$patmeasage)
#survival$age<-as.numeric(survival$age)
survival$DEATH_FLAG<-as.numeric(survival$DEATH_FLAG)
survival$tertile <- factor(survival$tertile, levels = c("1", "2", "3"))
#survival$smoking_history <- factor(survival$smoking_history, levels = c("Never", "Ex", "Current"))
#Univariant analysis
#making formulas
univ_formulas <- sapply(c("tertile","burdenlog","Sex","Age", "Smoking_history"),function(x)as.formula(paste('Surv(Survival,DEATH_FLAG)~',x)))
#making a list of models
univ_models <- lapply(univ_formulas, function(x){coxph(x,data=survival)})
#extract data (here I've gone for HR and confint and P Value)
univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)), pValue = summary(x)$coefficients[, 5]))})
#view univ_results for IPF (write down HR, CI, P cvalue)
univ_results

#Result: No association between bacterial burden and survival in CHP subjects ((hazard ratio, 1.35; 95% confidence interval, 0.92-1.99)

####################################################################################################################################################################################################################################################################################
#Burden vs cell counts correlations

df = read.csv("Cell-counts-burden.csv", header = TRUE, na.strings=c("","NA"))

df <- as.data.frame(df)
plotNormalHistogram(df$Burden) # Not normally distributed log transform
plotNormalHistogram(df$Burdenlog) #Normally distributed
plotNormalHistogram(df$Macrophages.percentage) # normally distributed do not transform 
plotNormalHistogram(df$Lymphocytes.percentage) # not normally distributed (Box-Cox transformation)
plotNormalHistogram(df$Neutrophils.percentage) # not normally distributed (Box-Cox transformation)

Box = boxcox(df$Neutrophils.percentage ~ 1, lambda = seq(-6,6,0.1)) 
Cox = data.frame(Box$x, Box$y) 
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]
lambda = Cox2[1, "Box.x"]
T_df = (df ^ lambda - 1)/lambda 
T_df <- as.data.frame(T_df)

ggscatter(df, x = "Neutrophils.percentage", y = "Burdenlog", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Neutrophils (%)", ylab = "Bacterial burden (Log 16S rRNA gene copies/ml)") 

ggscatter(df, x = "Lymphocytes.percentage", y = "Burdenlog", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Lymphocytes (%)", ylab = "Bacterial burden (Log 16S rRNA gene copies/ml)") + scale_y_continuous(trans='log10')

ggscatter(df, x = "Macrophages.percentage", y = "Burdenlog", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Macrophages (%)", ylab = "Bacterial burden (Log 16S rRNA gene copies/ml)") + scale_y_continuous(trans='log10')

####################################################################################################################################################################################################################################################################################
#Using files generated after decontamination script (0.005% filtered and converted to percentage)

df = read.csv("Subject_Phylum.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)

df1 = read.csv("Subject_Class.csv", header = TRUE, na.strings=c("","NA"))
df1 <- as.data.frame(df1)

df2 = read.csv("Subject_Genus.csv", header = TRUE, na.strings=c("","NA"))
df2 <- as.data.frame(df2)

#Phylum
#Stats for each group (medians, descriptive statistics, histograms, kruskal-wallis test, dunn's test)
#Can account for multiple comparisons using different methods (i.e. bonferroni, FDR)

library(dplyr)
library(FSA)
library(lattice)

df = mutate(df, Diagnosis = factor(Diagnosis, levels=unique(Diagnosis)))
Summarize(Firmicutes ~ Diagnosis, data = df)
histogram(~ Firmicutes | Diagnosis, data=df,layout=c(1,3))


#Phylum plots:
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))

#Proteobacteria
p <- ggbarplot(df, x = "Diagnosis", y = "Proteobacteria", size = 0.5, title = "Proteobacteria", ylim = c(0,100), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Proteobacteria ~ Diagnosis, data = df)
MC_Proteobacteria = dunnTest(Proteobacteria ~ Diagnosis, data = df, method = "bonferroni")
MC_Proteobacteria

#Firmicutes
p <- ggbarplot(df, x = "Diagnosis", y = "Firmicutes", size = 0.5, title = "Firmicutes", ylim = c(0,100), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Firmicutes ~ Diagnosis, data = df)
MC_Firmicutes = dunnTest(Firmicutes ~ Diagnosis, data = df, method = "bonferroni", table = TRUE)
MC_Firmicutes

#Bacteroidetes
p <- ggbarplot(df, x = "Diagnosis", y = "Bacteroidetes", size = 0.5, title = "Bacteroidetes", ylim = c(0,100), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Bacteroidetes ~ Diagnosis, data = df)
MC_Bacteroidetes = dunnTest(Bacteroidetes ~ Diagnosis, data = df, method = "bonferroni", table = TRUE)
MC_Bacteroidetes

#Actinobacteria
p <- ggbarplot(df, x = "Diagnosis", y = "Actinobacteria", size = 0.5, title = "Actinobacteria", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Actinobacteria ~ Diagnosis, data = df)
MC_Actinobacteria = dunnTest(Actinobacteria ~ Diagnosis, data = df, method = "bonferroni", table = TRUE)
MC_Actinobacteria

#Fusobacteria
p <- ggbarplot(df, x = "Diagnosis", y = "Fusobacteria", size = 0.5, title = "Fusobacteria", ylim = c(0,20), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Fusobacteria ~ Diagnosis, data = df)
MC_Fusobacteria = dunnTest(Fusobacteria ~ Diagnosis, data = df, method = "bonferroni", table = TRUE)
MC_Fusobacteria

#Class

df1$Diagnosis <- factor(df1$Diagnosis, levels = c("Controls", "CHP", "IPF"))

#Bacteroidia
p <- ggbarplot(df1, x = "Diagnosis", y = "Bacteroidia", size = 0.5, title = "Bacteroidia", ylim = c(0,80), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Bacilli
p <- ggbarplot(df1, x = "Diagnosis", y = "Bacilli", size = 0.5, title = "Bacilli", ylim = c(0,80), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Clostridia
p <- ggbarplot(df1, x = "Diagnosis", y = "Clostridia", size = 0.5, title = "Clostridia", ylim = c(0,50), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Actinobacteria
p <- ggbarplot(df1, x = "Diagnosis", y = "Actinobacteria", size = 0.5, title = "Actinobacteria", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Gammaproteobacteria
p <- ggbarplot(df1, x = "Diagnosis", y = "Gammaproteobacteria", size = 0.5, title = "Gammaproteobacteria", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Betaproteobacteria
p <- ggbarplot(df1, x = "Diagnosis", y = "Betaproteobacteria", size = 0.5, title = "Betaproteobacteria", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Fusobacteriia
p <- ggbarplot(df1, x = "Diagnosis", y = "Fusobacteriia", size = 0.5, title = "Fusobacteriia", ylim = c(0,20), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Alphaproteobacteria
p <- ggbarplot(df1, x = "Diagnosis", y = "Alphaproteobacteria", size = 0.5, title = "Alphaproteobacteria", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Coriobacteriia
p <- ggbarplot(df1, x = "Diagnosis", y = "Coriobacteriia", size = 0.5, title = "Coriobacteriia", ylim = c(0,5), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis") 

#Epsilonproteobacteria
p <- ggbarplot(df1, x = "Diagnosis", y = "Epsilonproteobacteria", size = 0.5, title = "Epsilonproteobacteria", ylim = c(0,5), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Flavobacteriia 
p <- ggbarplot(df1, x = "Diagnosis", y = "Flavobacteriia", size = 0.5, title = "Flavobacteriia", ylim = c(0,5), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

#Genus

df2$Diagnosis <- factor(df2$Diagnosis, levels = c("Controls", "CHP", "IPF"))

#Prevotella
p <- ggbarplot(df2, x = "Diagnosis", y = "Prevotella", size = 0.5, title = "Prevotella", ylim = c(0,100), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Prevotella ~ Diagnosis, data = df2)
MC_Prevotella = dunnTest(Prevotella ~ Diagnosis, data = df2, method = "bonferroni", table = TRUE)
MC_Prevotella

#Streptococcus
p <- ggbarplot(df2, x = "Diagnosis", y = "Streptococcus", size = 0.5, title = "Streptococcus", ylim = c(0,80), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Streptococcus ~ Diagnosis, data = df2)
MC_Streptococcus = dunnTest(Streptococcus ~ Diagnosis, data = df2, method = "bonferroni")
MC_Streptococcus

#Veillonella
p <- ggbarplot(df2, x = "Diagnosis", y = "Veillonella", size = 0.5, title = "Veillonella", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Veillonella ~ Diagnosis, data = df2)
MC_Veillonella = dunnTest(Veillonella ~ Diagnosis, data = df2, method = "bonferroni")
MC_Veillonella

#Neisseria
p <- ggbarplot(df2, x = "Diagnosis", y = "Neisseria", size = 0.5, title = "Neisseria", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Neisseria ~ Diagnosis, data = df2)
MC_Neisseria = dunnTest(Neisseria ~ Diagnosis, data = df2, method = "bonferroni")
MC_Neisseria

#Haemophilus
p <- ggbarplot(df2, x = "Diagnosis", y = "Haemophilus", size = 0.5, title = "Haemophilus", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Haemophilus ~ Diagnosis, data = df2)
MC_Haemophilus = dunnTest(Haemophilus ~ Diagnosis, data = df2, method = "bonferroni")
MC_Haemophilus

#Actinomyces
p <- ggbarplot(df2, x = "Diagnosis", y = "Actinomyces", size = 0.5, title = "Actinomyces", ylim = c(0,40), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Actinomyces ~ Diagnosis, data = df2)
MC_Actinomyces = dunnTest(Actinomyces ~ Diagnosis, data = df2, method = "bonferroni")
MC_Actinomyces

#Rothia
p <- ggbarplot(df2, x = "Diagnosis", y = "Rothia", size = 0.5, title = "Rothia", ylim = c(0,20), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Rothia ~ Diagnosis, data = df2)
MC_Rothia = dunnTest(Rothia ~ Diagnosis, data = df2, method = "bonferroni")
MC_Rothia

#Fusobacterium
p <- ggbarplot(df2, x = "Diagnosis", y = "Fusobacterium", size = 0.5, title = "Fusobacterium", ylim = c(0,20), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Fusobacterium ~ Diagnosis, data = df2)
MC_Fusobacterium = dunnTest(Fusobacterium ~ Diagnosis, data = df2, method = "bonferroni")
MC_Fusobacterium

#Actinobacillus
p <- ggbarplot(df2, x = "Diagnosis", y = "Actinobacillus", size = 0.5, title = "Actinobacillus", ylim = c(0,5), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  

kruskal.test(Actinobacillus ~ Diagnosis, data = df2)
MC_Actinobacillus = dunnTest(Actinobacillus ~ Diagnosis, data = df2, method = "bonferroni")
MC_Actinobacillus

#Staphylococcus
p <- ggbarplot(df2, x = "Diagnosis", y = "Staphylococcus", size = 0.5, title = "Staphylococcus", ylim = c(0,5), add = c("mean_sd", "jitter"), color = "Diagnosis")
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) +  labs(y="% of 16S rRNA reads", x = "Diagnosis", color = "Diagnosis")  
kruskal.test(Staphylococcus ~ Diagnosis, data = df2)
MC_Staphylococcus = dunnTest(Staphylococcus ~ Diagnosis, data = df2, method = "bonferroni")
MC_Staphylococcus

##############################################################################################################################
#Correlations between airway microbiota and leukocyte types (CHP N = 110)

df = read.csv("Microbiota-Leukocytes.csv", header = TRUE, na.strings=c("","NA"))

#Macrophages

ggscatter(df, x = "Prevotella", y = "Macrophages", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Prevotella (%)", ylab = "Macrophages (%)") 

ggscatter(df, x = "Streptococcus", y = "Macrophages", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Streptococcus (%)", ylab = "Macrophages (%)") 

ggscatter(df, x = "Veillonella", y = "Macrophages", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Veillonella (%)", ylab = "Macrophages (%)") 

ggscatter(df, x = "Neisserialog", y = "Macrophages", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Neisseria (%)", ylab = "Macrophages (%)") 

ggscatter(df, x = "Haemophiluslog", y = "Macrophages", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Haemophilus (%)", ylab = "Macrophages (%)") 

ggscatter(df, x = "Staphylococcuslog", y = "Macrophages", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Staphylococcus (%)", ylab = "Macrophages (%)") 

#Lymphocytes 

ggscatter(df, x = "Prevotella", y = "Lymphocytes", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Prevotella (%)", ylab = "Lymphocytes (%)") 

ggscatter(df, x = "Streptococcus", y = "Lymphocytes", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Streptococcus (%)", ylab = "Lymphocytes (%)") 

ggscatter(df, x = "Veillonellalog", y = "Lymphocytes", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Veillonella (%)", ylab = "Lymphocytes(%)") 

ggscatter(df, x = "Neisseria", y = "Lymphocytes", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Neisseria (%)", ylab = "Lymphocytes (%)") 

ggscatter(df, x = "Haemophilus", y = "Lymphocytes", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Haemophilus (%)", ylab = "Lymphocytes (%)") 

ggscatter(df, x = "Staphylococcuslog", y = "Lymphocytes", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Staphylococcus (%)", ylab = "Lymphocytes (%)") 

#Neutrophils

ggscatter(df, x = "Prevotella", y = "Neutrophilslog", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Prevotella (%)", ylab = "Neutrophils (%)") 

res2 <-cor.test(df$Prevotella, df$Neutrophilslog,  method = "spearman")
re2

ggscatter(df, x = "Streptococcus", y = "Neutrophils", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Streptococcus (%)", ylab = "Neutrophils (%)") 

ggscatter(df, x = "Veillonella", y = "Neutrophils", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Veillonella (%)", ylab = "Neutrophils (%)") 

res2 <-cor.test(df$Veillonella, df$Neutrophils,  method = "spearman")
re2

ggscatter(df, x = "Neisserialog", y = "Neutrophils", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Neisseria (%)", ylab = "Neutrophils (%)") 

ggscatter(df, x = "Haemophilus", y = "Neutrophils", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Haemophilus (%)", ylab = "Neutrophils (%)") 

ggscatter(df, x = "Staphylococcuslog", y = "Neutrophils", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "Staphylococcus (%)", ylab = "Neutrophils (%)") 

##############################################################################################################################
#Correlations between Proteobacteria abundance and FEV1 measures 

df = read.csv("Proteobacteria-Lung-Function-Correlation.csv", header = TRUE, na.strings=c("","NA"))

ggscatter(df, x = "FEV1", y = "Proteobacteria", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = TRUE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "FEV1 % predicted", ylab = "Proteobacteria %") 


ggscatter(df, x = "FVC_1", y = "Proteobacteria", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = FALSE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "FVC % predicted", ylab = "Proteobacteria %") 

ggscatter(df, x = "DLCO", y = "Proteobacteria", color = "black", size = 1,
          conf.int = FALSE, ggtheme = theme_classic(), 
          cor.coef = FALSE, cor.method = "spearman", add.params = list(color = "red"), 
          xlab = "DLCO % predicted", ylab = "Proteobacteria %") 


######################################################################################################
#Negative bronchoscopy controls (n =6) vs BAL of matched patients (n = 6)

df <- read_csv("Bronch-Neg-Matched-Average.csv")
df <- melt(df, id.vars = "Diagnosis", measure.vars = c("Cyanobacteria", "Acidobacteria", "Actinobacteria", "Armatimonadetes", "Spirochaetes", "Elusimicrobia", "Bacteroidetes", "Proteobacteria", "Firmicutes", "Fusobacteria"), factorsAsStrings = TRUE, na.rm = TRUE)
colnames(df) <- c("Diagnosis", "Phylum", "Abundance")
p <- ggplot(df, aes(df$Diagnosis, df$Abundance)) + 
  geom_bar(aes(y = df$Abundance, x = df$Diagnosis, fill = df$Phylum),
           stat="identity", position = position_stack())
df$Diagnosis <- ordered(df$Diagnosis, levels = c("BAL", "Neg"))
levels(df$Diagnosis)
p + scale_fill_manual(values = c("lightsteelblue", "yellow", "darkorchid4", "mediumseagreen", "darkgreen", "darkgoldenrod1", "lightskyblue", "darkblue", "firebrick", "purple")) + theme_classic() + labs(x = "Diagnosis", y = "Relative abundance (%)") + 
  guides(fill = guide_legend(title = "Phylum"))

######################################################################################################
# Staph/Strep abundance disease progression

#CHP Staphylococcus
df = read.csv("Staph-Strep-CHP-Disease-Progression.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Staphgroup <- factor(df$Staphgroup, levels = c("Low", "High"))
surv_object <- Surv(time = df$Survival, event = df$DEATH_FLAG)
fit1 <- survfit(surv_object ~ Staphgroup, data = df)
summary(fit1)
fit.coxph <- coxph(surv_object ~ Staphgroup + Age + Sex + Smoking_history + `FVC_%_predicted`, data = df)
summary(fit.coxph)
ggforest(fit.coxph, data = df)
ggsurvplot(fit1, data = df, pval = TRUE, risk.table = TRUE, break.x.by=365, xlab="Time (Days)", legend.labs=c("Low", "High"), palette = c("darkblue", "red"), xlim = c(0, 1460))

#CHP Streptococcus
df = read.csv("Staph-Strep-CHP-Disease-Progression.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Strepgroup <- factor(df$Strepgroup, levels = c("Low", "High"))
surv_object <- Surv(time = df$Survival, event = df$DEATH_FLAG)
fit1 <- survfit(surv_object ~ Strepgroup, data = df)
summary(fit1)
fit.coxph <- coxph(surv_object ~ Strepgroup + Age + Sex + Smoking_history + `FVC_%_predicted`, data = df)
summary(fit.coxph)
ggforest(fit.coxph, data = df)
ggsurvplot(fit1, data = df, pval = TRUE, risk.table = TRUE, break.x.by=365, xlab="Time (Days)", legend.labs=c("Low", "High"), palette = c("dodgerblue4", "cadetblue"), xlim = c(0, 1460))

#IPF Staphylococcus
df = read.csv("Staph-Strep-IPF2-Disease-Progression.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Staphgroup <- factor(df$Staphgroup, levels = c("Low", "High"))
surv_object <- Surv(time = df$Survival, event = df$DEATH_FLAG)
fit1 <- survfit(surv_object ~ Staphgroup, data = df)
summary(fit1)
fit.coxph <- coxph(surv_object ~ Staphgroup + Age + Sex + Smoking_history + `FVC_%_predicted`, data = df)
summary(fit.coxph)
ggforest(fit.coxph, data = df)
ggsurvplot(fit1, data = df, pval = TRUE, risk.table = TRUE, break.x.by=365, xlab="Time (Days)", legend.labs=c("Low", "High"), palette = c("darkblue", "red"), xlim = c(0, 1460))

#IPF Streptococcus
df = read.csv("Staph-Strep-IPF2-Disease-Progression.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Strepgroup <- factor(df$Strepgroup, levels = c("Low", "High"))
surv_object <- Surv(time = df$Survival, event = df$DEATH_FLAG)
fit1 <- survfit(surv_object ~ Strepgroup, data = df)
summary(fit1)
fit.coxph <- coxph(surv_object ~ Strepgroup + Age + Sex + Smoking_history + `FVC_%_predicted`, data = df)
summary(fit.coxph)
ggforest(fit.coxph, data = df)
ggsurvplot(fit1, data = df, pval = TRUE, risk.table = TRUE, break.x.by=365, xlab="Time (Days)", legend.labs=c("Low", "High"), palette = c("darkblue", "red"), xlim = c(0, 1460))


########################################################################################################################################
#Alpha diversity (Ctrl vs CHP vs IPF - chao1 and shannon index)


#Chao1 comparison CHP vs IPF (Chao1 index)
df = read.csv("Diversity-Ctrl-CHP-IPF", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=chao1, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("Controls", "CHP", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) + labs(y="Chao1 index") + geom_jitter(position=position_jitter(0.1), size = 1) 

#Shannon Index comparison CHP vs IPF (Shannon index)
df = read.csv("Shannon.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=shannon, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("Controls", "CHP", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) + labs(y="Shannon index") + geom_jitter(position=position_jitter(0.1), size = 1) 

#Mann whitney U test
#By default, the function pairwise.wilcox.test() reports the pairwise adjusted (Holm) p-values.
pairwise.wilcox.test(df$shannon, sample_data(df)$Diagnosis)


#Chao1 Index comparison Ctrl vs IPF (Chao1 index)
df = read.csv("Shannon.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=chao1, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("Controls", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + stat_compare_means(comparisons = my_comparisons) + scale_color_manual(values=c("firebrick","cadetblue")) + labs(y="Shannon index") + geom_jitter(position=position_jitter(0.1), size = 1) 

#Shannon Index comparison Ctrl vs IPF (Shannon index)
df = read.csv("Shannon.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "IPF"))
p <- ggplot(df, aes(x=Diagnosis, y=shannon, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("Controls", "IPF"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + stat_compare_means(comparisons = my_comparisons) + scale_color_manual(values=c("firebrick","cadetblue")) + labs(y="Shannon index") + geom_jitter(position=position_jitter(0.1), size = 1) 

#Chao1 Index comparison Ctrl vs IPF (Chao1 index)
df = read.csv("Shannon.csv", header = TRUE, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP"))
p <- ggplot(df, aes(x=Diagnosis, y=chao1, color=Diagnosis)) + 
  geom_boxplot()
my_comparisons <- list( c("Controls", "CHP"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)                  
p + theme_classic() + stat_compare_means(comparisons = my_comparisons) + scale_color_manual(values=c("firebrick","dodgerblue4")) + labs(y="Shannon index") + geom_jitter(position=position_jitter(0.1), size = 1) 

#PCA (Controls vs CHP vs IPF)
#Statistically test PCA with PERMANOVA
library(ggplot2)
library(ggbiplot)

df = read.csv("PCA-beta-diversity.csv", header = TRUE, row.names = 1, na.strings=c("","NA")) #contains CHP_relative_beta dataset + additional filtering in clustvis (removed 7 IPF subjects that didn't make the cut)
df <- as.data.frame(df)

df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))

head(df)

CHP_IPF_PCA <- prcomp(df[,1:2], scale. = TRUE)

summary(CHP_IPF_PCA)

plot(CHP_IPF_PCA, type = 'l') #can plot the PC's accoding to variance using a scre plot plots the variances sq of StDev vs each PC #PCA1 and 2 show most variances

#Variance is square of standard deviation
#interpret PCA --> use function called biplot

biplot(CHP_IPF_PCA, scale = 0)

#each number is observation from dataframe

#need to attach the relevent PCA analysis to the dataframe
#arrows --> represent eigen vectors for each variable 

str(CHP_IPF_PCA) # looking at the structure of the PCA analysis x variable is the PCA analysis whth each PCA in a list 1-4

str(CHP_IPF_PCA$x)

#want to bind PCA 1 and 2 to our dataframe in r this can be done within r

CHP_IPF_PCA2 <- cbind(df, CHP_IPF_PCA$x[,1:2]) # sticks them together

head(CHP_IPF_PCA2)

#plot with ggplot2

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

ggplot(CHP_IPF_PCA2, aes(PC1, PC2))+
  geom_point()+theme

#Perform k-means clustering on PCA results only 
#create dataframe with just PCA results
#K means clustering - works by identifyign K-centres - number of clusters within your data - looks at distance betwen K-mean and each sample

kmPCA <-- data.frame(CHP_IPF_PCA$x[,1:2])
plot(kmPCA)
fitK <- kmeans(kmPCA, 3)  #create arg using kmeans function (dataset, number of clusters)
fitK
k <- list()
for(i in 1:6){
  k[[i]]<- kmeans(kmPCA, i)
}
betweenss_totss <- list()

for(i in 1:6){
  betweenss_totss[[i]] <- k[[i]]$betweenss/k[[i]]$totss
}
k
plot(1:6, betweenss_totss, type = "b", ylab = "Between SS/ Total SS", xlab = "Clusters (k)") #where you see a shoulder in your data is the optimial k means to select for your data # in my data-set clear that K=3 is correct cluster

#Bind k-membership to CHP_IPF_PCA2 --> create a new data frame

str(fitK)

CHP_IPF_PCA3 <- cbind(CHP_IPF_PCA2, fitK$cluster)

#convert cluster data into factor data
str(CHP_IPF_PCA3)

CLUSTER <- as.factor(CHP_IPF_PCA3$`fitK$cluster`)

CHP_IPF_PCA4 <- cbind(CHP_IPF_PCA3, CLUSTER)

str(CHP_IPF_PCA4)

#plot with ggfortify --> k-means analysis

autoplot(prcomp(CHP_IPF_PCA4[,1:6], scale = TRUE),data=CHP_IPF_PCA4,loadings=TRUE, loadings.label= TRUE, loadings.colour= 'blue',loadings.label.size = 3 , loadings.label.colour = 'black', colour= 'CLUSTER', frame = TRUE) + theme

#plot with ggfortify --> ungrouped analysis

autoplot(prcomp(CHP_IPF_PCA4[,1:6], scale = TRUE),data=CHP_IPF_PCA4,loadings=FALSE, loadings.label= FALSE, loadings.colour= 'blue',loadings.label.size = 3 , loadings.label.colour = 'black', colour = "Disease_Group", frame = TRUE) + theme

#correlations between var and PC ---- # allows you to interpret relationships between your original variables and the PCA analysis

PCAVARCorr <- data.frame(cor(CHP_IPF_PCA4[,1:6], CHP_IPF_PCA4[,1:6]))

##Final Code for PCA
df = read.csv("PCA-beta-diversity.csv", header = TRUE, row.names = 1, na.strings=c("","NA")) #contains CHP_relative_beta dataset + additional filtering in clustvis (removed 7 IPF subjects that didn't make the cut)
df <- as.data.frame(df)

df$Diagnosis <- factor(df$Diagnosis, levels = c("Controls", "CHP", "IPF"))

head(df)

CHP_IPF_PCA <- prcomp(df[,1:2], scale. = TRUE)

ggbiplot(CHP_IPF_PCA, obs.scale = 0, var.scale = 0,
         groups = df$Diagnosis, ellipse = TRUE, circle = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') + theme_classic() + scale_color_manual(values=c("firebrick", "dodgerblue4", "cadetblue")) 

#Statistically test for differences using PERMANOVA adjusted for confounding variables
adonis(wunifrac_dist)

#PERMANOVA FOR BIPLOT ABOVE
#Compute a distance matrix 
#Excel file with only microbiome data as matrix (x)

m <- dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
as.dist(m, diag = FALSE, upper = FALSE)

print(x, diag = NULL, upper = NULL,
      digits = getOption("digits"), justify = "none",
      right = TRUE)

as.matrix(x)

#Create a separe spreasdheet with confounding variables in a dataframe (call it df2)

#Statistically test for differences using PERMANOVA adjusted for confounding variables

adonis2(x ~ Diagnosis*Age*Gender*Smoking*FVC*DLCO, data = df2)

#PCA (FVC thresholds CHP & IPF)

df = read.csv("FVC-PCA.csv", header = TRUE, row.names = 1, na.strings=c("","NA")) #contains CHP_relative_beta dataset + additional filtering in clustvis (removed 7 IPF subjects that didn't make the cut)
df <- as.data.frame(df)

df$FVCgroup <- factor(df$FVCgroup, levels = c("low", "middle", "high"))

head(df)

FVC_PCA <- prcomp(df[1:2], scale. = TRUE)

ggbiplot(FVC_PCA, obs.scale = 0, var.scale = 0, var.axes = FALSE,
         groups = df$FVCgroup, ellipse = FALSE, circle = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') + theme_classic() + scale_color_manual(values=c("red", "gainsboro", "darkblue")) 

#Bacterial burden (FVC thresholds CHP & IPF)

df = read.csv("FVC-burden.csv", header = TRUE, row.names = 1, na.strings=c("","NA"))
df <- as.data.frame(df)
df$FVCgroup <- factor(df$FVCgroup, levels = c("low", "middle", "high"))
p <- ggplot(df, aes(x=FVCgroup, y=Burdenlog, color=FVCgroup)) + 
  geom_boxplot()
my_comparisons <- list( c("low", "middle"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)  
p + theme_classic() + stat_compare_means() +scale_color_manual(values=c("red", "gainsboro", "darkblue")) + labs(y="Log 16S rRNA gene copies/ml", x = "FVC (% predicted)") + scale_y_continuous(trans='log10') + geom_jitter(position=position_jitter(0.1), size = 0.8) 

#low vs middle:0.24
#low vs high: 0.63
#middle vs high: 0.049
#overall Kruskal-Wallis p value: 0.11

#Exposure vs No-exposure comparing burden in CHP subjects 
df = read.csv("Exposure-burden-CHP.csv", header = TRUE, row.names = 1, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Exposure <- factor(df$Exposure, levels = c("0", "1"))
p <- ggplot(df, aes(x=Exposure, y=Burden, color=Exposure)) + 
  geom_boxplot()
my_comparisons <- list( c("0", "1"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)  
p + theme_classic() +scale_color_manual(values=c("darkblue", "red")) + labs(y="Log 16S rRNA gene copies/ml", x = "Exposure") + scale_y_continuous(trans='log10') + geom_jitter(position=position_jitter(0.1), size = 0.8) 

df = read.csv("Exposure-burden-CHP.csv", header = TRUE, row.names = 1, na.strings=c("","NA"))
df <- as.data.frame(df)
df$Exposuregroup <- factor(df$Exposuregroup, levels = c("Bird", "Down", "Mould", "Other", "None"))
p <- ggplot(df, aes(x=Exposuregroup, y=Burden, color=Exposuregroup)) + 
  geom_boxplot()
my_comparisons <- list( c("Bird", "Down"))
p + stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)  
p + theme_classic() +scale_color_manual(values=c("indianred3", "lightblue", "darkolivegreen4", "orange", "darkblue")) + labs(y="Log 16S rRNA gene copies/ml", x = "Exposure") + scale_y_continuous(trans='log10') + geom_jitter(position=position_jitter(0.1), size = 0.8) 

#Exposure vs No-exposure PCA comparing phyla in CHP subjects (groups)
df = read.csv("Exposure-PCA-CHP.csv", header = TRUE, row.names = 1, na.strings=c("","NA"))

df <- as.data.frame(df)

df$Exposuregroup <- factor(df$Exposuregroup, levels = c("Bird", "Down", "Mould", "Other", "None"))

head(df)

Exposure_PCA <- prcomp(df[1:2], scale. = TRUE)

ggbiplot(Exposure_PCA, obs.scale = 0, var.scale = 0, var.axes = FALSE,
         groups = df$Exposuregroup, ellipse = FALSE, circle = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') + theme_classic() + scale_color_manual(values=c("indianred3", "lightblue", "darkolivegreen4", "orange", "darkblue")) 

#Exposure vs No-exposure PCA comparing phyla in CHP subjects 
df = read.csv("Exposure-PCA-CHP.csv", header = TRUE, row.names = 1, na.strings=c("","NA"))

df <- as.data.frame(df)

df$Exposure <- factor(df$Exposure, levels = c("0", "1"))

head(df)

Exposure_PCA <- prcomp(df[1:2], scale. = TRUE)

ggbiplot(Exposure_PCA, obs.scale = 0, var.scale = 0, var.axes = FALSE,
         groups = df$Exposure, ellipse = FALSE, circle = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') + theme_classic() + scale_color_manual(values=c("darkblue", "red")) 


