## TOMATO & ARABIDOPSIS RHIZOSPHERE AMPLICON DATA ANALYSIS SCRIPT

# 16S data from rhizospheres of Solanum lycopersicum, Arabidopsis thaliana and unplanted soil

# Note: in the data, treatments are abbreviated by approximate soil moisture contents: "12" = well-watered = "WW" ; "9" = mild drought = "MDR" ; "6" = severe drought = "SDR"




setwd("C:/bryophytes_rhizosphere/")


# 0 - Loading libraries ----

library(phyloseq)
library(tidyr)
library(vegan)
library(dplyr)
library(ggplot2)
library(ggtext)



# 1 - Loading, filtering, rarefying, treebuilding

# 1.1 - Load 16S data ----
# MiSeq PE300 data was processed in qiime2-2022.2 with DADA2
# we load "seqtab_final.txt" (ASV table), "taxonomy.tsv" (taxonomy), "repset.fasta" (ASV sequences)

# define the path to the data. 
data_path_16S <- "./for_publication/code_for_github/angiosperms_experiment/data_angiosperms"


# loads the ASV table, immediately saving it as a phyloseq OTU table (lighter than a df)
raw_bac_otutab <- read.table(file = paste0(data_path_16S, "/feature-table.tsv"), 
                             header = TRUE,
                             row.names = 1, # first column has row names (ASV names)
                             check.names = FALSE,
                             quote = "")
raw_bac_otutab<-otu_table(object = raw_bac_otutab, taxa_are_rows = TRUE)
           


# taxonomy from qiime2-sklearn
raw_bac_taxtab<-read.table(file = paste0(data_path_16S, "/taxonomy.tsv"), 
                           header = TRUE,
                           sep = "\t",
                           row.names = 1, # first column has row names (ASV names)
                           check.names = FALSE) # prevents "X" to be added to column names, such as X49_16S,



# adjusts number and name of columns
raw_bac_taxtab<- separate(data = raw_bac_taxtab,
                          col = Taxon,
                          into =c("Kingdom", 
                                  "Phylum", 
                                  "Class", 
                                  "Order", 
                                  "Family",
                                  "Genus", 
                                  "Species"),
                          sep = ";")

# change from d__Bacteria to k__bacteria, matching fungi dataset
raw_bac_taxtab$Kingdom<-gsub("d__", "k__", raw_bac_taxtab$Kingdom)

# saves taxa as phyloseq object
raw_bac_taxtab<-tax_table(object = as.matrix(raw_bac_taxtab))



# loads the representative sequences table, immediately saving it as a phyloseq OTU table (lighter than a df)
raw_bac_refseq<-refseq(physeq = Biostrings::readDNAStringSet(filepath = paste0(data_path_16S, "/dna-sequences.fasta"), use.names = TRUE)) 
# taxa_names(raw_bac_refseq)<-gsub(" .*", "", taxa_names(raw_bac_refseq)) # drops taxonomy from ASV names

# loads the mapping file, immediately saving it as a phyloseq OTU table (lighter than a df)
raw_bac_metadata<-sample_data(object = read.table(file = paste0(data_path_16S, "/tomato_arabidopsis_mapping.txt"), 
                                                  header = TRUE,
                                                  sep = "\t",
                                                  row.names = 1, # first column has row names (ASV names)
                                                  check.names = FALSE)) # prevents "X" to be added to column names, such as X49_16S,

# build the main 16S phyloseq object by putting all these phyloseq-class objects (otu table, tax table, ref seq, metadata) into a single ps
raw_bac_ps<-merge_phyloseq(raw_bac_otutab,
                           raw_bac_taxtab,
                           raw_bac_refseq,
                           raw_bac_metadata)


# remove old objects to reduce memory use
rm(raw_bac_otutab,raw_bac_taxtab,raw_bac_refseq,raw_bac_metadata)

tax_table(raw_bac_ps)<-gsub(" ", "", tax_table(raw_bac_ps)) # drops a space character from taxa names. I don't know how that character ended up in there

# let's check the imported objects
otu_table(raw_bac_ps)[1:10,1:10]
sample_data(raw_bac_ps)[1:10,1:2]
tax_table(raw_bac_ps)[1:10,1:6]
refseq(raw_bac_ps)


# check dimensions of the different tables contained in the phyloseq object
dim(otu_table(raw_bac_ps))
dim(sample_data(raw_bac_ps))
dim(tax_table(raw_bac_ps))

# run garbage collection after creating large objects
gc()





# 1.2 - Filtering ----

# here we will remove sequences associated with the host, background prokaryotes, low abundance ASVs, poorly identified sequences, and other bad ASVs


dim(otu_table(raw_bac_ps)) # input data: 12516 ASVs

# we keep the raw phyloseq object and create a new one that will be filtered
filtered_bac_ps <- raw_bac_ps

# here we see that shorter sequences were still retained by the updated dada2 pipeline
hist(as.data.frame(refseq(filtered_bac_ps)@ranges)$width, breaks = 300)
summary(as.data.frame(refseq(filtered_bac_ps)@ranges)$width)

# here we get rid of sequences shorter than 380 bp
filtered_bac_ps<-prune_taxa(taxa = as.data.frame(refseq(filtered_bac_ps)@ranges)$width>379,
                       x = filtered_bac_ps)
dim(otu_table(filtered_bac_ps)) # 12486 ASVs

# histogram of sequence lengths after filtering <380 bp
hist(as.data.frame(refseq(filtered_bac_ps)@ranges)$width, breaks = 300)
summary(as.data.frame(refseq(filtered_bac_ps)@ranges)$width)

# keeps only ASVs identified as bacterial
filtered_bac_ps<-subset_taxa(filtered_bac_ps, Kingdom=="k__Bacteria")
dim(otu_table(filtered_bac_ps)) # 12425 ASVs


# removes taxa occuring less than 8 times in the dataset (changed compared to bryophytes experiment which had deeper sequencing, there the filter was for min. 20 counts)
otu_table(filtered_bac_ps) <- otu_table(filtered_bac_ps)[which (rowSums(otu_table(filtered_bac_ps)) > 7),]
dim(otu_table(filtered_bac_ps)) # 5883 ASVs

# filter ASV according samples (more than 0 reads, present in 3 samples)
filter <- phyloseq::genefilter_sample(filtered_bac_ps, filterfun_sample(function(x) x > 0), A = 3) 
filtered_bac_ps <- prune_taxa(filter, filtered_bac_ps)
rm(filter) # we no longer need this string
dim(otu_table(filtered_bac_ps)) # 3291 ASVs



# define plastid, mitochondria and host plant contamination ps objects
Mitochondria_ps<-subset_taxa(filtered_bac_ps, Family == "f__Mitochondria" | Family == "Mitochondria")
Plastid_ps<-subset_taxa(filtered_bac_ps, Order == "o__Chloroplast" | Order == "Chloroplast") 
host_plant_ps<-merge_phyloseq(Mitochondria_ps,Plastid_ps)

# quick histogram showing plant DNA contamination
hist(sample_sums(host_plant_ps)/sample_sums(filtered_bac_ps)*100, breaks = 100)

# define host plant 16S contamination as metadata
filtered_bac_ps@sam_data$Mitochondria_reads<-sample_sums(Mitochondria_ps)
filtered_bac_ps@sam_data$Plastid_reads-sample_sums(Plastid_ps)
filtered_bac_ps@sam_data$Host_DNA_n_reads<-sample_sums(host_plant_ps)
filtered_bac_ps@sam_data$Host_DNA_contamination_pct<-sample_sums(host_plant_ps)/sample_sums(filtered_bac_ps)*100



## remove plant host sequences (plastid and mitochondrial DNA)

# define function
remove_Chloroplast_Mitochondria<- function(physeq_object){
  
  #Removes Chloroplast
  #ATTENTION: if you just do >subset_taxa(physeq_object, Rank4!="o__Chloroplast") ; you will also remove NAs in the identification. thus you have to turn the ASV id into a factor before removing them, check https://hypocolypse.github.io/16s-data-prep.html
  
  # generate a df with Chloroplast ASVs
  CH1 <- subset_taxa(physeq_object, Order == "o__Chloroplast" | Order == "Chloroplast") # get all Chloroplasts...
  CH1 <-  as(tax_table(CH1), "matrix")
  CH1 <- row.names(CH1) # get ASV IDs...
  CH1df <- as.factor(CH1) # set IDs as factors
  goodTaxa <- setdiff(taxa_names(physeq_object), CH1df) #define taxa you should keep
  ps_no_chloro <- prune_taxa(goodTaxa, physeq_object) # your new physeq object is now chloroplast-free, but retains NA in identification
  
  
  #Removes Mitochondria
  #ATTENTION: if you just do >subset_taxa(physeq_chloro_mito, Rank5!="f__Mitochondria") ; you will also remove NAs in the identification. thus you have to turn the ASV id into a factor before removing them, check https://hypocolypse.github.io/16s-data-prep.html
  MT1 <- subset_taxa(physeq_object, Family == "f__Mitochondria" | Family == "Mitochondria")
  MT1 <-  as(tax_table(MT1), "matrix")
  MT1 <- row.names(MT1)
  MT1df <- as.factor(MT1)
  goodTaxa <- setdiff(taxa_names(physeq_object), MT1df)
  ps_no_chloro_mito <- prune_taxa(goodTaxa, ps_no_chloro)
  
  #excellent! let's save this a new phyloseq object
  physeq_clean<-ps_no_chloro_mito
  return(physeq_clean)
  
}

# use function to remove plastid and mitochondria ASVs
filtered_bac_ps<-remove_Chloroplast_Mitochondria(filtered_bac_ps)
dim(otu_table(filtered_bac_ps)) # 3253 ASVs



# add library sizes as part of metadata
sample_data(filtered_bac_ps)$library_size<-sample_sums(filtered_bac_ps)

# check library size distribution
hist(sample_data(filtered_bac_ps)$library_size)
summary(sample_data(filtered_bac_ps)$library_size)

# samples ordered by library size
sample_data(filtered_bac_ps)[order(sample_data(filtered_bac_ps)$library_size),]

# run garbage collection after creating large objects
gc()





## 1.3 - Rarefaction ----


# converts column of interest from numeric into factor
sample_data(filtered_bac_ps)$treatment<-as.factor(sample_data(filtered_bac_ps)$treatment)
sample_data(raw_bac_ps)$treatment<-as.factor(sample_data(raw_bac_ps)$treatment)

#draw rarefaction curve
filtered_bac_ps_df <- as.data.frame(otu_table(filtered_bac_ps))
rarecurve(t(filtered_bac_ps_df), 
          col = sample_data(filtered_bac_ps)$treatment, 
          label = FALSE, 
          step = 200,
          main="16S rarefaction curve, filtered dataset", ylab = "Number of ASVs", xlab = "Number of DNA reads",
          abline(v = 15741, col="red", lwd=3, lty=2))
# the red vertical dotted line indicates the depth at which we decide to rarefy (see below)


# samples ordered by library size
sample_data(filtered_bac_ps)[order(sample_data(filtered_bac_ps)$library_size),]

# rarefy to an even depth of 15741 reads (=lowest library size of all samples)
set.seed(100)
rarefied_bac_ps = rarefy_even_depth(filtered_bac_ps, 
                                             sample.size = 15741, 
                                             rngseed = FALSE, replace = TRUE, 
                                             trimOTUs = TRUE, verbose = TRUE)
# by rarefying at 15741 reads we lose 2 ASVs, but no samples

# check the data
rarefied_bac_ps
colSums(otu_table(rarefied_bac_ps))


# results of rarefaction - how many reads left
sum(otu_table(filtered_bac_ps)) # 1476746
sum(otu_table(rarefied_bac_ps)) # 771309 = 49 * 15741
# 2 ASVs lost, no samples lost. 48% of reads discarded


# now lets check rarefaction of the whole, unfiltered, dataset 
# because rarefaction of the filtered dataset is biased, less low-abundance ASVs

raw_bac_ps_df <- as.data.frame(otu_table(raw_bac_ps))
rarecurve(t(raw_bac_ps_df), 
          col = sample_data(raw_bac_ps)$treatment, 
          label = FALSE, 
          step = 200,
          main="16S rarefaction curve, unfiltered dataset", ylab = "Number of ASVs", xlab = "Number of DNA reads",
          abline(v = 15741, col="red", lwd=3, lty=2))




# we will continue working with the filtered and rarefied dataset: rarefied_bac_ps
rarefied_bac_ps  # it has 3251 ASVs and 49 samples



## 1.4 - Transform to relative abundance ----
rarefied_relative_bac_ps <- transform_sample_counts(rarefied_bac_ps, function(OTU) OTU/sum(OTU) )



## 1.5 - Save physeq objects ----

save(filtered_bac_ps, file = paste0(data_path_16S, "./data_angiosperms/filtered_bac_ps.RData"))
save(rarefied_bac_ps, file = paste0(data_path_16S, "./data_angiosperms/rarefied_bac_ps.RData"))
save(rarefied_relative_bac_ps, file = paste0(data_path_16S, "./data_angiosperms/rarefied_relative_bac_ps.RData"))




## 2 - Plots ----

# load physeq object (the one with relative abundances)
load(paste0(data_path_16S, "/rarefied_relative_bac_ps.RData"))


## 2.0 - Making melted dataframe with metadata and taxonomy ====

# use psmelt to combine the ASV table with metadata and taxonomy
bac_melted <- psmelt(rarefied_relative_bac_ps)




## 2.1 - Stacked bar plots - composition at phylum level ====

# define color-blind friendly palette
pal <- c("#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")



# group phyla
bac_phylum <- bac_melted %>% 
  group_by(Phylum, Sample, plant_species, treatment) %>%  # group by phylum and sample. also added plant_species and treatment as grouping factors in order to keep those columns
  summarise(Abundance = sum(Abundance))  # sum up the abundance of ASVs at phylum level

# remove "p__" in the phyla names
bac_phylum$Phylum <- gsub("p__", "", bac_phylum$Phylum)

# get the list of phyla by decreasing abundance, to order them
bac_phylum_order <- bac_melted %>% 
  group_by(Phylum) %>%
  summarise(Abundance = sum(Abundance))  # sum up the abundance for each phylum

# remove "p__" in the phyla names
bac_phylum_order$Phylum <- gsub("p__", "", bac_phylum_order$Phylum)
bac_phylum_order <- bac_phylum_order[order(bac_phylum_order$Abundance, decreasing=TRUE),]  # sort phyla by decreasing abundance


top_n_phyla <- 12

# replace phyla below the top N with "Other"
for (phylum_discard in bac_phylum_order$Phylum[(top_n_phyla+1):nrow(bac_phylum_order)]){
  bac_phylum$Phylum[bac_phylum$Phylum == phylum_discard] <- "Other"
}
# replace NA with "Other
bac_phylum$Phylum[is.na(bac_phylum$Phylum)] <- "Other"

# preparing the ordered list of phyla with only the top N and "Other"
bac_phyla_list <- bac_phylum_order$Phylum[1:top_n_phyla]
bac_phyla_list <- append(bac_phyla_list, "Other")
# ordering the phyla
bac_phylum$Phylum <- factor(bac_phylum$Phylum, levels = bac_phyla_list)

# ordering the treatments for the plot
bac_phylum$treatment <- factor(bac_phylum$treatment, levels = c("12", "9", "6"))

# make species labeller for the plot, to have entire names instead of abbreviations
# we use asterisks for the species so that ggtext::element_markdown() will then interprete them as italicization (see theme(strip.text)) in the plot
species_labeller <- as_labeller(c(`soil` = "unplanted soil", `tomato` = "*S. lycopersicum*",`arabidopsis` = "*A. thaliana*"))


# stacked bar plot by treatment for each plant species (as facets)
ggplot(bac_phylum, aes(fill=Phylum, y=Abundance, x=treatment)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(~factor(plant_species, levels=c("tomato", "arabidopsis", "soil")), labeller = species_labeller, scales="free_x", space="free_x") +
  scale_x_discrete(labels = c("6" = "severe drought","9" = "mild drought","12" = "well-watered"))+
  theme_bw()+
  xlab(NULL)+ # remove x axis label "treatment"
  ylab("relative abundance")+
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements ; angle x-axis labels
  scale_fill_manual(values=pal) +  # color-blind friendly palette
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = ggtext::element_markdown()) # facet label in italic for the 2 species, not unplanted soil. ggtext interpretes asterisks in the labeller as italicization

# export plot for Supp. Fig. 9
ggsave("R_plots/angiosperms_bac_phylum_barplot.png", width=2300, height=1300, unit="px")





## 3 - Differential abundance at phylum level ----

# we test differential abundance of phyla for the rhizopshere of tomato and Arabidopsis, but not for the unplanted soil control, as it only has 2 replicates.

alpha <- 0.05


## 3.1 - Wilcoxon test ====


# separate to each plant species
rarefied_relative_bac_ps_species_l <- metagMisc::phyloseq_sep_variable(rarefied_relative_bac_ps, variable = "plant_species")

# remove unplanted soil control
rarefied_relative_bac_ps_species_l <- within(rarefied_relative_bac_ps_species_l, rm(soil))

# agglomerate at phylum level
# note that the DESeq2 object will still have ASV IDs (of a random ASV corresponding to the phylum) as names, instead of phylum names.. we change this after getting the results tables from DESeq2
bac_phyla_rarefied_relative_ps_species_l <- lapply(rarefied_relative_bac_ps_species_l, function(x) tax_glom(x, taxrank = "Phylum"))




# lapply goes through the species
phyla_results <- lapply(bac_phyla_rarefied_relative_ps_species_l, function(x){
  
  
  # extract otu_table, tax_table and sample_data and make them dataframes
  phyla_table <- as.data.frame(otu_table(x))
  phyla_taxonomy <- as.data.frame(tax_table(x))
  sample_metadata <- as(sample_data(x), "data.frame") # for some reason it needs to be transformed with this function to not give an error with the merge() function later
  
  # replace ASV names with phyla names
  for (i in 1:nrow(phyla_table)){ # for loop goes through all the row numbers of our results table
    ASV_ID <- rownames(phyla_table)[i] # get the ASV ID of row i
    index <- which(rownames(phyla_taxonomy)==ASV_ID) # get the row index of the tax_table of filtered_bac_ps containing this ASV
    rownames(phyla_table)[i] <- phyla_taxonomy$Phylum[index] # put the phylum name instead of the ASV ID
  }
  
  # remove "p__" in phyla names
  rownames(phyla_table) <- gsub("p__", "", rownames(phyla_table))
  
  # move phyla names to new column and name that column "phylum"
  phyla_table <- cbind(rownames(phyla_table), phyla_table)
  colnames(phyla_table)[1] <- "phylum"
  
  # convert wide to long
  phyla_table <- phyla_table %>% pivot_longer(cols=colnames(phyla_table)[-1],
                                              names_to='sample_code',
                                              values_to='abundance')
  
  # subset columns plant_species and treatment from sample_metadata and move sample names to new column
  sample_metadata <- sample_metadata[, c("plant_species", "treatment")]
  sample_metadata$sample_code <- rownames(sample_metadata)
  
  # add sample metadata to phyla_table
  phyla_table <- merge(phyla_table, sample_metadata, by = "sample_code")
  
  
  # calculate average of each phylum, per treatment
  phyla_avg <- phyla_table %>% 
    group_by(phylum, treatment) %>%
    summarise(abundance = mean(abundance))
  
  # make new dataframe for the results, with a column with all phyla names
  results <- data.frame(phylum = unique(phyla_table$phylum))
  
  # for loop goes through all phyla and fills results table with LFC and p-values, both for MDR and SDR
  for (p in results$phylum){
    # LFC
    abundance_WW <- phyla_avg$abundance[which(phyla_avg$treatment == "12" & phyla_avg$phylum == p)]
    abundance_MDR <- phyla_avg$abundance[which(phyla_avg$treatment == "9" & phyla_avg$phylum == p)]
    results[which(results$phylum == p), "LFC_MDR_WW"] <- log2(abundance_MDR / abundance_WW)
    
    # p-value by Wilcoxon test
    values_WW <- subset(phyla_table, (treatment == "12" & phylum == p))$abundance
    values_MDR <- subset(phyla_table, (treatment == "9" & phylum == p))$abundance
    results[which(results$phylum == p), "pvalue_MDR_WW"] <- wilcox.test(values_MDR, values_WW)$p.value
    
    
    # calculate SDR vs. WW LFC and p-value
      # LFC
      abundance_SDR <- phyla_avg$abundance[which(phyla_avg$treatment == "6" & phyla_avg$phylum == p)]
      results[which(results$phylum == p), "LFC_SDR_WW"] <- log2(abundance_SDR / abundance_WW)
      
      # p-value by Wilcoxon test
      values_SDR <- subset(phyla_table, (treatment == "6" & phylum == p))$abundance
      results[which(results$phylum == p), "pvalue_SDR_WW"] <- wilcox.test(values_SDR, values_WW)$p.value
    
    
    
  }
  
  return(results)
  
})





# subset top 12 phyla from the bryophytes experiment. note: these may not be the top 12 phyla from the angiosperms experiment
top_12_phyla <- c("Actinobacteriota", "Proteobacteria", "Firmicutes", "Acidobacteriota", "Chloroflexi", "Planctomycetota", "Gemmatimonadota", "Verrucomicrobiota", "Bacteroidota", "Myxococcota", "Methylomirabilota", "Cyanobacteria" )
phyla_results <- lapply(phyla_results, function(x) x[which(x$phylum %in% top_12_phyla),])


## apply FDR-correction to the p-values
phyla_results <- lapply(phyla_results, function(x){
  x$pvalue_MDR_WW <- stats::p.adjust(p = x$pvalue_MDR_WW, method = "fdr")
  x$pvalue_SDR_WW <- stats::p.adjust(p = x$pvalue_SDR_WW, method = "fdr")
  return(x)
})




## 3.2 - Merged dataframe for heatmap ====


## MDR (mild drought)
res_MDR_WW_tomato <- data.frame(phylum = phyla_results$tomato$phylum, log2FoldChange = phyla_results$tomato$LFC_MDR_WW, padj = phyla_results$tomato$pvalue_MDR_WW)
res_MDR_WW_tomato$plant_species <- "S. lycopersicum"

res_MDR_WW_arabidopsis <- data.frame(phylum = phyla_results$arabidopsis$phylum, log2FoldChange = phyla_results$arabidopsis$LFC_MDR_WW, padj = phyla_results$arabidopsis$pvalue_MDR_WW)
res_MDR_WW_arabidopsis$plant_species <- "A. thaliana"

# merge the results tables
res_MDR_WW <- rbind(res_MDR_WW_tomato, res_MDR_WW_arabidopsis)

# add column "treatment"
res_MDR_WW$treatment <- "mild drought"



## SDR (severe drought)
res_SDR_WW_tomato <- data.frame(phylum = phyla_results$tomato$phylum, log2FoldChange = phyla_results$tomato$LFC_SDR_WW, padj = phyla_results$tomato$pvalue_SDR_WW)
res_SDR_WW_tomato$plant_species <- "S. lycopersicum"

res_SDR_WW_arabidopsis <- data.frame(phylum = phyla_results$arabidopsis$phylum, log2FoldChange = phyla_results$arabidopsis$LFC_SDR_WW, padj = phyla_results$arabidopsis$pvalue_SDR_WW)
res_SDR_WW_arabidopsis$plant_species <- "A. thaliana"

# merge all the results tables
res_SDR_WW <- rbind(res_SDR_WW_tomato, res_SDR_WW_arabidopsis)

# add column "treatment"
res_SDR_WW$treatment <- "severe drought"



## merge MDR and SDR
res_DR_WW_merged <- rbind(res_MDR_WW, res_SDR_WW)



# add column "significance" with asterisks for different significance levels (based on adjusted p-value)
# note that this is the adjusted p-value from testing one treatment (MDR or SDR) vs. WW in one plant species. so it is only adjusted for multiple testing with several phyla, not several treatments or species
res_DR_WW_merged$significance <- ""
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.05)] <- "*"
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.01)] <- "**"
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.001)] <- "***"




## 3.3 - Heatmap ====

# set LFC limits to -1 and +1 manually by replacing values under -1 with -1 and over 1 with 1
res_DR_WW_merged$log2FoldChange[which(res_DR_WW_merged$log2FoldChange < -1)] <- -1
res_DR_WW_merged$log2FoldChange[which(res_DR_WW_merged$log2FoldChange > 1)] <- 1

# ggplot heatmap
ggplot(res_DR_WW_merged, aes(x = treatment, y = phylum)) + 
  geom_tile(aes(fill = log2FoldChange))+ # LFC as colour
  geom_text(aes(label = significance), color = "black", size = 4, vjust = 0.8)+ # significance asterisks
  facet_grid(~factor(plant_species, levels=c("S. lycopersicum", "A. thaliana")), scales="free_x", space="free_x")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
  scale_y_discrete(limits = rev(top_12_phyla))+ # phyla are ordered from top to bottom by decreasing order of abundance in the whole dataset
  xlab(NULL)+
  ylab(NULL)+
  labs(fill= "Log2-fold Enrichment")+
  theme_bw()+
  theme(legend.position="bottom", legend.key.width = unit(1, 'cm')) +
  theme(text = element_text(size = 11), strip.text = element_text(size = 12)) +  # increase text size for all text elements, and facet title size even bigger
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic

# export plot for Fig. 4A
ggsave("R_plots/angiosperms_heatmap_bac_phylum_level.png", width=1500, height=1000, unit="px")





