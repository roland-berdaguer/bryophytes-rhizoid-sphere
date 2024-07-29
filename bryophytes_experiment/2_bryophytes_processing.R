

# 0 - Loading packages ----


library(phyloseq)
library(metagenomeSeq)



# 1 - Loading, filtering, rarefying, treebuilding ----

# 1.1 - Load microbiome data ----
# our NovaSeq data was processed in R with DADA2/Ernakovich lab pipeline, and then in Qiime for assigning taxonomy
# we load, for both 16S and ITS: "seqtab_final.txt" (ASV table), "taxonomy.tsv" (taxonomy), "repset.fasta" (ASV sequences)

# define the paths to the data. 
data_path_16S <- "./data_bryophytes/16S"
data_path_ITS <- "./data_bryophytes/ITS"

## 1.1.1 - Loading 16S data ====

# loads the ASV table, immediately saving it as a phyloseq OTU table (lighter than a df)
raw_bac_otutab <- read.table(file = paste0(data_path_16S, "/seqtab_final.txt"), 
                             header = TRUE,
                             row.names = 1, # first column has row names (ASV names)
                             check.names = FALSE)
colnames(raw_bac_otutab) <- gsub("B-", "", colnames(raw_bac_otutab))  # remove "B-" in sample names
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
raw_bac_refseq<-refseq(physeq = Biostrings::readDNAStringSet(filepath = paste0(data_path_16S, "/repset.fasta"), use.names = TRUE)) 
taxa_names(raw_bac_refseq)<-gsub(" .*", "", taxa_names(raw_bac_refseq)) # drops taxonomy from ASV names

# loads the mapping file, immediately saving it as a phyloseq OTU table (lighter than a df)
raw_bac_metadata<-sample_data(object = read.table(file = "./data_bryophytes/bryophytes_mapping.txt", 
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



## 1.1.2 - Loading ITS data ====
# This will load all the essential data (OTU frequencies,  metadata, taxonomy) for a single phyloseq object

# loads the OTU table, immediatetly saving it as a phyloseq OTU table (lighter than a df)
raw_fun_otutab <- read.table(file = paste0(data_path_ITS, "/seqtab_final.txt"),
                             header = TRUE,
                             row.names = 1, # first column has row names (ASV names)
                             check.names = FALSE)
colnames(raw_fun_otutab) <- gsub("B-", "", colnames(raw_fun_otutab))  # remove "B-" in sample names
raw_fun_otutab<-otu_table(object = raw_fun_otutab, taxa_are_rows = TRUE)



# taxonomy from qiime2-sklearn
raw_fun_taxtab<-read.table(file = paste0(data_path_ITS, "/taxonomy.tsv"),
                           header = TRUE,
                           sep = "\t",
                           row.names = 1, # first column has row names (ASV names)
                           check.names = FALSE) # prevents "X" to be added to column names, such as X49_16S,

# adjusts number and name of columns
raw_fun_taxtab<- separate(data = raw_fun_taxtab,
                          col = Taxon,
                          into =c("Kingdom",
                                  "Phylum",
                                  "Class",
                                  "Order",
                                  "Family",
                                  "Genus",
                                  "Species"),
                          sep = ";")

# saves taxa as phyloseq object
raw_fun_taxtab<-tax_table(object = as.matrix(raw_fun_taxtab))




# loads the representative sequences table, immediatetly saving it as a phyloseq OTU table (lighter than a df)
raw_fun_refseq<-refseq(physeq = Biostrings::readDNAStringSet(filepath = paste0(data_path_ITS, "/repset.fasta"), use.names = TRUE))
taxa_names(raw_fun_refseq)<-gsub(" .*", "", taxa_names(raw_fun_refseq)) # drops taxonomy from ASV names

# loads the mapping tile, immediatetly saving it as a phyloseq OTU table (lighter than a df)
raw_fun_metadata<-sample_data(object = read.table(file = "./data_bryophytes/bryophytes_mapping.txt",
                                                  header = TRUE,
                                                  sep = "\t",
                                                  row.names = 1, # first column has row names (ASV names)
                                                  check.names = FALSE)) # prevents "X" to be added to column names, such as X49_ITS,

# build the main ITS phyloseq object by putting all these phyloseq-class objects (otu table, tax table, ref seq, metadata) into a single ps
raw_fun_ps<-merge_phyloseq(raw_fun_otutab,
                           raw_fun_taxtab,
                           raw_fun_refseq,
                           raw_fun_metadata)

# remove old objects to reduce memory use
rm(raw_fun_otutab,raw_fun_taxtab,raw_fun_refseq,raw_fun_metadata)

# change names from ASV to fASV, so they can be distinguished from bacterial ASVs
taxa_names(raw_fun_ps)<-paste("f", taxa_names(raw_fun_ps), sep = "")

# let's check the imported objects. Often errors will arise from typos when filling up the data sheets
otu_table(raw_fun_ps)[1:10,1:10]
sample_data(raw_fun_ps)[1:10,1:2]
tax_table(raw_fun_ps)[1:10,1:7]
refseq(raw_fun_ps)

# check dimensions of the different tables contained in the phyloseq object
dim(otu_table(raw_fun_ps)) # 17165 ASVs, 44 samples
dim(sample_data(raw_fun_ps))
dim(tax_table(raw_fun_ps))

# run garbage collection after creating large objects
gc()





# 1.2 - 16S filtering ----

# here we will remove sequences associated with the host, background prokaryotes, low abundance ASVs, poorly identified sequences, and other bad ASVs


## 1.2.1 - 16S "pre-filtering": remove non-bacterial ASVs ====

# remove ASVs shorter than 380 bp, ASVs that are not taxonomically assigned as bacteria, ASVs assigned as plastids/mitochondria

dim(otu_table(raw_bac_ps)) # input data: 28580 ASVs

# we keep the raw phyloseq object and create a new one that will be filtered
prefiltered_bac_ps <- raw_bac_ps

# here we see that shorter sequences were still retained by the updated dada2 pipeline
hist(as.data.frame(refseq(prefiltered_bac_ps)@ranges)$width, breaks = 300)
summary(as.data.frame(refseq(prefiltered_bac_ps)@ranges)$width)

# here we get rid of sequences shorter than 380 bp
prefiltered_bac_ps<-prune_taxa(taxa = as.data.frame(refseq(prefiltered_bac_ps)@ranges)$width>379,
                               x = prefiltered_bac_ps)
dim(otu_table(prefiltered_bac_ps)) # 28201 ASVs

# histogram of sequence lengths after filtering <380 bp
hist(as.data.frame(refseq(prefiltered_bac_ps)@ranges)$width, breaks = 300)
summary(as.data.frame(refseq(prefiltered_bac_ps)@ranges)$width)

# keeps only ASVs identified as bacterial
prefiltered_bac_ps<-subset_taxa(prefiltered_bac_ps, Kingdom=="k__Bacteria")
dim(otu_table(prefiltered_bac_ps)) # 28179 ASVs




# define plastid, mitochondria and host plant contamination ps objects
Mitochondria_ps<-subset_taxa(prefiltered_bac_ps, Family == "f__Mitochondria" | Family == "Mitochondria")
Plastid_ps<-subset_taxa(prefiltered_bac_ps, Order == "o__Chloroplast" | Order == "Chloroplast") 
host_plant_ps<-merge_phyloseq(Mitochondria_ps,Plastid_ps)

# quick histogram showing plant DNA contamination
hist(sample_sums(host_plant_ps)/sample_sums(prefiltered_bac_ps)*100)
# check average plant DNA contamination
mean(sample_sums(host_plant_ps)/sample_sums(prefiltered_bac_ps)*100)
# plastid/mitochondrial contamination is around 2.9 %.. we used PNA clamps to suppress it. very good.




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
prefiltered_bac_ps<-remove_Chloroplast_Mitochondria(prefiltered_bac_ps)
dim(otu_table(prefiltered_bac_ps)) # 27981 ASVs


# we now removed non-bacterial ASVs from "prefiltered_bac_ps" but did not remove rare ASVs. this dataset will be used for alpha-diversity, as the phyloseq::estimate_richness() function should be fed unfiltered data.



## 1.2.2 - Remove rare ASVs ====

# for all analyses other than alpha-diversity, we filter out rare ASVs as they might be artifacts
# note that we later repeat this filter after removing some samples

# first let's copy the phyloseq object
filtered_bac_ps <- prefiltered_bac_ps


# removes ASVs with less than 20 counts in the dataset
otu_table(filtered_bac_ps) <- otu_table(filtered_bac_ps)[which (rowSums(otu_table(filtered_bac_ps)) > 19),]
dim(otu_table(filtered_bac_ps)) # 10441 ASVs

# filter ASV according samples (more than 0 reads, present in 3 samples)
filter <- phyloseq::genefilter_sample(filtered_bac_ps, filterfun_sample(function(x) x > 0), A = 3) 
filtered_bac_ps <- prune_taxa(filter, filtered_bac_ps)
rm(filter) # we no longer need this string
dim(otu_table(filtered_bac_ps)) # 7674 ASVs







## 1.2.3 - Filter out extraction blanks ====

# check library sizes
sort(colSums(otu_table(filtered_bac_ps)))

# filter out extraction blank samples from filtered phyloseq object
filtered_bac_ps <- subset_samples(filtered_bac_ps, treatment != "-")
dim(otu_table(filtered_bac_ps)) # 7674 ASVs in 42 samples

# filter out extraction blank samples from "prefiltered" phyloseq object
prefiltered_bac_ps <- subset_samples(prefiltered_bac_ps, treatment != "-")
dim(otu_table(prefiltered_bac_ps)) # 27981 ASVs in 42 samples






## 1.2.4 - Filter out rare ASVs again ====

# filter the samples again with the same filter as in section 1.2.2 (min. 20 reads, present in min. 3 samples) as we removed invalid samples. this only removes 3 ASVs from bacteria and 11 ASVs from fungi.
filter <- phyloseq::genefilter_sample(filtered_bac_ps, filterfun_sample(function(x) x > 0), A = 3)
filtered_bac_ps <- prune_taxa(filter, filtered_bac_ps)
otu_table(filtered_bac_ps) <- otu_table(filtered_bac_ps)[which (rowSums(otu_table(filtered_bac_ps)) > 19),]


# check number of samples and ASVs left
dim(otu_table(filtered_bac_ps)) # 7672 ASVs, 42 samples




## 1.2.5 - Plot library sizes and host DNA contamination ====

# add library sizes as part of metadata
sample_data(filtered_bac_ps)$library_size <- sample_sums(filtered_bac_ps)

# check library size distribution
library_sizes <- sample_data(filtered_bac_ps)$library_size
hist(library_sizes)

summary(library_sizes) # quartiles, min. and max.

# samples ordered by library size
sample_data(filtered_bac_ps)[order(sample_data(filtered_bac_ps)$library_size),]

# run garbage collection after creating large objects
gc()




# plot of library sizes by host plant species and treatment. only library sizees below 1E6.
ggplot(data = sample_data(filtered_bac_ps), mapping = aes(x=library_size, y=plant_species, color = treatment, fill = treatment))+
  geom_jitter()+
  theme_classic()+
  geom_boxplot()+
  labs(x = "library size", y = "plant species")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))






## 1.2.6 - Cumulative sum scaling (CSS) transformation ====

# We use the Metagenomeseq package to be able to normalize library sizes without using rarefaction and also accounting for sparsity (high number of zeros in the dataset). This is done by considering the counts up to a certain quantile (comulative sum scaling, CSS). We will perform this with the  metagenomeSeq package.
# We will use CSS-normalized data for beta diversity analysis: ordinations, permanovas and beta dispersion


# first, let's transform the phyloseq object into an MR experiment object
MRexp_objt <- phyloseq_to_metagenomeSeq(filtered_bac_ps)

# # can be skipped because we can retrieve CSS-normalized data with MRcounts() below anyway:
# # normalize the object by cumulative sum scaling
# MRexp_objt_l <- lapply(MRexp_objt_l, function(x)
# cumNorm(x)) 

# here we can access the abundance matrix normalized by cumulative sum scaling. we use this table to fill in the new phyloseq object (see below)
CSS_matrix <- MRcounts(MRexp_objt, norm = TRUE, log = TRUE) # using a log scale will essentially reduce the impact of common species and increase the impact of rare species


# copy filtered phyloseq object to be transformed to CSS
CSS_bac_ps <- filtered_bac_ps

# and now change its taxa table
otu_table(CSS_bac_ps) <- otu_table(CSS_matrix, taxa_are_rows = TRUE)

# this is our CSS-transformed phyloseq object
CSS_bac_ps





## 1.2.7 - Exporting sequences for building phylogenetic tree of ASVs ====


# get refseq data into a dataframe in order to be imported in MAFFT to make the alignment and later on build a tree
# the resulting tree can then be used to compute UniFrac distances
# the tree will be built using qiime2, in the conda environment qiime-2022.2. see bryophytes_processing_3.


# sequences
refseq(filtered_bac_ps)

# export refseq df to folder "./data_bryophytes/16S"

filtered_bac_ps %>%
  refseq() %>%
  Biostrings::writeXStringSet("./data_bryophytes/16S/refseq_Bac.fna",
                              append=FALSE,compress=FALSE,
                              compression_level=NA, format="fasta")









## 1.2.9 - Save phyloseq objects ====

# make directory for phyloseq objects
dir.create("./data_bryophytes/16S/phyloseq_objects")

# save phyloseq objects
save(prefiltered_bac_ps, file = "./data_bryophytes/16S/phyloseq_objects/prefiltered_bac_ps.RData")
save(filtered_bac_ps, file = "./data_bryophytes/16S/phyloseq_objects/filtered_bac_ps.RData")
save(filtered_relative_bac_ps, file = "./data_bryophytes/16S/phyloseq_objects/filtered_relative_bac_ps.RData")
save(CSS_bac_ps, file = "./data_bryophytes/16S/phyloseq_objects/CSS_bac_ps.RData")











## 1.3 - ITS filtering ----


## 1.3.1 - ITS "pre-filtering": remove non-fungal ASVs ====


dim(otu_table(raw_fun_ps)) # input data: 17165 ASVs


# we keep the raw phyloseq object and create a new one that will be filtered
prefiltered_fun_ps <- raw_fun_ps

# histogram of amplicon lengths
hist(as.data.frame(refseq(prefiltered_fun_ps)@ranges)$width, breaks = 300)
summary(as.data.frame(refseq(prefiltered_fun_ps)@ranges)$width)



# take  all plant ITS sequences here to calculate host contamination
raw_PlantITS_ps<-subset_taxa(prefiltered_fun_ps, Kingdom=="k__Viridiplantae")
# note: no plant ASVs in filtered phyloseq object. there were some in the raw phyloseq object but they were filtered out already

# keeps only ASVs identified as fungal
prefiltered_fun_ps<-subset_taxa(prefiltered_fun_ps, Kingdom=="k__Fungi")
dim(otu_table(prefiltered_fun_ps)) # 11984 ASVs



# quick histogram showing plant DNA contamination
hist(sample_sums(raw_PlantITS_ps)/sample_sums(prefiltered_fun_ps)*100)
# check average plant DNA contamination
mean(sample_sums(raw_PlantITS_ps)/sample_sums(prefiltered_fun_ps)*100)
# plant DNA contamination is around 0.01 %. very good.


# we now removed non-fungal ASVs from "prefiltered_fun_ps" but did not remove rare ASVs. this dataset will be used for alpha-diversity, as the phyloseq::estimate_richness() function should be fed unfiltered data.




## 1.3.2 - Remove rare ASVs ====

# for all analyses other than alpha-diversity, we filter out rare ASVs as they might be artifacts
# note that we later repeat this filter after removing some samples (section 1.5)

# first let's copy the phyloseq object
filtered_fun_ps <- prefiltered_fun_ps


# removes taxa occuring less than 20 times in the dataset
otu_table(filtered_fun_ps) <- otu_table(filtered_fun_ps)[which (rowSums(otu_table(filtered_fun_ps)) > 19),]
dim(otu_table(filtered_fun_ps)) # 5951 ASVs

# filter  samples (min 0 reads, present in 3 samples)
filter <- phyloseq::genefilter_sample(filtered_fun_ps, filterfun_sample(function(x) x > 0), A = 3)
filtered_fun_ps <- prune_taxa(filter, filtered_fun_ps)
rm(filter) # we no longer need this string
dim(otu_table(filtered_fun_ps)) # 2684 ASVs



# check histogram of library sizes
hist(sample_sums(filtered_fun_ps))

# add library sizes as part of metadata
sample_data(filtered_fun_ps)$library_size<-sample_sums(filtered_fun_ps)


# remove unnecessary objects and collect garbage
rm(raw_PlantITS_ps)
gc()





## 1.3.3 - Filter out extraction blanks ====

# check library sizes
sort(colSums(otu_table(filtered_fun_ps)))

# filter out extraction blank samples from filtered phyloseq object
filtered_fun_ps <- subset_samples(filtered_fun_ps, treatment != "-")
dim(otu_table(filtered_fun_ps)) # 2684 ASVs in 42 samples

# filter out extraction blank samples from "prefiltered" phyloseq object
prefiltered_fun_ps <- subset_samples(prefiltered_fun_ps, treatment != "-")
dim(otu_table(prefiltered_fun_ps)) # 11984 ASVs in 42 samples






## 1.3.4 - Filter out rare ASVs again ====

# filter the samples again with the same filter as in section 1.2.2 (min. 20 reads, present in min. 3 samples) as we removed invalid samples. this only removes 3 ASVs from funteria and 11 ASVs from fungi.
filter <- phyloseq::genefilter_sample(filtered_fun_ps, filterfun_sample(function(x) x > 0), A = 3)
filtered_fun_ps <- prune_taxa(filter, filtered_fun_ps)
otu_table(filtered_fun_ps) <- otu_table(filtered_fun_ps)[which (rowSums(otu_table(filtered_fun_ps)) > 19),]


# check number of samples and ASVs left
dim(otu_table(filtered_fun_ps)) # 2660 ASVs, 42 samples




## 1.3.5 - Plot library sizes and host DNA contamination ====

# add library sizes as part of metadata
sample_data(filtered_fun_ps)$library_size <- sample_sums(filtered_fun_ps)

# check library size distribution
library_sizes <- sample_data(filtered_fun_ps)$library_size
hist(library_sizes)

summary(library_sizes) # quartiles, min. and max.

# samples ordered by library size
sample_data(filtered_fun_ps)[order(sample_data(filtered_fun_ps)$library_size),]

# run garbage collection after creating large objects
gc()




# plot of library sizes by host plant species and treatment. only library sizees below 1E6.
ggplot(data = sample_data(filtered_fun_ps), mapping = aes(x=library_size, y=plant_species, color = treatment, fill = treatment))+
  geom_jitter()+
  theme_classic()+
  geom_boxplot()+
  labs(x = "library size", y = "plant species")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))










## 1.3.6 - Cumulative sum scaling (CSS) transformation ====

# We use the Metagenomeseq package to be able to normalize library sizes without using rarefaction and also accounting for sparsity (high number of zeros in the dataset). This is done by considering the counts up to a certain quantile (comulative sum scaling, CSS). We will perform this with the  metagenomeSeq package.
# We will use CSS-normalized data for beta diversity analysis: ordinations, permanovas and beta dispersion


# first, let's transform the phyloseq object into an MR experiment object
MRexp_objt <- phyloseq_to_metagenomeSeq(filtered_fun_ps)

# # can be skipped because we can retrieve CSS-normalized data with MRcounts() below anyway:
# # normalize the object by cumulative sum scaling
# MRexp_objt_l <- lapply(MRexp_objt_l, function(x)
# cumNorm(x)) 

# here we can access the abundance matrix normalized by cumulative sum scaling. we use this table to fill in the new phyloseq object (see below)
CSS_matrix <- MRcounts(MRexp_objt, norm = TRUE, log = TRUE) # using a log scale will essentially reduce the impact of common species and increase the impact of rare species


# copy filtered phyloseq object to be transformed to CSS
CSS_fun_ps <- filtered_fun_ps

# and now change its taxa table
otu_table(CSS_fun_ps) <- otu_table(CSS_matrix, taxa_are_rows = TRUE)

# this is our CSS-transformed phyloseq object
CSS_fun_ps




## 1.3.7 - Exporting sequences for building phylogenetic tree of ASVs ====


# get refseq data into a dataframe in order to be imported in MAFFT to make the alignment and later on build a tree
# the resulting tree can then be used to compute UniFrac distances
# the tree will be built using qiime2, in the conda environment qiime-2022.2. see bryophytes_processing_3.


# sequences
refseq(filtered_fun_ps)

# export refseq df to folder "./data_bryophytes/ITS"

filtered_fun_ps %>%
  refseq() %>%
  Biostrings::writeXStringSet("./data_bryophytes/ITS/refseq_Fun.fna",
                              append=FALSE,compress=FALSE,
                              compression_level=NA, format="fasta")








## 1.3.9 - Save phyloseq objects ====

# make directory for phyloseq objects
dir.create("./data_bryophytes/ITS/phyloseq_objects")

# save phyloseq objects
save(prefiltered_fun_ps, file = "./data_bryophytes/ITS/phyloseq_objects/prefiltered_fun_ps.RData")
save(filtered_fun_ps, file = "./data_bryophytes/ITS/phyloseq_objects/filtered_fun_ps.RData")
save(filtered_relative_fun_ps, file = "./data_bryophytes/ITS/phyloseq_objects/filtered_relative_fun_ps.RData")
save(CSS_fun_ps, file = "./data_bryophytes/ITS/phyloseq_objects/CSS_fun_ps.RData")








