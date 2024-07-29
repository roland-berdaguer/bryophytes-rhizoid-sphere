## BRYOPHYTE RHIZOSPHERE AMPLICON DATA ANALYSIS SCRIPT

# 16S and ITS data from rhizoid-spheres of Marchantia polymorpha, Marchantia paleacea and Physcomitrium patens

# Note: the following abbreviations are used in the data:
# species: "mpol" = Marchantia polymorpha ; "mpal" = Marchantia paleacea ; "physco" = Physcomitrium patens
# treatments are abbreviated by approximate soil moisture contents:
# "20" = well-watered = 16.7% (w/w) soil moisture
# "13" = mild drought = 11.5% (w/w) soil moisture
# "10" = severe drought = 9.1% (w/w) soil moisture



setwd("C:/bryophytes_rhizosphere")


# 0 - Loading libraries ----

library(phyloseq)
library(vegan)
library(umap)
library(ggplot2)
library(pairwiseAdonis)
library(dplyr)
library(car)
library(tidyr)
library(DESeq2)
library(ape)
library(pheatmap)
library(dunn.test)
library(multcompView)



## 1 - Beta-diversity ----
# Fig. 1B-C, Supp. Table 1, Supp. Table 2


# load CSS-transformed physeq objects
load("./data_bryophytes/16S/phyloseq_objects/CSS_bac_ps.RData")
load("./data_bryophytes/ITS/phyloseq_objects/CSS_fun_ps.RData")


## 1.1 - 16S ----


## 1.1.1 - 16S UMAP with all samples (Fig. 1B) ====

# calculate Bray-Curtis dissimilarity on transposed otu table
bray_dist_all <- as.matrix(vegdist(t(otu_table(CSS_bac_ps)), method = "bray"))

# set seed for reproducibility
set.seed(100)

# run UMAP on Bray-Curtis dissimilarity matrix
umap_bray_all <- umap(bray_dist_all)

# extract "layout", which is the matrix in which we are interested
layout_bray_all <- as.data.frame((umap_bray_all)$layout)
colnames(layout_bray_all) <- c("UMAP1", "UMAP2")

# add sample metadata
umap_meta_bray_all <- cbind(layout_bray_all, sample_data(CSS_bac_ps))


# rename treatments for plot
umap_meta_bray_all$treatment <- gsub("20", "well-watered", umap_meta_bray_all$treatment)
umap_meta_bray_all$treatment <- gsub("13", "mild drought", umap_meta_bray_all$treatment)
umap_meta_bray_all$treatment <- gsub("10", "severe drought", umap_meta_bray_all$treatment)

# rename species for plot
umap_meta_bray_all$plant_species <- gsub("physco", "P. patens", umap_meta_bray_all$plant_species)
umap_meta_bray_all$plant_species <- gsub("mpal", "M. paleacea", umap_meta_bray_all$plant_species)
umap_meta_bray_all$plant_species <- gsub("mpol", "M. polymorpha", umap_meta_bray_all$plant_species)

# order treatments and plant species for plot
umap_meta_bray_all$treatment <- factor(umap_meta_bray_all$treatment, levels=c("well-watered", "mild drought", "severe drought"))
umap_meta_bray_all$plant_species <- factor(umap_meta_bray_all$plant_species, levels=c("P. patens", "M. paleacea", "M. polymorpha"))


# UMAP plot
ggplot(umap_meta_bray_all, aes(x=UMAP1,y=UMAP2, color=treatment, shape=plant_species)) +
  geom_point(aes(color=treatment), size=3, alpha=0.9) +
  theme_classic() +
  labs(color="Treatment",shape="Plant species")+
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18), axis.title=element_text(size=18)) + # text size for axis labels and legend
  # scale_color_manual(values = c("well-watered" = "#619CFF", "mild drought" = "#00BA38", "severe drought" = "#F8766D")) + # blue/green/red
  scale_color_manual(values = c("well-watered" = "#E69F00", "mild drought" = "#56B4E9", "severe drought" = "#009E73")) + # colour-blind friendly palette (from Wong 2011: https://www.nature.com/articles/nmeth.1618)
  scale_shape_discrete(labels = c(expression(italic("P. patens")), expression(italic("M. paleacea")), expression(italic("M. polymorpha")))) + # make species labels italic
  theme(legend.text.align = 0) # align legend text to the left (alignment automatically changed from using expression())

# export plot for Fig. 1B
ggsave("R_plots/bac_UMAP.png", width=2000, height=1300, unit="px")





## 1.1.2 - PERMANOVAs (Bray-Curtis) ----

## 1.1.2.1 - Overall PERMANOVAs (Supp. Table 1) ====

# first we do the PERMANOVA with all samples. However the data is unbalanced as the severe drought (SDR, 10%) only exists for M. polymorpha
# therefore we then also subset the data to only MDR (13%) & WW (20%) and only SDR (10%) & WW (20%), and run the PERMANOVA again with those

# rhizosphere all samples
# not a good test because the data is unbalanced
metadata <- as(sample_data(CSS_bac_ps), "data.frame")
distance_matrix <- phyloseq::distance(t(otu_table(CSS_bac_ps)), method="bray")
adonis2(distance_matrix ~ treatment*plant_species, permutations=999, data = metadata)
# adonis2(distance_matrix ~ plant_species*treatment, permutations=999, data = metadata) # with factors swapped, it barely changes the R2 values



## 1.1.2.2 - Pairwise Adonis (Supp. Table 2) ====

# make factors out of the variables plant_species and treatment
physeq_adonis <- CSS_bac_ps
metadata_physeq_adonis <- as(sample_data(physeq_adonis), "data.frame")
metadata_physeq_adonis$plant_species<-as.factor(metadata_physeq_adonis$plant_species)
metadata_physeq_adonis$treatment<-as.factor(metadata_physeq_adonis$treatment)


## Contrasting plant_species, for each treatment

# mild drought
pairwise.adonis2(phyloseq::distance((otu_table(subset_samples(physeq_adonis, treatment=="13"))), method="bray") ~ plant_species, data = metadata_physeq_adonis[which(metadata_physeq_adonis$treatment=="13"),])

# well-watered
pairwise.adonis2(phyloseq::distance((otu_table(subset_samples(physeq_adonis, treatment=="20"))), method="bray") ~ plant_species, data = metadata_physeq_adonis[which(metadata_physeq_adonis$treatment=="20"),])







## 1.2 - ITS ----



## - 1.2.1 UMAP with all samples (Fig. 1C) ====

# calculate Bray-Curtis dissimilarity on transposed otu table
bray_dist_all <- as.matrix(vegdist(t(otu_table(CSS_fun_ps)), method = "bray"))

# set seed for reproducibility
set.seed(100)

# run UMAP on Bray-Curtis dissimilarity matrix
umap_bray_all <- umap(bray_dist_all)

# extract "layout", which is the matrix in which we are interested
layout_bray_all <- as.data.frame((umap_bray_all)$layout)
colnames(layout_bray_all) <- c("UMAP1", "UMAP2")

# add sample metadata
umap_meta_bray_all <- cbind(layout_bray_all, sample_data(CSS_fun_ps))


# rename treatments for plot
umap_meta_bray_all$treatment <- gsub("20", "well-watered", umap_meta_bray_all$treatment)
umap_meta_bray_all$treatment <- gsub("13", "mild drought", umap_meta_bray_all$treatment)
umap_meta_bray_all$treatment <- gsub("10", "severe drought", umap_meta_bray_all$treatment)

# rename species for plot
umap_meta_bray_all$plant_species <- gsub("physco", "P. patens", umap_meta_bray_all$plant_species)
umap_meta_bray_all$plant_species <- gsub("mpal", "M. paleacea", umap_meta_bray_all$plant_species)
umap_meta_bray_all$plant_species <- gsub("mpol", "M. polymorpha", umap_meta_bray_all$plant_species)

# order treatments and plant species for plot
umap_meta_bray_all$treatment <- factor(umap_meta_bray_all$treatment, levels=c("well-watered", "mild drought", "severe drought"))
umap_meta_bray_all$plant_species <- factor(umap_meta_bray_all$plant_species, levels=c("P. patens", "M. paleacea", "M. polymorpha"))


# UMAP plot
ggplot(umap_meta_bray_all, aes(x=UMAP1,y=UMAP2, color=treatment, shape=plant_species)) +
  geom_point(aes(color=treatment), size=3, alpha=0.9) +
  theme_classic() +
  labs(color="Treatment",shape="Plant species")+
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18), axis.title=element_text(size=18)) + # text size for axis labels and legend
  # scale_color_manual(values = c("well-watered" = "#619CFF", "mild drought" = "#00BA38", "severe drought" = "#F8766D")) + # blue/green/red
  scale_color_manual(values = c("well-watered" = "#E69F00", "mild drought" = "#56B4E9", "severe drought" = "#009E73")) + # colour-blind friendly palette (from Wong 2011: https://www.nature.com/articles/nmeth.1618)
  scale_shape_discrete(labels = c(expression(italic("P. patens")), expression(italic("M. paleacea")), expression(italic("M. polymorpha")))) + # make species labels italic
  theme(legend.text.align = 0) # align legend text to the left (alignment automatically changed from using expression())


# export plot for Fig. 1C
ggsave("R_plots/fun_UMAP.png", width=2000, height=1300, unit="px")






## 1.2.2 - PERMANOVAs (Bray-Curtis) ----

## 1.2.2.1 - Overall PERMANOVA (Supp. Table 1) ====

# first we do the PERMANOVA with all samples. However the data is unbalanced as the severe drought (SDR, 10%) only exists for M. polymorpha
# therefore we then also subset the data to only MDR (13%) & WW (20%) and only SDR (10%) & WW (20%), and run the PERMANOVA again with those

# rhizosphere all samples
# not a good test because the data is unbalanced
metadata <- as(sample_data(CSS_fun_ps), "data.frame")
distance_matrix <- phyloseq::distance(t(otu_table(CSS_fun_ps)), method="bray")
adonis2(distance_matrix ~ treatment*plant_species, permutations=999, data = metadata)
# adonis2(distance_matrix ~ plant_species*treatment, permutations=999, data = metadata) # with factors swapped, it barely changes the R2 values



## 1.2.2.2 - Pairwise Adonis (Supp. Table 2) ====

# make factors out of the variables plant_species and treatment
physeq_adonis <- CSS_fun_ps
metadata_physeq_adonis <- as(sample_data(physeq_adonis), "data.frame")
metadata_physeq_adonis$plant_species<-as.factor(metadata_physeq_adonis$plant_species)
metadata_physeq_adonis$treatment<-as.factor(metadata_physeq_adonis$treatment)



## Contrasting plant_species, for each treatment

# mild drought
pairwise.adonis2(phyloseq::distance((otu_table(subset_samples(physeq_adonis, treatment=="13"))), method="bray") ~ plant_species, data = metadata_physeq_adonis[which(metadata_physeq_adonis$treatment=="13"),])

# well-watered
pairwise.adonis2(phyloseq::distance((otu_table(subset_samples(physeq_adonis, treatment=="20"))), method="bray") ~ plant_species, data = metadata_physeq_adonis[which(metadata_physeq_adonis$treatment=="20"),])






## 2 - Stacked bar plots at phylum level ----
# Fig. 1D-E

# define color-blind friendly palette
pal <- c("#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")


## 2.1 - 16S (Fig. 1D) ----

# load physeq object (the one with relative abundances)
load("./data_bryophytes/16S/phyloseq_objects/filtered_relative_bac_ps.RData")

# use psmelt to combine the ASV table with metadata and taxonomy
bac_melted <- psmelt(filtered_relative_bac_ps)


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
bac_phylum$treatment <- factor(bac_phylum$treatment, levels = c("20", "13", "10"))

# make species labeller for the plot, to have entire names instead of abbreviations
species_labeller <- as_labeller(c(`mpol` = "M. polymorpha",`mpal` = "M. paleacea", `physco` = "P. patens"))

# stacked bar plot by treatment for each plant species (as facets)
ggplot(bac_phylum, aes(fill=Phylum, y=Abundance, x=treatment)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(~factor(plant_species, levels=c("physco", "mpal", "mpol")), labeller = species_labeller, scales="free_x", space="free_x") +
  scale_x_discrete(labels = c("10" = "severe drought","13" = "mild drought","20" = "well-watered"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  xlab(NULL)+ # remove x axis label "treatment"
  ylab("relative abundance")+
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements ; angle x-axis labels
  scale_fill_manual(values=pal) +  # color-blind friendly palette
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic


# export plot for Fig. 1D
ggsave("R_plots/bac_phylum_barplot.png", width=2000, height=1300, unit="px")




## 2.2 - ITS (Fig. 1E) ----


# load physeq object (the one with relative abundances)
load("./data_bryophytes/ITS/phyloseq_objects/filtered_relative_fun_ps.RData")

# use psmelt to combine the ASV table with metadata and taxonomy
fun_melted <- psmelt(filtered_relative_fun_ps)


# group phyla
fun_phylum <- fun_melted %>% 
  group_by(Phylum, Sample, plant_species, treatment) %>%  # group by phylum and sample. also added plant_species and treatment as grouping factors in order to keep those columns
  summarise(Abundance = sum(Abundance))  # sum up the abundance of ASVs at phylum level

# remove "p__" in the phyla names
fun_phylum$Phylum <- gsub("p__", "", fun_phylum$Phylum)

# get the list of phyla by decreasing abundance, to order them
fun_phylum_order <- fun_melted %>% 
  group_by(Phylum) %>%
  summarise(Abundance = sum(Abundance))  # sum up the abundance for each phylum

# remove "p__" in the phyla names
fun_phylum_order$Phylum <- gsub("p__", "", fun_phylum_order$Phylum)
fun_phylum_order <- fun_phylum_order[order(fun_phylum_order$Abundance, decreasing=TRUE),]  # sort phyla by decreasing abundance


# let's move NA and "unidentified" to the bottom of the list so they don't get considered as one of the top 12 phyla (it should instead be in "others")
fun_phylum_order <- rbind(fun_phylum_order[!is.na(fun_phylum_order$Phylum),], fun_phylum_order[is.na(fun_phylum_order$Phylum),])
fun_phylum_order <- rbind(fun_phylum_order[!(fun_phylum_order$Phylum == "unidentified"),], fun_phylum_order[(fun_phylum_order$Phylum == "unidentified"),])


top_n_phyla <- 12


# replace phyla below the top N with "Other"
for (phylum_discard in fun_phylum_order$Phylum[(top_n_phyla+1):nrow(fun_phylum_order)]){
  fun_phylum$Phylum[fun_phylum$Phylum == phylum_discard] <- "Other"
}
# replace NA with "Other"
fun_phylum$Phylum[is.na(fun_phylum$Phylum)] <- "Other"

# preparing the ordered list of phyla with only the top N and "Other"
fun_phyla_list <- fun_phylum_order$Phylum[1:top_n_phyla]
fun_phyla_list <- append(fun_phyla_list, "Other")
# ordering the phyla
fun_phylum$Phylum <- factor(fun_phylum$Phylum, levels = fun_phyla_list)

# ordering the treatments for the plot
fun_phylum$treatment <- factor(fun_phylum$treatment, levels = c("20", "13", "10"))

# make species labeller for the plot, to have entire names instead of abbreviations
species_labeller <- as_labeller(c(`mpol` = "M. polymorpha",`mpal` = "M. paleacea", `physco` = "P. patens"))

# stacked bar plot by treatment for each plant species (as facets)
ggplot(fun_phylum, aes(fill=Phylum, y=Abundance, x=treatment)) + 
  geom_bar(position="fill", stat="identity") +
  facet_grid(~factor(plant_species, levels=c("physco", "mpal", "mpol")), labeller = species_labeller, scales="free_x", space="free_x") +
  scale_x_discrete(labels = c("10" = "severe drought","13" = "mild drought","20" = "well-watered"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  xlab(NULL)+ # remove x axis label "treatment"
  ylab("relative abundance")+
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements ; angle x-axis labels
  scale_fill_manual(values=pal) +  # color-blind friendly palette
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic


# export plot for Fig. 1E
ggsave("R_plots/fun_phylum_barplot.png", width=2000, height=1300, unit="px")





## 3 - Alpha diversity ----
# Supp. Fig. 1, Supp. Table 3


# load prefiltered physeq objects
load("./data_bryophytes/16S/phyloseq_objects/prefiltered_bac_ps.RData")
load("./data_bryophytes/ITS/phyloseq_objects/prefiltered_fun_ps.RData")


# make species labeller for the plots, to have entire names instead of abbreviations
species_labeller <- as_labeller(c(`mpol` = "M. polymorpha",`mpal` = "M. paleacea", `physco` = "P. patens"))


## 3.1 - 16S ----

## 3.1.1 - Create metadata table with alpha diversity metrics ====

# calculate different alpha diversity metrics
diversity<-estimate_richness(prefiltered_bac_ps)

# remove the X that appeared in the sample names for no reason
rownames(diversity) <- sub("X", "", rownames(diversity))

# add diversity metrics to sample_data of phyloseq object
sample_data(prefiltered_bac_ps) <- merge_phyloseq(prefiltered_bac_ps, sample_data(diversity)) # merge the new phyloseq object with the old phyloseq object

# extract diversity metrics incl. metadata as data.frame for statistical testing
diversity<-as(sample_data(prefiltered_bac_ps),"data.frame") # forces sample data of updated phyloseq object into a dataframe



## 3.1.2 - Boxplots of Shannon index (Supp. Fig. 1A) ====

# rename treatments for plot
diversity$treatment <- gsub("20", "well-watered", diversity$treatment)
diversity$treatment <- gsub("13", "mild drought", diversity$treatment)
diversity$treatment <- gsub("10", "severe drought", diversity$treatment)

# rename plant species for plot
diversity$plant_species <- gsub("physco", "P. patens", diversity$plant_species)
diversity$plant_species <- gsub("mpal", "M. paleacea", diversity$plant_species)
diversity$plant_species <- gsub("mpol", "M. polymorpha", diversity$plant_species)

# order treatments and plant species for plot
diversity$treatment <- factor(diversity$treatment, levels=c("well-watered", "mild drought", "severe drought"))
diversity$plant_species <- factor(diversity$plant_species, levels=c("P. patens", "M. paleacea", "M. polymorpha"))

# boxplots of Shannon index
ggplot(diversity, aes(x=treatment, y=Shannon)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = .2) +
  facet_wrap(~plant_species) +
  labs(title=NULL, x =NULL, y = "Shannon index") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements and angle x-axis labels
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic


# export plot for Supp. Fig. 1A
ggsave(filename="R_plots/bryophytes_bac_shannon.png", width=1400, height=1000, unit="px")



## 3.1.3 - Shannon index statistical testing (Supp. Table 3) ====

# before ANOVA we do a Levene's test for homoscedasticity = whether the variances are equal
# homogeneity of variances is an assumption of ANOVA
# about the ANOVA: considering that we have an unbalanced design, i.e. no severe drought for M. paleacea and P. patens, we will use a "type II" ANOVA. see: https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/


# test assumption for ANOVA: equal variances (homoscedasticity)
# Levene's test
leveneTest((Shannon) ~ treatment * plant_species, data = diversity)
# no significant heteroscedasticity


# test assumption for ANOVA: normal distribution of residuals
lm1 <- lm((Shannon) ~ treatment * plant_species, data = diversity)
residuals <- lm1$residuals # extract residuals
hist(residuals) # histogram of residuals
shapiro.test(residuals) # normality test for residuals
qqnorm(residuals) # QQ plot of residuals
qqline(residuals) # add line to QQ plot
# residuals are normally distributed

# two-way anova, type II because we have an unbalanced design
Anova(lm1, type = "2")
# no significant effects







## 3.2 - ITS ----


## 3.2.1 Create metadata table with alpha diversity metrics ====

# calculate different alpha diversity metrics
diversity<-estimate_richness(prefiltered_fun_ps)

# remove the X that appeared in the sample names for no reason
rownames(diversity) <- sub("X", "", rownames(diversity))

# add diversity metrics to sample_data of phyloseq object
sample_data(prefiltered_fun_ps) <- merge_phyloseq(prefiltered_fun_ps, sample_data(diversity)) # merge the new phyloseq object with the old phyloseq object

# extract diversity metrics incl. metadata as data.frame for statistical testing
diversity<-as(sample_data(prefiltered_fun_ps),"data.frame") # forces sample data of updated phyloseq object into a dataframe


## 3.2.2 - Boxplots of Shannon index (Supp. Fig. 1B) ====

# rename treatments for plot
diversity$treatment <- gsub("20", "well-watered", diversity$treatment)
diversity$treatment <- gsub("13", "mild drought", diversity$treatment)
diversity$treatment <- gsub("10", "severe drought", diversity$treatment)

# rename plant species for plot
diversity$plant_species <- gsub("physco", "P. patens", diversity$plant_species)
diversity$plant_species <- gsub("mpal", "M. paleacea", diversity$plant_species)
diversity$plant_species <- gsub("mpol", "M. polymorpha", diversity$plant_species)

# order treatments and plant species for plot
diversity$treatment <- factor(diversity$treatment, levels=c("well-watered", "mild drought", "severe drought"))
diversity$plant_species <- factor(diversity$plant_species, levels=c("P. patens", "M. paleacea", "M. polymorpha"))

# boxplots of Shannon index
ggplot(diversity, aes(x=treatment, y=Shannon)) + 
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~plant_species) +
  labs(title=NULL, x =NULL, y = "Shannon index") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements and angle x-axis labels
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic





## 3.2.3 - Shannon index statistical testing (Supp. Table 3) ====

## 3.2.3.1 - ANOVA ####

# before ANOVA we do a Levene's test for homoscedasticity = whether the variances are equal
# homogeneity of variances and is an assumption of ANOVA
# we also test for normal distribution of residuals, another assumption of ANOVA
# about the ANOVA: considering that we have an unbalanced design, i.e. no severe drought for M. paleacea and P. patens, we will use a "type II" ANOVA. see: https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/


# test assumption for ANOVA: equal variances (homoscedasticity)
# Levene's test
leveneTest((Shannon) ~ treatment * plant_species, data = diversity)
# no significant heteroscedasticity

# test assumption for ANOVA: normal distribution of residuals
lm1 <- lm((Shannon) ~ treatment * plant_species, data = diversity)
residuals <- lm1$residuals # extract residuals
hist(residuals) # histogram of residuals
shapiro.test(residuals) # normality test for residuals
qqnorm(residuals) # QQ plot of residuals
qqline(residuals) # add line to QQ plot
# residuals are normally distributed

# two-way anova, type II because we have an unbalanced design
Anova(lm1, type = "2")
# significant treatment, plant_species and interaction effects





## 3.2.2 - Boxplots of Shannon index with significance group letters based on Tukey test (Supp. Fig. 1B) ====

# make column species/treatment
diversity$species_treatment <- paste0(diversity$plant_species, " / ", diversity$treatment)


# replace species/treatments with letter to avoid problems with an extra space with the multcompLetters() function
# alternatively, you can just make sure your group names don't have any spaces or hyphens
diversity$species_treatment <- gsub("P. patens / well-watered", "A", diversity$species_treatment)
diversity$species_treatment <- gsub("P. patens / mild drought", "B", diversity$species_treatment)
diversity$species_treatment <- gsub("M. paleacea / well-watered", "C", diversity$species_treatment)
diversity$species_treatment <- gsub("M. paleacea / mild drought", "D", diversity$species_treatment)
diversity$species_treatment <- gsub("M. polymorpha / well-watered", "E", diversity$species_treatment)
diversity$species_treatment <- gsub("M. polymorpha / mild drought", "F", diversity$species_treatment)
diversity$species_treatment <- gsub("M. polymorpha / severe drought", "G", diversity$species_treatment)


# one-way ANOVA
anova <- aov(Shannon ~ species_treatment, data = diversity)
# check p-value of ANOVA
summary(anova)[[1]][["Pr(>F)"]][1]
# Tukey HSD test
tukey <- TukeyHSD(anova)

# put p-values in a named vector for multcompLetters()
pvec <- tukey$species_treatment[,"p adj"]

## re-order treatment combinations to get the letters starting with A
group_order <- c("A", "B", "C", "D", "E", "F", "G")
# list all combinations of groups based on our previously defined group order. this will be useful to get significance letters going alphabetically in the order of our groups (so the first group starts with significance letter "a")
group_combinations <- combn(group_order, 2, paste0, collapse = "-")
# some group combinations are written the other way around in names(pvec) compared to group_combinations .. let's flip those around
for (i in which(names(pvec) %in% setdiff(names(pvec), group_combinations))){
  names(pvec)[i] <- paste0(strsplit(names(pvec)[i], "-")[[1]][2], "-", strsplit(names(pvec)[i], "-")[[1]][1])
}
# re-order pvec in the defined order
pvec <- pvec[group_combinations]

# compute significance letters
cld <- multcompLetters(pvec, threshold=0.05)$Letters


## make a dataframe with significance letters. it's a bit tricky because it's a faceted boxplot
# first we turn the named vector with significance letters into a dataframe
letters <- data.frame(cld)
# make column "group"
letters$group <- rownames(letters)
# re-order the dataframe so from A to G
letters <- letters[c("A", "B", "C", "D", "E", "F", "G"),]
# add columns plant_species and treatment
letters$plant_species <- factor(c("P. patens", "P. patens", "M. paleacea", "M. paleacea", "M. polymorpha", "M. polymorpha", "M. polymorpha"))
letters$treatment <- factor(c("well-watered", "mild drought", "well-watered", "mild drought", "well-watered", "mild drought", "severe drought"))


# boxplots of Shannon index
p <- ggplot(diversity, aes(x=treatment, y=Shannon)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = .2) +
  facet_wrap(~plant_species) +
  labs(title=NULL, x =NULL, y = "Shannon index") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements and angle x-axis labels
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic


# define y-positions for labels: extract maximum y of each boxplot from the ggplot object
letters$ypos <- layer_data(p)$ymax
# get y-axis range
yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
# add 5% of the y axis height to have the label slightly higher
letters$ypos <- letters$ypos + 0.1 * yrange


# plot plus significance labels
p + geom_text(data = letters, aes(x = treatment, y = ypos, label = cld)) # significance letters

# export plot for Supp. Fig. 1B
ggsave(filename="R_plots/bryophytes_fun_shannon.png", width=1400, height=1000, unit="px")



## 4 - Boxplots of relative abundance of specific ASVs ----
# Fig. 2A-B, Fig. 3A&C, Supp. Fig. 7


## 4.1 - 16S (Fig. 2A-B) ====

# load physeq object (the one with relative abundances)
load("./data_bryophytes/16S/phyloseq_objects/filtered_relative_bac_ps.RData")

# use psmelt to combine the ASV table with metadata and taxonomy
bac_melted <- psmelt(filtered_relative_bac_ps)


# function to plot a chosen taxon
# boxplots of relative abundance by treatment, for each species (facets)
# first argument: taxonomic level -> can be "Phylum", "Class", "Order", "Family", "Genus", "Species" or "ASV"
# second argument: the chosen taxon. examples: "p__Actinobacteriota", "o__Streptomycetales", etc.. note there are 2 underscores
# if the taxonomic level is "ASV", the chosen taxon should be for example "ASV_1" (16S) or "fASV_1", (ITS)
# third argument: chosen plant species -> use abbreviation (species code). for all species, use "all"


# change plant_species and treatment to character (needed for following code)
bac_melted$plant_species <- as.character(bac_melted$plant_species)
bac_melted$treatment <- as.character(bac_melted$treatment)

# change species names to full names
bac_melted$plant_species[which(bac_melted$plant_species=="mpol")] <- "M. polymorpha"
bac_melted$plant_species[which(bac_melted$plant_species=="mpal")] <- "M. paleacea"
bac_melted$plant_species[which(bac_melted$plant_species=="physco")] <- "P. patens"

# change treatment names to words
bac_melted$treatment <- as.character(bac_melted$treatment)
bac_melted$treatment[which(bac_melted$treatment=="20")] <- "well-watered"
bac_melted$treatment[which(bac_melted$treatment=="13")] <- "mild drought"
bac_melted$treatment[which(bac_melted$treatment=="10")] <- "severe drought"

# order treatments and species to have them in chosen order in boxplots
bac_melted$treatment <- factor(bac_melted$treatment, levels = c("well-watered", "mild drought", "severe drought"))
bac_melted$plant_species <- factor(bac_melted$plant_species, levels = c("P. patens", "M. paleacea", "M. polymorpha"))


# # define function for faceted boxplots
# plot_taxon_bac <- function(taxonomic_level, chosen_taxon, chosen_species="all"){
#   x <- bac_melted  # use melted dataframe of 16S
#   
#   if(taxonomic_level=="ASV"){x$taxon <- x$OTU} else{  # if "ASV" is chosen, we don't need the group_by function, so we copy the "OTU" column in "taxon" column
#     x <- x %>% 
#       group_by(taxon=get(taxonomic_level), Sample, plant_species, treatment) %>%  # group by taxonomic level and sample. also added plant_species and treatment as grouping factors in order to keep those columns
#       summarise(Abundance = sum(Abundance))  # sum up the abundance of all ASVs at a chosen taxonomic level
#   }
#   
#   x <- subset(x, taxon==chosen_taxon)  # subset only rows with the chosen taxon
#   if(chosen_species!="all"){
#     x <- subset(x, plant_species==chosen_species)}  # subset only rows with the chosen taxon
#   
#   ggplot(data = x, aes(x = treatment, y = Abundance)) +
#     geom_boxplot(outlier.shape = NA) +  # hide the outliers from the boxplot because we see them with geom_jitter anyway
#     geom_jitter(height = 0, width = .2) +  # geom_jitter shows all the data points
#     labs(x = "", y = "relative abundance") +
#     facet_grid(~ plant_species, scales = "free") +  # scales="free_y" indicates that the y-axis scale can be different in every plot. scales="fixed" gives a fixed scale
#     theme_bw() +
#     theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements ; angle x-axis labels
#     theme(legend.position = "none") + # remove legend
#     expand_limits(y = 0) + # force y-axis to show 0
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
#     theme(strip.background = element_blank()) + # no box around facet label
#     theme(strip.text = element_text(face = "italic")) # facet label in italic
#   
#   
# }



# define function for faceted boxplots with significance letters
plot_taxon_bac <- function(taxonomic_level, chosen_taxon, chosen_species="all"){
  x <- bac_melted  # use melted dataframe of 16S
  
  if(taxonomic_level=="ASV"){x$taxon <- x$OTU} else{  # if "ASV" is chosen, we don't need the group_by function, so we copy the "OTU" column in "taxon" column
    x <- x %>% 
      group_by(taxon=get(taxonomic_level), Sample, plant_species, treatment) %>%  # group by taxonomic level and sample. also added plant_species and treatment as grouping factors in order to keep those columns
      summarise(Abundance = sum(Abundance))  # sum up the abundance of all ASVs at a chosen taxonomic level
  }
  
  x <- subset(x, taxon==chosen_taxon)  # subset only rows with the chosen taxon
  if(chosen_species!="all"){
    x <- subset(x, plant_species==chosen_species)}  # subset only rows with the chosen taxon
  
  # make column species/treatment
  x$species_treatment <- paste0(x$plant_species, " / ", x$treatment)
  
  
  # replace species/treatments with letter to avoid problems with an extra space with the multcompLetters() function
  # alternatively, you can just make sure your group names don't have any spaces or hyphens
  x$species_treatment <- gsub("P. patens / well-watered", "A", x$species_treatment)
  x$species_treatment <- gsub("P. patens / mild drought", "B", x$species_treatment)
  x$species_treatment <- gsub("M. paleacea / well-watered", "C", x$species_treatment)
  x$species_treatment <- gsub("M. paleacea / mild drought", "D", x$species_treatment)
  x$species_treatment <- gsub("M. polymorpha / well-watered", "E", x$species_treatment)
  x$species_treatment <- gsub("M. polymorpha / mild drought", "F", x$species_treatment)
  x$species_treatment <- gsub("M. polymorpha / severe drought", "G", x$species_treatment)
  
  
  
  # Kruskal-Wallis test
  kw <- kruskal.test(Abundance ~ species_treatment, data = x)
  # check p-value of Kruskal-Wallis
  print(kw$p.value)
  print(kw$p.value < 0.05)
  # Dunn's test
  dunn <- dunn.test(x$Abundance, x$species_treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)
  # put p-values in a named vector for multcompLetters()
  dunnvec <- dunn$altP.adjusted
  names(dunnvec) <- dunn$comparisons
  # remove spaces because an extra space is messing up multcompLetters()
  names(dunnvec) <- gsub(" ", "", names(dunnvec))
  # compute significance letters
  cld <- multcompLetters(dunnvec, threshold=0.05)$Letters
  
  
  ## make a dataframe with significance letters. it's a bit tricky because it's a faceted boxplot
  # first we turn the named vector with significance letters into a dataframe
  letters <- data.frame(cld)
  # make column "group"
  letters$group <- rownames(letters)
  # re-order the dataframe so from A to G
  letters <- letters[c("A", "B", "C", "D", "E", "F", "G"),]
  # add columns plant_species and treatment
  letters$plant_species <- factor(c("P. patens", "P. patens", "M. paleacea", "M. paleacea", "M. polymorpha", "M. polymorpha", "M. polymorpha"))
  letters$treatment <- factor(c("well-watered", "mild drought", "well-watered", "mild drought", "well-watered", "mild drought", "severe drought"))
  
  # make ggplot object
  p <- ggplot(data = x, aes(x = treatment, y = Abundance)) +
    geom_boxplot(outlier.shape = NA) +  # hide the outliers from the boxplot because we see them with geom_jitter anyway
    geom_jitter(height = 0, width = .2) +  # geom_jitter shows all the data points
    labs(x = "", y = "relative abundance") +
    facet_grid(~ plant_species, scales = "free") +  # scales="free_y" indicates that the y-axis scale can be different in every plot. scales="fixed" gives a fixed scale
    # theme_bw() +
    theme_classic() +
    # theme(panel.border = element_blank()) +
    # theme(panel.background = element_rect(color = "black")) +
    theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements ; angle x-axis labels
    theme(legend.position = "none") + # remove legend
    expand_limits(y = 0) + # force y-axis to show 0
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
    theme(strip.background = element_blank()) + # no box around facet label
    theme(strip.text = element_text(face = "italic")) # facet label in italic
  
  # define y-positions for labels: extract maximum y of each boxplot from the ggplot object
  letters$ypos <- layer_data(p)$ymax
  # get y-axis range
  yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
  # add 5% of the y axis height to have the label slightly higher
  letters$ypos <- letters$ypos + 0.05 * yrange
  
  # if the Kruskal-Wallis test was significant
  if(kw$p.value < 0.05){
    # plot plus significance labels
    p + geom_text(data = letters, aes(x = treatment, y = ypos, label = cld)) # significance letters
  } else{
    # plot without significance labels
    p
  }
  
}





## FOR FIGURES

# ASV_133 (Paenibacillus)
plot_taxon_bac("ASV", "ASV_133")
# export plot for Fig. 2A
ggsave("R_plots/Paenibacillus_ASV_133.png", width=1600, height=1300, unit="px")

# ASV_150 (Paenibacillus)
plot_taxon_bac("ASV", "ASV_150")
# export plot for Fig. 2B
ggsave("R_plots/Paenibacillus_ASV_150.png", width=1600, height=1300, unit="px")









## 4.2 - ITS (Fig. 3A&C, Supp. Fig. 7) ====

# load physeq object (the one with relative abundances)
load("./data_bryophytes/ITS/phyloseq_objects/filtered_relative_fun_ps.RData")

# use psmelt to combine the ASV table with metadata and taxonomy
fun_melted <- psmelt(filtered_relative_fun_ps)


# function to plot a chosen taxon
# boxplots of relative abundance by treatment, for each species (facets)
# first argument: taxonomic level -> can be "Phylum", "Class", "Order", "Family", "Genus", "Species" or "ASV"
# second argument: the chosen taxon. examples: "p__Actinobacteriota", "o__Streptomycetales", etc.. note there are 2 underscores
# if the taxonomic level is "ASV", the chosen taxon should be for example "ASV_1" (16S) or "fASV_1", (ITS)
# third argument: chosen plant species -> use abbreviation (species code). for all species, use "all"


# change plant_species and treatment to character (needed for following code)
fun_melted$plant_species <- as.character(fun_melted$plant_species)
fun_melted$treatment <- as.character(fun_melted$treatment)

# change species names to full names
fun_melted$plant_species[which(fun_melted$plant_species=="mpol")] <- "M. polymorpha"
fun_melted$plant_species[which(fun_melted$plant_species=="mpal")] <- "M. paleacea"
fun_melted$plant_species[which(fun_melted$plant_species=="physco")] <- "P. patens"

# change treatment names to words
fun_melted$treatment <- as.character(fun_melted$treatment)
fun_melted$treatment[which(fun_melted$treatment=="20")] <- "well-watered"
fun_melted$treatment[which(fun_melted$treatment=="13")] <- "mild drought"
fun_melted$treatment[which(fun_melted$treatment=="10")] <- "severe drought"

# order treatments and species to have them in chosen order in boxplots
fun_melted$treatment <- factor(fun_melted$treatment, levels = c("well-watered", "mild drought", "severe drought"))
fun_melted$plant_species <- factor(fun_melted$plant_species, levels = c("P. patens", "M. paleacea", "M. polymorpha"))


# # define function for faceted boxplots
# plot_taxon_fun <- function(taxonomic_level, chosen_taxon, chosen_species="all"){
#   x <- fun_melted  # use melted dataframe of 16S
#   
#   if(taxonomic_level=="ASV"){x$taxon <- x$OTU} else{  # if "ASV" is chosen, we don't need the group_by function, so we copy the "OTU" column in "taxon" column
#     x <- x %>% 
#       group_by(taxon=get(taxonomic_level), Sample, plant_species, treatment) %>%  # group by taxonomic level and sample. also added plant_species and treatment as grouping factors in order to keep those columns
#       summarise(Abundance = sum(Abundance))  # sum up the abundance of all ASVs at a chosen taxonomic level
#   }
#   
#   x <- subset(x, taxon==chosen_taxon)  # subset only rows with the chosen taxon
#   if(chosen_species!="all"){
#     x <- subset(x, plant_species==chosen_species)}  # subset only rows with the chosen taxon
#   
#   ggplot(data = x, aes(x = treatment, y = Abundance)) +
#     geom_boxplot(outlier.shape = NA) +  # hide the outliers from the boxplot because we see them with geom_jitter anyway
#     geom_jitter(height = 0, width = .2) +  # geom_jitter shows all the data points
#     labs(x = "", y = "relative abundance") +
#     facet_grid(~ plant_species, scales = "free") +  # scales="free_y" indicates that the y-axis scale can be different in every plot. scales="fixed" gives a fixed scale
#     theme_bw() +
#     theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements ; angle x-axis labels
#     theme(legend.position = "none") + # remove legend
#     expand_limits(y = 0) + # force y-axis to show 0
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
#     theme(strip.background = element_blank()) + # no box around facet label
#     theme(strip.text = element_text(face = "italic")) # facet label in italic
#   
# }





# define function for faceted boxplots with significance letters
plot_taxon_fun <- function(taxonomic_level, chosen_taxon, chosen_species="all"){
  x <- fun_melted  # use melted dataframe of 16S
  
  if(taxonomic_level=="ASV"){x$taxon <- x$OTU} else{  # if "ASV" is chosen, we don't need the group_by function, so we copy the "OTU" column in "taxon" column
    x <- x %>% 
      group_by(taxon=get(taxonomic_level), Sample, plant_species, treatment) %>%  # group by taxonomic level and sample. also added plant_species and treatment as grouping factors in order to keep those columns
      summarise(Abundance = sum(Abundance))  # sum up the abundance of all ASVs at a chosen taxonomic level
  }
  
  x <- subset(x, taxon==chosen_taxon)  # subset only rows with the chosen taxon
  if(chosen_species!="all"){
    x <- subset(x, plant_species==chosen_species)}  # subset only rows with the chosen taxon
  
  # make column species/treatment
  x$species_treatment <- paste0(x$plant_species, " / ", x$treatment)
  
  
  # replace species/treatments with letter to avoid problems with an extra space with the multcompLetters() function
  # alternatively, you can just make sure your group names don't have any spaces or hyphens
  x$species_treatment <- gsub("P. patens / well-watered", "A", x$species_treatment)
  x$species_treatment <- gsub("P. patens / mild drought", "B", x$species_treatment)
  x$species_treatment <- gsub("M. paleacea / well-watered", "C", x$species_treatment)
  x$species_treatment <- gsub("M. paleacea / mild drought", "D", x$species_treatment)
  x$species_treatment <- gsub("M. polymorpha / well-watered", "E", x$species_treatment)
  x$species_treatment <- gsub("M. polymorpha / mild drought", "F", x$species_treatment)
  x$species_treatment <- gsub("M. polymorpha / severe drought", "G", x$species_treatment)
  
  
  
  # Kruskal-Wallis test
  kw <- kruskal.test(Abundance ~ species_treatment, data = x)
  # check p-value of Kruskal-Wallis
  print(kw$p.value)
  print(kw$p.value < 0.05)
  # Dunn's test
  dunn <- dunn.test(x$Abundance, x$species_treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)
  # put p-values in a named vector for multcompLetters()
  dunnvec <- dunn$altP.adjusted
  names(dunnvec) <- dunn$comparisons
  # remove spaces because an extra space is messing up multcompLetters()
  names(dunnvec) <- gsub(" ", "", names(dunnvec))
  # compute significance letters
  cld <- multcompLetters(dunnvec, threshold=0.05)$Letters
  
  
  ## make a dataframe with significance letters. it's a bit tricky because it's a faceted boxplot
  # first we turn the named vector with significance letters into a dataframe
  letters <- data.frame(cld)
  # make column "group"
  letters$group <- rownames(letters)
  # re-order the dataframe so from A to G
  letters <- letters[c("A", "B", "C", "D", "E", "F", "G"),]
  # add columns plant_species and treatment
  letters$plant_species <- factor(c("P. patens", "P. patens", "M. paleacea", "M. paleacea", "M. polymorpha", "M. polymorpha", "M. polymorpha"))
  letters$treatment <- factor(c("well-watered", "mild drought", "well-watered", "mild drought", "well-watered", "mild drought", "severe drought"))
  
  # make ggplot object
  p <- ggplot(data = x, aes(x = treatment, y = Abundance)) +
    geom_boxplot(outlier.shape = NA) +  # hide the outliers from the boxplot because we see them with geom_jitter anyway
    geom_jitter(height = 0, width = .2) +  # geom_jitter shows all the data points
    labs(x = "", y = "relative abundance") +
    facet_grid(~ plant_species, scales = "free") +  # scales="free_y" indicates that the y-axis scale can be different in every plot. scales="fixed" gives a fixed scale
    theme_classic() +
    theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # increase text size for all text elements ; angle x-axis labels
    theme(legend.position = "none") + # remove legend
    expand_limits(y = 0) + # force y-axis to show 0
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
    theme(strip.background = element_blank()) + # no box around facet label
    theme(strip.text = element_text(face = "italic")) # facet label in italic
  
  # define y-positions for labels: extract maximum y of each boxplot from the ggplot object
  letters$ypos <- layer_data(p)$ymax
  # get y-axis range
  yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
  # add 5% of the y axis height to have the label slightly higher
  letters$ypos <- letters$ypos + 0.05 * yrange
  
  
  # if the Kruskal-Wallis test was significant
  if(kw$p.value < 0.05){
    # plot plus significance labels
    p + geom_text(data = letters, aes(x = treatment, y = ypos, label = cld)) # significance letters
  } else{
    # plot without significance labels
    p
  }

  
}







## FOR FIGURES

# Olpidium
plot_taxon_fun("Genus", "g__Olpidium")
# export plot for Fig. 3A
ggsave("R_plots/Olpidium.png", width=1600, height=1300, unit="px")

# Glomeromycota
plot_taxon_fun("Phylum", "p__Glomeromycota")
# export plot for Fig. 3C
ggsave("R_plots/Glomeromycota.png", width=1600, height=1300, unit="px")



## Fusarium (Supp. Fig. 7)

# fASV_19
plot_taxon_fun("ASV", "fASV_19")
# export plot for Fig. 7A
ggsave("R_plots/Fusarium_fASV_19.png", width=1600, height=1300, unit="px")

# fASV_45
plot_taxon_fun("ASV", "fASV_45")
# export plot for Fig. 7B
ggsave("R_plots/Fusarium_fASV_45.png", width=1600, height=1300, unit="px")

# fASV_50
plot_taxon_fun("ASV", "fASV_50")
# export plot for Fig. 7C
ggsave("R_plots/Fusarium_fASV_50.png", width=1600, height=1300, unit="px")

# fASV_22
plot_taxon_fun("ASV", "fASV_22")
# export plot for Fig. 7D
ggsave("R_plots/Fusarium_fASV_22.png", width=1600, height=1300, unit="px")

# fASV_59
plot_taxon_fun("ASV", "fASV_59")
# export plot for Fig. 7E
ggsave("R_plots/Fusarium_fASV_59.png", width=1600, height=1300, unit="px")








## 5 - Differential abundance at phylum level ----
# Fig. 4B, Supp. Fig. 2


# using the Mann-Whitney (Wilcoxon) test

## 5.1 - 16S (Fig. 4B) ----

alpha <- 0.05


## 5.1.1 - Wilcoxon test ====


# separate to each plant species
filtered_relative_bac_ps_species_l <- metagMisc::phyloseq_sep_variable(filtered_relative_bac_ps, variable = "plant_species")


# agglomerate at phylum level
# note that the phyla still have ASV IDs (of a random ASV corresponding to the phylum) as names, instead of phylum names.. we change this after getting the results tables from DESeq2
bac_phyla_rarefied_relative_ps_species_l <- lapply(filtered_relative_bac_ps_species_l, function(x) tax_glom(x, taxrank = "Phylum"))



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
    abundance_WW <- phyla_avg$abundance[which(phyla_avg$treatment == "20" & phyla_avg$phylum == p)]
    abundance_MDR <- phyla_avg$abundance[which(phyla_avg$treatment == "13" & phyla_avg$phylum == p)]
    results[which(results$phylum == p), "LFC_MDR_WW"] <- log2(abundance_MDR / abundance_WW)
  
    # p-value by Wilcoxon test
    values_WW <- subset(phyla_table, (treatment == "20" & phylum == p))$abundance
    values_MDR <- subset(phyla_table, (treatment == "13" & phylum == p))$abundance
    results[which(results$phylum == p), "pvalue_MDR_WW"] <- wilcox.test(values_MDR, values_WW)$p.value
  
  
    # calculate SDR vs. WW LFC and p-value only for mpol
    if(phyla_table$plant_species[1] %in% c("mpol")){
      # LFC
      abundance_SDR <- phyla_avg$abundance[which(phyla_avg$treatment == "10" & phyla_avg$phylum == p)]
      results[which(results$phylum == p), "LFC_SDR_WW"] <- log2(abundance_SDR / abundance_WW)
  
      # p-value by Wilcoxon test
      values_SDR <- subset(phyla_table, (treatment == "10" & phylum == p))$abundance
      results[which(results$phylum == p), "pvalue_SDR_WW"] <- wilcox.test(values_SDR, values_WW)$p.value
    }
  
  
  }
  
  return(results)

})



# subset top 12 phyla.  note that they are listed by order decreasing abundance
top_12_phyla <- c("Actinobacteriota", "Proteobacteria", "Firmicutes", "Acidobacteriota", "Chloroflexi", "Planctomycetota", "Gemmatimonadota", "Verrucomicrobiota", "Bacteroidota", "Myxococcota", "Methylomirabilota", "Cyanobacteria" )
phyla_results <- lapply(phyla_results, function(x) x[which(x$phylum %in% top_12_phyla),])

## apply FDR-correction to the p-values
# MDR
phyla_results <- lapply(phyla_results, function(x){
  x$pvalue_MDR_WW <- stats::p.adjust(p = x$pvalue_MDR_WW, method = "fdr")
  return(x)
})
# SDR (only M. polymorpha)
phyla_results$mpol$pvalue_SDR_WW <- stats::p.adjust(p = phyla_results$mpol$pvalue_SDR_WW, method = "fdr")



## 5.1.2 - Merged dataframe for heatmap ====

## MDR (mild drought)
res_MDR_WW_mpol <- data.frame(phylum = phyla_results$mpol$phylum, log2FoldChange = phyla_results$mpol$LFC_MDR_WW, padj = phyla_results$mpol$pvalue_MDR_WW)
res_MDR_WW_mpol$plant_species <- "M. polymorpha"

res_MDR_WW_mpal <- data.frame(phylum = phyla_results$mpal$phylum, log2FoldChange = phyla_results$mpal$LFC_MDR_WW, padj = phyla_results$mpal$pvalue_MDR_WW)
res_MDR_WW_mpal$plant_species <- "M. paleacea"

res_MDR_WW_physco <- data.frame(phylum = phyla_results$physco$phylum, log2FoldChange = phyla_results$physco$LFC_MDR_WW, padj = phyla_results$physco$pvalue_MDR_WW)
res_MDR_WW_physco$plant_species <- "P. patens"

# merge all the results tables
res_MDR_WW <- rbind(res_MDR_WW_mpol, res_MDR_WW_mpal, res_MDR_WW_physco)

# add column "treatment"
res_MDR_WW$treatment <- "mild drought"



## SDR (severe drought)
res_SDR_WW_mpol <- data.frame(phylum = phyla_results$mpol$phylum, log2FoldChange = phyla_results$mpol$LFC_SDR_WW, padj = phyla_results$mpol$pvalue_SDR_WW)
res_SDR_WW_mpol$plant_species <- "M. polymorpha"

# add column "treatment"
res_SDR_WW_mpol$treatment <- "severe drought"



## merge MDR and SDR
res_DR_WW_merged <- rbind(res_MDR_WW, res_SDR_WW_mpol)



# add column "significance" with asterisks for different significance levels (based on adjusted p-value)
# note that this is the adjusted p-value from testing one treatment (MDR or SDR) vs. WW in one plant species. so it is only adjusted for multiple testing with several phyla, not several treatments or species
res_DR_WW_merged$significance <- ""
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.05)] <- "*"
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.01)] <- "**"
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.001)] <- "***"





## 5.1.3 - Heatmap ====

# set LFC limits to -1 and +1 manually by replacing values under -1 with -1 and over +1 with 1
res_DR_WW_merged$log2FoldChange[which(res_DR_WW_merged$log2FoldChange < -1)] <- -1
res_DR_WW_merged$log2FoldChange[which(res_DR_WW_merged$log2FoldChange > 1)] <- 1


# ggplot heatmap
ggplot(res_DR_WW_merged, aes(x = treatment, y = phylum)) + 
  geom_tile(aes(fill = log2FoldChange))+ # LFC as colour
  geom_text(aes(label = significance), color = "black", size = 4, vjust = 0.8)+ # significance asterisks
  facet_grid(~factor(plant_species, levels=c("P. patens", "M. paleacea", "M. polymorpha")), scales="free_x", space="free_x")+
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

# export plot for Fig. 4B
ggsave("R_plots/bryophytes_heatmap_bac_phylum_level.png", width=1600, height=1000, unit="px")







## 5.2 - ITS (Supp. Fig. 2) ----

alpha <- 0.05


## 5.2.1 - Wilcoxon test ====


# separate to each plant species
filtered_relative_fun_ps_species_l <- metagMisc::phyloseq_sep_variable(filtered_relative_fun_ps, variable = "plant_species")

# remove blank
filtered_relative_fun_ps_species_l <- within(filtered_relative_fun_ps_species_l, rm(blank))

# agglomerate at phylum level
# note that the phyla still have ASV IDs (of a random ASV corresponding to the phylum) as names, instead of phylum names.. we change this after getting the results tables from DESeq2
fun_phyla_rarefied_relative_ps_species_l <- lapply(filtered_relative_fun_ps_species_l, function(x) tax_glom(x, taxrank = "Phylum"))



# lapply goes through the species
phyla_results <- lapply(fun_phyla_rarefied_relative_ps_species_l, function(x){


  # extract otu_table, tax_table and sample_data and make them dataframes
  phyla_table <- as.data.frame(otu_table(x))
  phyla_taxonomy <- as.data.frame(tax_table(x))
  sample_metadata <- as(sample_data(x), "data.frame") # for some reason it needs to be transformed with this function to not give an error with the merge() function later

  # replace ASV names with phyla names
  for (i in 1:nrow(phyla_table)){ # for loop goes through all the row numbers of our results table
    ASV_ID <- rownames(phyla_table)[i] # get the ASV ID of row i
    index <- which(rownames(phyla_taxonomy)==ASV_ID) # get the row index of the tax_table of filtered_fun_ps containing this ASV
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
    abundance_WW <- phyla_avg$abundance[which(phyla_avg$treatment == "20" & phyla_avg$phylum == p)]
    abundance_MDR <- phyla_avg$abundance[which(phyla_avg$treatment == "13" & phyla_avg$phylum == p)]
    results[which(results$phylum == p), "LFC_MDR_WW"] <- log2(abundance_MDR / abundance_WW)

    # p-value by Wilcoxon test
    values_WW <- subset(phyla_table, (treatment == "20" & phylum == p))$abundance
    values_MDR <- subset(phyla_table, (treatment == "13" & phylum == p))$abundance
    results[which(results$phylum == p), "pvalue_MDR_WW"] <- wilcox.test(values_MDR, values_WW)$p.value


    # calculate SDR vs. WW LFC and p-value only for mpol
    if(phyla_table$plant_species[1] %in% c("mpol")){
      # LFC
      abundance_SDR <- phyla_avg$abundance[which(phyla_avg$treatment == "10" & phyla_avg$phylum == p)]
      results[which(results$phylum == p), "LFC_SDR_WW"] <- log2(abundance_SDR / abundance_WW)

      # p-value by Wilcoxon test
      values_SDR <- subset(phyla_table, (treatment == "10" & phylum == p))$abundance
      results[which(results$phylum == p), "pvalue_SDR_WW"] <- wilcox.test(values_SDR, values_WW)$p.value
    }


  }

  return(results)

})



# subset top 12 phyla.  note that they are listed by order decreasing abundance
top_12_phyla <- c("Ascomycota", "Basidiomycota", "Rozellomycota", "Mortierellomycota", "Olpidiomycota", "Aphelidiomycota", "Chytridiomycota", "Glomeromycota", "Zoopagomycota", "Basidiobolomycota", "Monoblepharomycota", "Kickxellomycota" )
phyla_results <- lapply(phyla_results, function(x) x[which(x$phylum %in% top_12_phyla),])



## apply FDR-correction to the p-values
# MDR
phyla_results <- lapply(phyla_results, function(x){
  x$pvalue_MDR_WW <- stats::p.adjust(p = x$pvalue_MDR_WW, method = "fdr")
  return(x)
})
# SDR (only M. polymorpha)
phyla_results$mpol$pvalue_SDR_WW <- stats::p.adjust(p = phyla_results$mpol$pvalue_SDR_WW, method = "fdr")



## 5.2.2 - Merged dataframe for heatmap ====

## MDR (mild drought)
res_MDR_WW_mpol <- data.frame(phylum = phyla_results$mpol$phylum, log2FoldChange = phyla_results$mpol$LFC_MDR_WW, padj = phyla_results$mpol$pvalue_MDR_WW)
res_MDR_WW_mpol$plant_species <- "M. polymorpha"

res_MDR_WW_mpal <- data.frame(phylum = phyla_results$mpal$phylum, log2FoldChange = phyla_results$mpal$LFC_MDR_WW, padj = phyla_results$mpal$pvalue_MDR_WW)
res_MDR_WW_mpal$plant_species <- "M. paleacea"

res_MDR_WW_physco <- data.frame(phylum = phyla_results$physco$phylum, log2FoldChange = phyla_results$physco$LFC_MDR_WW, padj = phyla_results$physco$pvalue_MDR_WW)
res_MDR_WW_physco$plant_species <- "P. patens"

# merge all the results tables
res_MDR_WW <- rbind(res_MDR_WW_mpol, res_MDR_WW_mpal, res_MDR_WW_physco)

# add column "treatment"
res_MDR_WW$treatment <- "mild drought"



## SDR (severe drought)
res_SDR_WW_mpol <- data.frame(phylum = phyla_results$mpol$phylum, log2FoldChange = phyla_results$mpol$LFC_SDR_WW, padj = phyla_results$mpol$pvalue_SDR_WW)
res_SDR_WW_mpol$plant_species <- "M. polymorpha"

# add column "treatment"
res_SDR_WW_mpol$treatment <- "severe drought"



## merge MDR and SDR
res_DR_WW_merged <- rbind(res_MDR_WW, res_SDR_WW_mpol)



# add column "significance" with asterisks for different significance levels (based on adjusted p-value)
# note that this is the adjusted p-value from testing one treatment (MDR or SDR) vs. WW in one plant species. so it is only adjusted for multiple testing with several phyla, not several treatments or species
res_DR_WW_merged$significance <- ""
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.05)] <- "*"
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.01)] <- "**"
res_DR_WW_merged$significance[which(res_DR_WW_merged$padj < 0.001)] <- "***"





## 5.2.3 - Heatmap ====

# set LFC limits to -2 and +2 manually by replacing values under -2 with -2 and over +2 with 2
res_DR_WW_merged$log2FoldChange[which(res_DR_WW_merged$log2FoldChange < -2)] <- -2
res_DR_WW_merged$log2FoldChange[which(res_DR_WW_merged$log2FoldChange > 2)] <- 2


# ggplot heatmap
ggplot(res_DR_WW_merged, aes(x = treatment, y = phylum)) + 
  geom_tile(aes(fill = log2FoldChange))+ # LFC as colour
  geom_text(aes(label = significance), color = "black", size = 4, vjust = 0.8)+ # significance asterisks
  facet_grid(~factor(plant_species, levels=c("P. patens", "M. paleacea", "M. polymorpha")), scales="free_x", space="free_x")+
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

# export plot for Supp. Fig. 2
ggsave("R_plots/bryophytes_heatmap_fun_phylum_level.png", width=1600, height=1000, unit="px")








## 6 - Differential abundance at ASV level: drought vs. well-watered ----


# load filtered physeq object with counts
load("./data_bryophytes/16S/phyloseq_objects/filtered_bac_ps.RData")
load("./data_bryophytes/ITS/phyloseq_objects/filtered_fun_ps.RData")


## 6.1 - mild drought vs. well-watered (all 3 species) ====
# for Supp. Dataset S1 and phylogenetic heatmaps (sections 7.1 and 7.2, for Supp. Fig. 3 and 4)


## 6.1.1 - 16S ####

# separate to each plant species
filtered_bac_ps_species_l <- metagMisc::phyloseq_sep_variable(filtered_bac_ps, variable = "plant_species")

# exclude severe drought
filtered_bac_ps_species_l <- lapply(filtered_bac_ps_species_l, function(x) subset_samples(x, treatment!="10"))

# convert to DESeq object
dds_bac_species_l <- lapply(filtered_bac_ps_species_l, function(x) phyloseq_to_deseq2(x, ~ treatment))

# DESeq with LRT test
dds2_LRT_Bac <- lapply(dds_bac_species_l, function(x)
  DESeq2::DESeq(x,
                test="LRT",
                reduced=~1,
                sfType="poscounts",
                minmu=1e-6,
                minReplicatesForReplace=Inf,
                fitType = "local"))


# Results dataframe generation
res_dds2_LRT_Bac_MDR_WW <- lapply(dds2_LRT_Bac, function(x) DESeq2::results(x, contrast = c("treatment", "13", "20")))  # 13% is numerator and 20% is denominator for LFC




## 6.1.2 - ITS ####

# separate to each plant species
filtered_fun_ps_species_l <- metagMisc::phyloseq_sep_variable(filtered_fun_ps, variable = "plant_species")

# exclude severe drought
filtered_fun_ps_species_l <- lapply(filtered_fun_ps_species_l, function(x) subset_samples(x, treatment!="10"))

# convert to DESeq object
dds_fun_species_l <- lapply(filtered_fun_ps_species_l, function(x) phyloseq_to_deseq2(x, ~ treatment))

# DESeq with LRT test
dds2_LRT_Fun <- lapply(dds_fun_species_l, function(x)
  DESeq2::DESeq(x,
                test="LRT",
                reduced=~1,
                sfType="poscounts",
                minmu=1e-6,
                minReplicatesForReplace=Inf,
                fitType = "local"))


# Results dataframe generation
res_dds2_LRT_Fun_MDR_WW <- lapply(dds2_LRT_Fun, function(x) DESeq2::results(x, contrast = c("treatment", "13", "20")))  # 13% is numerator and 20% is denominator for LFC




## 6.2 - severe drought vs. well-watered (M. polymorpha) ====
# for Supp. Dataset S1


## 6.1.2 - 16S ####

# subset M. polymorpha, severe drought vs. well-watered
filtered_bac_ps_mpol_SDR_WW <- subset_samples(filtered_bac_ps, plant_species == "mpol" & treatment %in% c("20", "10"))

# convert to DESeq object
dds_bac_species_mpol_SDR_WW <- phyloseq_to_deseq2(filtered_bac_ps_mpol_SDR_WW, ~ treatment)

# DESeq with LRT test
dds2_LRT_Bac_mpol_SDR_WW <- DESeq2::DESeq(dds_bac_species_mpol_SDR_WW,
                test="LRT",
                reduced=~1,
                sfType="poscounts",
                minmu=1e-6,
                minReplicatesForReplace=Inf,
                fitType = "local")


# Results dataframe generation
res_dds2_LRT_Bac_mpol_SDR_WW <- DESeq2::results(dds2_LRT_Bac_mpol_SDR_WW, contrast = c("treatment", "10", "20"))  # 10% is numerator and 20% is denominator for LFC




## 6.2.2 - ITS ####

# subset M. polymorpha, severe drought vs. well-watered
filtered_fun_ps_mpol_SDR_WW <- subset_samples(filtered_fun_ps, plant_species == "mpol" & treatment %in% c("20", "10"))

# convert to DESeq object
dds_fun_species_mpol_SDR_WW <- phyloseq_to_deseq2(filtered_fun_ps_mpol_SDR_WW, ~ treatment)

# DESeq with LRT test
dds2_LRT_Fun_mpol_SDR_WW <- DESeq2::DESeq(dds_fun_species_mpol_SDR_WW,
                                          test="LRT",
                                          reduced=~1,
                                          sfType="poscounts",
                                          minmu=1e-6,
                                          minReplicatesForReplace=Inf,
                                          fitType = "local")


# Results dataframe generation
res_dds2_LRT_Fun_mpol_SDR_WW <- DESeq2::results(dds2_LRT_Fun_mpol_SDR_WW, contrast = c("treatment", "10", "20"))  # 10% is numerator and 20% is denominator for LFC




## 6.3 - Format DESeq results for Supp. Dataset S1 ====

## 6.3.1 - extract results tables and convert to data frames ####

# 16S
res_df_bac_mpol_MDR_WW <- as.data.frame(res_dds2_LRT_Bac_MDR_WW$mpol)
res_df_bac_mpal_MDR_WW <- as.data.frame(res_dds2_LRT_Bac_MDR_WW$mpal)
res_df_bac_physco_MDR_WW <- as.data.frame(res_dds2_LRT_Bac_MDR_WW$physco)
res_df_bac_mpol_SDR_WW <- as.data.frame(res_dds2_LRT_Bac_mpol_SDR_WW)

# ITS
res_df_fun_mpol_MDR_WW <- as.data.frame(res_dds2_LRT_Fun_MDR_WW$mpol)
res_df_fun_mpal_MDR_WW <- as.data.frame(res_dds2_LRT_Fun_MDR_WW$mpal)
res_df_fun_physco_MDR_WW <- as.data.frame(res_dds2_LRT_Fun_MDR_WW$physco)
res_df_fun_mpol_SDR_WW <- as.data.frame(res_dds2_LRT_Fun_mpol_SDR_WW)



## 6.3.2 - prepare taxonomy ####

# extract taxonomy from phyloseq objects and change to data frame
bac_tax_df <- as.data.frame(tax_table(filtered_bac_ps))
fun_tax_df <- as.data.frame(tax_table(filtered_fun_ps))

# move ASV IDs from rownames to new column (for merging)
bac_tax_df <- cbind(rownames(bac_tax_df), bac_tax_df)
colnames(bac_tax_df)[1] <- "ASV ID"
fun_tax_df <- cbind(rownames(fun_tax_df), fun_tax_df)
colnames(fun_tax_df)[1] <- "ASV ID"

# remove prefixes "k__", "p__", etc.
bac_tax_df$Kingdom <- gsub("k__", "", bac_tax_df$Kingdom)
bac_tax_df$Phylum <- gsub("p__", "", bac_tax_df$Phylum)
bac_tax_df$Class <- gsub("c__", "", bac_tax_df$Class)
bac_tax_df$Order <- gsub("o__", "", bac_tax_df$Order)
bac_tax_df$Family <- gsub("f__", "", bac_tax_df$Family)
bac_tax_df$Genus <- gsub("g__", "", bac_tax_df$Genus)
bac_tax_df$Species <- gsub("s__", "", bac_tax_df$Species)

fun_tax_df$Kingdom <- gsub("k__", "", fun_tax_df$Kingdom)
fun_tax_df$Phylum <- gsub("p__", "", fun_tax_df$Phylum)
fun_tax_df$Class <- gsub("c__", "", fun_tax_df$Class)
fun_tax_df$Order <- gsub("o__", "", fun_tax_df$Order)
fun_tax_df$Family <- gsub("f__", "", fun_tax_df$Family)
fun_tax_df$Genus <- gsub("g__", "", fun_tax_df$Genus)
fun_tax_df$Species <- gsub("s__", "", fun_tax_df$Species)



## 6.3.3 - prepare DESeq results dataframes for merging ####

# DESeq results dataframes: move ASV IDs from rownames to new column (for merging)
res_df_bac_mpol_MDR_WW <- cbind(rownames(res_df_bac_mpol_MDR_WW), res_df_bac_mpol_MDR_WW)
colnames(res_df_bac_mpol_MDR_WW)[1] <- "ASV ID"
res_df_bac_mpal_MDR_WW <- cbind(rownames(res_df_bac_mpal_MDR_WW), res_df_bac_mpal_MDR_WW)
colnames(res_df_bac_mpal_MDR_WW)[1] <- "ASV ID"
res_df_bac_physco_MDR_WW <- cbind(rownames(res_df_bac_physco_MDR_WW), res_df_bac_physco_MDR_WW)
colnames(res_df_bac_physco_MDR_WW)[1] <- "ASV ID"
res_df_bac_mpol_SDR_WW <- cbind(rownames(res_df_bac_mpol_SDR_WW), res_df_bac_mpol_SDR_WW)
colnames(res_df_bac_mpol_SDR_WW)[1] <- "ASV ID"

res_df_fun_mpol_MDR_WW <- cbind(rownames(res_df_fun_mpol_MDR_WW), res_df_fun_mpol_MDR_WW)
colnames(res_df_fun_mpol_MDR_WW)[1] <- "ASV ID"
res_df_fun_mpal_MDR_WW <- cbind(rownames(res_df_fun_mpal_MDR_WW), res_df_fun_mpal_MDR_WW)
colnames(res_df_fun_mpal_MDR_WW)[1] <- "ASV ID"
res_df_fun_physco_MDR_WW <- cbind(rownames(res_df_fun_physco_MDR_WW), res_df_fun_physco_MDR_WW)
colnames(res_df_fun_physco_MDR_WW)[1] <- "ASV ID"
res_df_fun_mpol_SDR_WW <- cbind(rownames(res_df_fun_mpol_SDR_WW), res_df_fun_mpol_SDR_WW)
colnames(res_df_fun_mpol_SDR_WW)[1] <- "ASV ID"



## order by ascending ASV ID (16S)

# put ASV ID number in a new column for sorting
res_df_bac_mpol_MDR_WW$ASV_ID_number <- as.numeric(gsub("ASV_", "", res_df_bac_mpol_MDR_WW$`ASV ID`))
res_df_bac_mpal_MDR_WW$ASV_ID_number <- as.numeric(gsub("ASV_", "", res_df_bac_mpal_MDR_WW$`ASV ID`))
res_df_bac_physco_MDR_WW$ASV_ID_number <- as.numeric(gsub("ASV_", "", res_df_bac_physco_MDR_WW$`ASV ID`))
res_df_bac_mpol_SDR_WW$ASV_ID_number <- as.numeric(gsub("ASV_", "", res_df_bac_mpol_SDR_WW$`ASV ID`))

# sort by ASV ID number
res_df_bac_mpol_MDR_WW <- res_df_bac_mpol_MDR_WW[order(res_df_bac_mpol_MDR_WW$ASV_ID_number, decreasing=FALSE),]
res_df_bac_mpal_MDR_WW <- res_df_bac_mpal_MDR_WW[order(res_df_bac_mpal_MDR_WW$ASV_ID_number, decreasing=FALSE),]
res_df_bac_physco_MDR_WW <- res_df_bac_physco_MDR_WW[order(res_df_bac_physco_MDR_WW$ASV_ID_number, decreasing=FALSE),]
res_df_bac_mpol_SDR_WW <- res_df_bac_mpol_SDR_WW[order(res_df_bac_mpol_SDR_WW$ASV_ID_number, decreasing=FALSE),]

# drop column "ASV_ID_number"
res_df_bac_mpol_MDR_WW <- res_df_bac_mpol_MDR_WW[,which(colnames(res_df_bac_mpol_MDR_WW) != "ASV_ID_number")]
res_df_bac_mpal_MDR_WW <- res_df_bac_mpal_MDR_WW[,which(colnames(res_df_bac_mpal_MDR_WW) != "ASV_ID_number")]
res_df_bac_physco_MDR_WW <- res_df_bac_physco_MDR_WW[,which(colnames(res_df_bac_physco_MDR_WW) != "ASV_ID_number")]
res_df_bac_mpol_SDR_WW <- res_df_bac_mpol_SDR_WW[,which(colnames(res_df_bac_mpol_SDR_WW) != "ASV_ID_number")]


## order by ascending ASV ID (ITS)

# put ASV ID number in a new column for sorting
res_df_fun_mpol_MDR_WW$ASV_ID_number <- as.numeric(gsub("fASV_", "", res_df_fun_mpol_MDR_WW$`ASV ID`))
res_df_fun_mpal_MDR_WW$ASV_ID_number <- as.numeric(gsub("fASV_", "", res_df_fun_mpal_MDR_WW$`ASV ID`))
res_df_fun_physco_MDR_WW$ASV_ID_number <- as.numeric(gsub("fASV_", "", res_df_fun_physco_MDR_WW$`ASV ID`))
res_df_fun_mpol_SDR_WW$ASV_ID_number <- as.numeric(gsub("fASV_", "", res_df_fun_mpol_SDR_WW$`ASV ID`))

# sort by ASV ID number
res_df_fun_mpol_MDR_WW <- res_df_fun_mpol_MDR_WW[order(res_df_fun_mpol_MDR_WW$ASV_ID_number, decreasing=FALSE),]
res_df_fun_mpal_MDR_WW <- res_df_fun_mpal_MDR_WW[order(res_df_fun_mpal_MDR_WW$ASV_ID_number, decreasing=FALSE),]
res_df_fun_physco_MDR_WW <- res_df_fun_physco_MDR_WW[order(res_df_fun_physco_MDR_WW$ASV_ID_number, decreasing=FALSE),]
res_df_fun_mpol_SDR_WW <- res_df_fun_mpol_SDR_WW[order(res_df_fun_mpol_SDR_WW$ASV_ID_number, decreasing=FALSE),]

# drop column "ASV_ID_number"
res_df_fun_mpol_MDR_WW <- res_df_fun_mpol_MDR_WW[,which(colnames(res_df_fun_mpol_MDR_WW) != "ASV_ID_number")]
res_df_fun_mpal_MDR_WW <- res_df_fun_mpal_MDR_WW[,which(colnames(res_df_fun_mpal_MDR_WW) != "ASV_ID_number")]
res_df_fun_physco_MDR_WW <- res_df_fun_physco_MDR_WW[,which(colnames(res_df_fun_physco_MDR_WW) != "ASV_ID_number")]
res_df_fun_mpol_SDR_WW <- res_df_fun_mpol_SDR_WW[,which(colnames(res_df_fun_mpol_SDR_WW) != "ASV_ID_number")]




## 6.3.4 - add taxonomy to DESeq results dataframes ####

# 16S
res_df_bac_mpol_MDR_WW <- left_join(res_df_bac_mpol_MDR_WW, bac_tax_df, by = "ASV ID")
res_df_bac_mpal_MDR_WW <-  left_join(res_df_bac_mpal_MDR_WW, bac_tax_df, by = "ASV ID")
res_df_bac_physco_MDR_WW <-  left_join(res_df_bac_physco_MDR_WW, bac_tax_df, by = "ASV ID")
res_df_bac_mpol_SDR_WW <-  left_join(res_df_bac_mpol_SDR_WW, bac_tax_df, by = "ASV ID")

# ITS
res_df_fun_mpol_MDR_WW <- left_join(res_df_fun_mpol_MDR_WW, fun_tax_df, by = "ASV ID")
res_df_fun_mpal_MDR_WW <-  left_join(res_df_fun_mpal_MDR_WW, fun_tax_df, by = "ASV ID")
res_df_fun_physco_MDR_WW <-  left_join(res_df_fun_physco_MDR_WW, fun_tax_df, by = "ASV ID")
res_df_fun_mpol_SDR_WW <-  left_join(res_df_fun_mpol_SDR_WW, fun_tax_df, by = "ASV ID")





## 6.3.5 - assemble results tables to an Excel sheet for Supp. Dataset ####

# make a first sheet "abbreviations" explaining what "WW", "MDR" and "SDR" mean. also "16S" and "ITS"
abbreviations_df <- data.frame(c("", "WW: well-watered", "MDR: mild drought", "SDR: severe drought", "16S: 16S amplicon sequencing (bacteria)", "ITS: ITS amplicon sequencing (fungi)"))
colnames(abbreviations_df) <- "List of abbreviations used in this dataset"


# each pairwise comparison becomes one tab in the Excel sheet
# note: different order than the code above. first physco MDR, then mpal MDR, then mpol MDR, then mpol SDR

names <- list("Abbreviations" = abbreviations_df,
              "16S, P. patens, MDR vs. WW" = res_df_bac_physco_MDR_WW,
              "16S, M. paleacea, MDR vs. WW" = res_df_bac_mpal_MDR_WW,
              "16S, M. polymorpha, MDR vs. WW" = res_df_bac_mpol_MDR_WW,
              "16S, M. polymorpha, SDR vs. WW" = res_df_bac_mpol_SDR_WW,
              "ITS, P. patens, MDR vs. WW" = res_df_fun_physco_MDR_WW,
              "ITS, M. paleacea, MDR vs. WW" = res_df_fun_mpal_MDR_WW,
              "ITS, M. polymorpha, MDR vs. WW" = res_df_fun_mpol_MDR_WW,
              "ITS, M. polymorpha, SDR vs. WW" = res_df_fun_mpol_SDR_WW)

openxlsx::write.xlsx(names, file = "supp_dataset_1.xlsx")






## 7 - Phylogenetic heatmaps of top 200 bacterial ASVs ----
# Supp. Fig. 3 & 5

# load filtered physeq objects with relative abundances
load("./data_bryophytes/16S/phyloseq_objects/filtered_relative_bac_ps.RData")
load("./data_bryophytes/ITS/phyloseq_objects/filtered_relative_fun_ps.RData")


## 7.1 - Subset top N ASVs ====

# how many ASVs?
bac_top_N <- 200

# list top N ASVs
bac_top_N_ASVs <- names(sort(taxa_sums(filtered_relative_bac_ps), decreasing=TRUE)[1:bac_top_N])

# filter top N ASVs
bac_ps_top_N <- prune_taxa(bac_top_N_ASVs, filtered_relative_bac_ps)



## 7.2 - Prepare phylogeny ====

# load bacterial phylogeny
bac_phylogeny <- read.tree("./data_bryophytes/16S/tree_Bac.nwk")
# plot.phylo(bac_phylogeny) # plot tree

# remove extra apostrophes from tip labels
bac_phylogeny$tip.label <- gsub("'", "", bac_phylogeny$tip.label)

# prune phylogeny to ASVs present in table
bac_phylogeny_top_N <- keep.tip(bac_phylogeny, bac_top_N_ASVs)
plot.phylo(bac_phylogeny_top_N) # plot tree

# align tree tips to make dendogram later
bac_phylogeny_top_N <- chronos(bac_phylogeny_top_N)

# turn the tree to a hclust object
# bac_phylogeny_top_N <- RRphylo::fix.poly(bac_phylogeny_top_N, type="resolve") # resolve polytomies to non-zero branch lengths, as the function as.hclust() only takes binary trees
bac_phylogeny_top_N <- as.hclust(bac_phylogeny_top_N) # convert to hclust object
plot(bac_phylogeny_top_N)



## prepare Genus names for rownames

# first let's load the tax table (output from qiime2-sklearn)
raw_bac_taxtab <- read.table(file = "C:/bryophytes_rhizosphere/data_bryophytes/16S/taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
raw_bac_taxtab <- separate(data = raw_bac_taxtab, col = Taxon, into =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

# get the genus names of the ASVs of our data frames (via rownames)
genus_labels_bac_top_N <- raw_bac_taxtab[bac_phylogeny_top_N$labels, "Genus"]

# remove "g__" before genus names
genus_labels_bac_top_N <- gsub(" g__", "", genus_labels_bac_top_N)

# change NA in genus labels to "NA" as text
genus_labels_bac_top_N[is.na(genus_labels_bac_top_N)] <- "NA"

# find which labels should be italicized: all genus labels except "uncultured" and "NA"
label_numbers_to_italicize <- which(!(genus_labels_bac_top_N %in% c("uncultured", "NA")))

# copy genus labels list to italicized genus labels list
genus_labels_bac_top_N_italicized <- genus_labels_bac_top_N

# italicize genus labels except "uncultured" and "NA"
for(i in label_numbers_to_italicize){
  genus_labels_bac_top_N_italicized[i] <- as.expression(bquote(italic(.(genus_labels_bac_top_N[i]))))
}




## 7.3 - Heatmap of relative abundance in each species/treatment ====

# use psmelt to combine the ASV table with metadata and taxonomy
bac_top_N_melted <- psmelt(bac_ps_top_N)

# calculate average relative abundance per species/treatment
bac_top_N_species_treatment <- bac_top_N_melted %>% 
  group_by(OTU, plant_species, treatment) %>%
  summarise(Abundance = mean(Abundance))



# make concatenated species/treatment column
bac_top_N_species_treatment$species_treatment <- paste0(bac_top_N_species_treatment$plant_species, "_", bac_top_N_species_treatment$treatment)

# remove columns plant_species and treatment
bac_top_N_species_treatment <- bac_top_N_species_treatment[, which(colnames(bac_top_N_species_treatment) %in% c("OTU", "species_treatment", "Abundance"))]

# pivot to wide format: ASVs as rows, species/treatment as columns
bac_top_N_species_treatment <- bac_top_N_species_treatment %>% pivot_wider(names_from='species_treatment',
                                                                           values_from='Abundance')
# change to dataframe
bac_top_N_species_treatment <- as.data.frame(bac_top_N_species_treatment)

# move ASV names to rownames
rownames(bac_top_N_species_treatment) <- bac_top_N_species_treatment$OTU

# remove column "OTU"
bac_top_N_species_treatment <- bac_top_N_species_treatment[,which(colnames(bac_top_N_species_treatment) != "OTU")]

# order and rename columns for heatmap
columns_order <- c("mpal_20", "mpal_13", "mpol_20", "mpol_13", "mpol_10", "physco_20", "physco_13")
bac_top_N_species_treatment <- bac_top_N_species_treatment[, columns_order]
colnames(bac_top_N_species_treatment) <- c("M. paleacea - well-watered", "M. paleacea - mild drought", "M. polymorpha - well-watered", "M. polymorpha - mild drought", "M. polymorpha - severe drought", "P. patens - well-watered", "P. patens - mild drought")




# order rows for heatmap
# rows should be in the same order as the ASV phylogeny (hclust object)
bac_top_N_species_treatment <- bac_top_N_species_treatment[bac_phylogeny_top_N$labels,]



# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 0.05
myBreaks <- c(seq(min1, max1, length=paletteLength))



# make heatmap
pheatmap(bac_top_N_species_treatment, cluster_cols=FALSE, cluster_rows=bac_phylogeny_top_N, labels_row=genus_labels_bac_top_N, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")



## 7.4 - Heatmap of normalized relative abundance in each species/treatment (Supp. Fig. 5) ====

# normalized by ASV-average, to see the host species and treatment preference of each ASV

# normalize rel. abundance by average of each ASV across plant species
bac_top_N_species_treatment_normalized <- bac_top_N_species_treatment*ncol(bac_top_N_species_treatment)/rowSums(bac_top_N_species_treatment)



## define colours and breaks for heatmaps - for normalized relative abundance

# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 5
myBreaks <- c(seq(min1, max1, length=paletteLength))


# re-order species/treatments for consistency (first P. patens, then M. paleacea, then M. polymorpha)
bac_top_N_species_treatment_normalized <- bac_top_N_species_treatment_normalized[, c(6, 7, 1:5)]

# make column labels with italicized species names
hm_column_labels_italicized <- expression(italic("P. patens") ~ "- well-watered", italic("P. patens") ~ "- mild drought", italic("M. paleacea") ~ "- well-watered", italic("M. paleacea") ~ "- mild drought", italic("M. polymorpha") ~ "- well-watered", italic("M. polymorpha") ~ "- mild drought", italic("M. polymorpha") ~ "- severe drought")

# make heatmap
hm <- pheatmap(bac_top_N_species_treatment_normalized, cluster_cols=FALSE, cluster_rows=bac_phylogeny_top_N, labels_row=genus_labels_bac_top_N_italicized, labels_col=hm_column_labels_italicized, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")

# export plot for Supp. Fig. 5
ggsave(plot=hm, filename="R_plots/bac_top_200_host_treatment_preferences.png", width=2000, height=8000, unit="px")




## 7.5 - Heatmap of LFC mild drought vs. well-watered in each species (Supp. Fig. 3) ====


## extract DESeq2 results for top N ASVs
res_dds2_LRT_Bac_MDR_WW_top_N <- lapply(res_dds2_LRT_Bac_MDR_WW, function(x) x[rownames(x) %in% bac_top_N_ASVs,])



## get LFCs from DESeq output

# extract log2FoldChange MDR vs. WW from each species (DESeq2 output), to a list of data frames
LFC_top_N_MDR_WW <- lapply(res_dds2_LRT_Bac_MDR_WW_top_N, function(x) {
  x <- as(x, "data.frame")
  x$ASV <- rownames(x)
  x <- x[c("ASV", "log2FoldChange")]
  x <- subset(x, ASV %in% bac_top_N_ASVs)
  return(x)
})

# rename column "log2FoldChange" to the species name (hard-coded, don't know how to access list element name in lapply)
LFC_top_N_MDR_WW$mpol <- LFC_top_N_MDR_WW$mpol %>% dplyr::rename(mpol = log2FoldChange)
LFC_top_N_MDR_WW$mpal <- LFC_top_N_MDR_WW$mpal %>% dplyr::rename(mpal = log2FoldChange)
LFC_top_N_MDR_WW$physco <- LFC_top_N_MDR_WW$physco %>% dplyr::rename(physco = log2FoldChange)

# merge the dataframes to one
LFC_top_N_MDR_WW_merged <- merge(LFC_top_N_MDR_WW$mpol, LFC_top_N_MDR_WW$mpal, by = "ASV", all = TRUE)
LFC_top_N_MDR_WW_merged <- merge(LFC_top_N_MDR_WW_merged, LFC_top_N_MDR_WW$physco, by = "ASV", all = TRUE)
LFC_top_N_MDR_WW_merged
rownames(LFC_top_N_MDR_WW_merged) <- LFC_top_N_MDR_WW_merged$ASV # put ASV as rownames
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[, which(colnames(LFC_top_N_MDR_WW_merged) != "ASV")] # remove column ASV




## get unadjusted p-values from DESeq output and FDR-correct them for N tests (the number of ASVs that are shown in the heatmap)

# extract unadjusted p-values for MDR vs. WW from each species (DESeq2 output), to a list of data frames
pvalue_top_N_MDR_WW <- lapply(res_dds2_LRT_Bac_MDR_WW_top_N, function(x) {
  x <- as(x, "data.frame")
  x$ASV <- rownames(x)
  x <- x[c("ASV", "pvalue")]
  x <- subset(x, ASV %in% bac_top_N_ASVs)
  return(x)
})

# rename column "pvalue" to the species name (hard-coded, don't know how to access list element name in lapply)
pvalue_top_N_MDR_WW$mpol <- pvalue_top_N_MDR_WW$mpol %>% dplyr::rename(mpol = pvalue)
pvalue_top_N_MDR_WW$mpal <- pvalue_top_N_MDR_WW$mpal %>% dplyr::rename(mpal = pvalue)
pvalue_top_N_MDR_WW$physco <- pvalue_top_N_MDR_WW$physco %>% dplyr::rename(physco = pvalue)

# merge the dataframes to one
pvalue_top_N_MDR_WW_merged <- merge(pvalue_top_N_MDR_WW$mpol, pvalue_top_N_MDR_WW$mpal, by = "ASV", all = TRUE)
pvalue_top_N_MDR_WW_merged <- merge(pvalue_top_N_MDR_WW_merged, pvalue_top_N_MDR_WW$physco, by = "ASV", all = TRUE)
pvalue_top_N_MDR_WW_merged
rownames(pvalue_top_N_MDR_WW_merged) <- pvalue_top_N_MDR_WW_merged$ASV # put ASV as rownames
pvalue_top_N_MDR_WW_merged <- pvalue_top_N_MDR_WW_merged[, which(colnames(pvalue_top_N_MDR_WW_merged) != "ASV")] # remove column ASV




## FDR-correct p-values for N tests (the number of ASVs that are shown in the heatmap)
pvalue_top_N_MDR_WW_merged$mpol <- stats::p.adjust(p = pvalue_top_N_MDR_WW_merged$mpol, method = "fdr") # "fdr" = Benjamini-Hochberg
pvalue_top_N_MDR_WW_merged$mpal <- stats::p.adjust(p = pvalue_top_N_MDR_WW_merged$mpal, method = "fdr") # "fdr" = Benjamini-Hochberg
pvalue_top_N_MDR_WW_merged$physco <- stats::p.adjust(p = pvalue_top_N_MDR_WW_merged$physco, method = "fdr") # "fdr" = Benjamini-Hochberg




## turn p-values into significance asterisks

# now transform p-values into asterisks
asterisks_top_N_MDR_WW_merged <- pvalue_top_N_MDR_WW_merged
asterisks_top_N_MDR_WW_merged[pvalue_top_N_MDR_WW_merged > 0.05] <- ""
asterisks_top_N_MDR_WW_merged[pvalue_top_N_MDR_WW_merged <= 0.05] <- "*"
asterisks_top_N_MDR_WW_merged[is.na(pvalue_top_N_MDR_WW_merged)] <- ""






# order and rename columns for heatmap
columns_order <- c("mpal", "mpol", "physco")
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[, columns_order]
asterisks_top_N_MDR_WW_merged <- asterisks_top_N_MDR_WW_merged[, columns_order]
colnames(LFC_top_N_MDR_WW_merged) <- c("M. paleacea", "M. polymorpha", "P. patens")
colnames(asterisks_top_N_MDR_WW_merged) <- c("M. paleacea", "M. polymorpha", "P. patens")


# order rows for heatmap
# rows should be in the same order as the ASV phylogeny (hclust object)
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[bac_phylogeny_top_N$labels,]
asterisks_top_N_MDR_WW_merged <- asterisks_top_N_MDR_WW_merged[bac_phylogeny_top_N$labels,]




# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- -3
mid1 <- 0
max1 <- 3
myBreaks <- c(seq(min1, mid1, length.out=ceiling(paletteLength/2) + 1),
              seq(mid1+0.1, max1, length.out=floor(paletteLength/2)))



# re-order species for consistency (first P. patens, then M. paleacea, then M. polymorpha)
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[, c(3, 1, 2)]
asterisks_top_N_MDR_WW_merged <- asterisks_top_N_MDR_WW_merged[, c(3, 1, 2)]
colnames(LFC_top_N_MDR_WW_merged) == colnames(asterisks_top_N_MDR_WW_merged)

# make column labels with italicized species names
hm_column_labels_italicized <- expression(italic("P. patens"), italic("M. paleacea"), italic("M. polymorpha"))

# make heatmap
hm <- pheatmap(LFC_top_N_MDR_WW_merged, cluster_cols=FALSE, cluster_rows=bac_phylogeny_top_N, labels_row=genus_labels_bac_top_N_italicized, labels_col=hm_column_labels_italicized, display_numbers=asterisks_top_N_MDR_WW_merged, color = myColor, breaks = myBreaks, na_col = "grey")

# export plot for Supp. Fig. 3
ggsave(plot=hm, filename="R_plots/bac_top_200_LFC_mild_drought.png", width=1800, height=8000, unit="px")




## 8 - Phylogenetic heatmaps of top 100 fungal ASVs ----
# Supp. Fig. 4 & 6


## 8.1 - Subset top N ASVs ====

# how many ASVs?
fun_top_N <- 100

# list top N ASVs
fun_top_N_ASVs <- names(sort(taxa_sums(filtered_relative_fun_ps), decreasing=TRUE)[1:fun_top_N])

# filter top N ASVs
fun_ps_top_N <- prune_taxa(fun_top_N_ASVs, filtered_relative_fun_ps)



## 8.2 - Prepare phylogeny ====

# load fungal phylogeny
fun_phylogeny <- read.tree("./data_bryophytes/ITS/tree_Fun.nwk")
# plot.phylo(fun_phylogeny) # plot tree

# remove extra apostrophes from tip labels
fun_phylogeny$tip.label <- gsub("'", "", fun_phylogeny$tip.label)

# prune phylogeny to ASVs present in table
fun_phylogeny_top_N <- keep.tip(fun_phylogeny, fun_top_N_ASVs)
plot.phylo(fun_phylogeny_top_N) # plot tree

# align tree tips to make dendogram later
fun_phylogeny_top_N <- chronos(fun_phylogeny_top_N)

# turn the tree to a hclust object
fun_phylogeny_top_N <- RRphylo::fix.poly(fun_phylogeny_top_N, type="resolve") # resolve polytomies to non-zero branch lengths, as the function as.hclust() only takes binary trees
fun_phylogeny_top_N <- as.hclust(fun_phylogeny_top_N) # convert to hclust object
plot(fun_phylogeny_top_N)



## prepare Genus names for rownames

# first let's load the tax table (output from qiime2-sklearn)
raw_fun_taxtab<-read.table(file = "C:/bryophytes_rhizosphere/data_bryophytes/ITS/taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
raw_fun_taxtab<- separate(data = raw_fun_taxtab, col = Taxon, into =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
rownames(raw_fun_taxtab) <- paste0("f", rownames(raw_fun_taxtab))

# get the genus names of the ASVs of our data frames (via rownames)
genus_labels_fun_top_N <- raw_fun_taxtab[fun_phylogeny_top_N$labels, "Genus"]

# remove "g__" before genus names
genus_labels_fun_top_N <- gsub("g__", "", genus_labels_fun_top_N)

# change NA in genus labels to "NA" as text
genus_labels_fun_top_N[is.na(genus_labels_fun_top_N)] <- "NA"

# find which labels should be italicized: all genus labels except "unidentified" and "NA"
label_numbers_to_italicize <- which(!(genus_labels_fun_top_N %in% c("unidentified", "NA")))

# copy genus labels list to italicized genus labels list
genus_labels_fun_top_N_italicized <- genus_labels_fun_top_N

# italicize genus labels except "unidentified" and "NA"
for(i in label_numbers_to_italicize){
  genus_labels_fun_top_N_italicized[i] <- as.expression(bquote(italic(.(genus_labels_fun_top_N[i]))))
}





## 8.3 - Heatmap of relative abundance in each species/treatment ====

# use psmelt to combine the ASV table with metadata and taxonomy
fun_top_N_melted <- psmelt(fun_ps_top_N)

# calculate average relative abundance per species/treatment
fun_top_N_species_treatment <- fun_top_N_melted %>% 
  group_by(OTU, plant_species, treatment) %>%
  summarise(Abundance = mean(Abundance))



# make concatenated species/treatment column
fun_top_N_species_treatment$species_treatment <- paste0(fun_top_N_species_treatment$plant_species, "_", fun_top_N_species_treatment$treatment)

# remove columns plant_species and treatment
fun_top_N_species_treatment <- fun_top_N_species_treatment[, which(colnames(fun_top_N_species_treatment) %in% c("OTU", "species_treatment", "Abundance"))]

# pivot to wide format: ASVs as rows, species/treatment as columns
fun_top_N_species_treatment <- fun_top_N_species_treatment %>% pivot_wider(names_from='species_treatment',
                                                                           values_from='Abundance')
# change to dataframe
fun_top_N_species_treatment <- as.data.frame(fun_top_N_species_treatment)

# move ASV names to rownames
rownames(fun_top_N_species_treatment) <- fun_top_N_species_treatment$OTU

# remove column "OTU"
fun_top_N_species_treatment <- fun_top_N_species_treatment[,which(colnames(fun_top_N_species_treatment) != "OTU")]

# order and rename columns for heatmap
columns_order <- c("mpal_20", "mpal_13", "mpol_20", "mpol_13", "mpol_10", "physco_20", "physco_13")
fun_top_N_species_treatment <- fun_top_N_species_treatment[, columns_order]
colnames(fun_top_N_species_treatment) <- c("M. paleacea - well-watered", "M. paleacea - mild drought", "M. polymorpha - well-watered", "M. polymorpha - mild drought", "M. polymorpha - severe drought", "P. patens - well-watered", "P. patens - mild drought")



# order rows for heatmap
# rows should be in the same order as the ASV phylogeny (hclust object)
fun_top_N_species_treatment <- fun_top_N_species_treatment[fun_phylogeny_top_N$labels,]



# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 0.05
myBreaks <- c(seq(min1, max1, length=paletteLength))




# make heatmap
pheatmap(fun_top_N_species_treatment, cluster_cols=FALSE, cluster_rows=fun_phylogeny_top_N, labels_row=genus_labels_fun_top_N, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")





## 8.4 - Heatmap of normalized relative abundance in each species/treatment (Supp. Fig. 6) ====

# normalized by ASV-average, to see the host species and treatment preference of each ASV

# normalize rel. abundance by average of each ASV across plant species
fun_top_N_species_treatment_normalized <- fun_top_N_species_treatment*ncol(fun_top_N_species_treatment)/rowSums(fun_top_N_species_treatment)



## define colours and breaks for heatmaps - for normalized relative abundance

# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 5
myBreaks <- c(seq(min1, max1, length=paletteLength))



# re-order species/treatments for consistency (first P. patens, then M. paleacea, then M. polymorpha)
fun_top_N_species_treatment_normalized <- fun_top_N_species_treatment_normalized[, c(6, 7, 1:5)]

# make column labels with italicized species names
hm_column_labels_italicized <- expression(italic("P. patens") ~ "- well-watered", italic("P. patens") ~ "- mild drought", italic("M. paleacea") ~ "- well-watered", italic("M. paleacea") ~ "- mild drought", italic("M. polymorpha") ~ "- well-watered", italic("M. polymorpha") ~ "- mild drought", italic("M. polymorpha") ~ "- severe drought")

# make heatmap
hm <- pheatmap(fun_top_N_species_treatment_normalized, cluster_cols=FALSE, cluster_rows=fun_phylogeny_top_N, labels_row=genus_labels_fun_top_N_italicized, labels_col=hm_column_labels_italicized, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")

# export plot for Supp. Fig. 6
ggsave(plot=hm, filename="R_plots/fun_top_100_host_treatment_preferences.png", width=1300, height=5000, unit="px")





## 8.5 - Heatmap of LFC mild drought vs. well-watered in each species (Supp. Fig. 4) ====


## extract DESeq2 results for top N ASVs
res_dds2_LRT_Fun_MDR_WW_top_N <- lapply(res_dds2_LRT_Fun_MDR_WW, function(x) x[rownames(x) %in% fun_top_N_ASVs,])



## get LFCs from DESeq output

# extract log2FoldChange MDR vs. WW from each species (DESeq2 output), to a list of data frames
LFC_top_N_MDR_WW <- lapply(res_dds2_LRT_Fun_MDR_WW_top_N, function(x) {
  x <- as(x, "data.frame")
  x$ASV <- rownames(x)
  x <- x[c("ASV", "log2FoldChange")]
  x <- subset(x, ASV %in% fun_top_N_ASVs)
  return(x)
})

# rename column "log2FoldChange" to the species name (hard-coded, don't know how to access list element name in lapply)
LFC_top_N_MDR_WW$mpol <- LFC_top_N_MDR_WW$mpol %>% dplyr::rename(mpol = log2FoldChange)
LFC_top_N_MDR_WW$mpal <- LFC_top_N_MDR_WW$mpal %>% dplyr::rename(mpal = log2FoldChange)
LFC_top_N_MDR_WW$physco <- LFC_top_N_MDR_WW$physco %>% dplyr::rename(physco = log2FoldChange)

# merge the dataframes to one
LFC_top_N_MDR_WW_merged <- merge(LFC_top_N_MDR_WW$mpol, LFC_top_N_MDR_WW$mpal, by = "ASV", all = TRUE)
LFC_top_N_MDR_WW_merged <- merge(LFC_top_N_MDR_WW_merged, LFC_top_N_MDR_WW$physco, by = "ASV", all = TRUE)
LFC_top_N_MDR_WW_merged
rownames(LFC_top_N_MDR_WW_merged) <- LFC_top_N_MDR_WW_merged$ASV # put ASV as rownames
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[, which(colnames(LFC_top_N_MDR_WW_merged) != "ASV")] # remove column ASV




## get unadjusted p-values from DESeq output and FDR-correct them for N tests (the number of ASVs that are shown in the heatmap)

# extract unadjusted p-values for MDR vs. WW from each species (DESeq2 output), to a list of data frames
pvalue_top_N_MDR_WW <- lapply(res_dds2_LRT_Fun_MDR_WW_top_N, function(x) {
  x <- as(x, "data.frame")
  x$ASV <- rownames(x)
  x <- x[c("ASV", "pvalue")]
  x <- subset(x, ASV %in% fun_top_N_ASVs)
  return(x)
})

# rename column "pvalue" to the species name (hard-coded, don't know how to access list element name in lapply)
pvalue_top_N_MDR_WW$mpol <- pvalue_top_N_MDR_WW$mpol %>% dplyr::rename(mpol = pvalue)
pvalue_top_N_MDR_WW$mpal <- pvalue_top_N_MDR_WW$mpal %>% dplyr::rename(mpal = pvalue)
pvalue_top_N_MDR_WW$physco <- pvalue_top_N_MDR_WW$physco %>% dplyr::rename(physco = pvalue)

# merge the dataframes to one
pvalue_top_N_MDR_WW_merged <- merge(pvalue_top_N_MDR_WW$mpol, pvalue_top_N_MDR_WW$mpal, by = "ASV", all = TRUE)
pvalue_top_N_MDR_WW_merged <- merge(pvalue_top_N_MDR_WW_merged, pvalue_top_N_MDR_WW$physco, by = "ASV", all = TRUE)
pvalue_top_N_MDR_WW_merged
rownames(pvalue_top_N_MDR_WW_merged) <- pvalue_top_N_MDR_WW_merged$ASV # put ASV as rownames
pvalue_top_N_MDR_WW_merged <- pvalue_top_N_MDR_WW_merged[, which(colnames(pvalue_top_N_MDR_WW_merged) != "ASV")] # remove column ASV




## FDR-correct p-values for N tests (the number of ASVs that are shown in the heatmap)
pvalue_top_N_MDR_WW_merged$mpol <- stats::p.adjust(p = pvalue_top_N_MDR_WW_merged$mpol, method = "fdr") # "fdr" = Benjamini-Hochberg
pvalue_top_N_MDR_WW_merged$mpal <- stats::p.adjust(p = pvalue_top_N_MDR_WW_merged$mpal, method = "fdr") # "fdr" = Benjamini-Hochberg
pvalue_top_N_MDR_WW_merged$physco <- stats::p.adjust(p = pvalue_top_N_MDR_WW_merged$physco, method = "fdr") # "fdr" = Benjamini-Hochberg




## turn p-values into significance asterisks

# now transform p-values into asterisks
asterisks_top_N_MDR_WW_merged <- pvalue_top_N_MDR_WW_merged
asterisks_top_N_MDR_WW_merged[pvalue_top_N_MDR_WW_merged > 0.05] <- ""
asterisks_top_N_MDR_WW_merged[pvalue_top_N_MDR_WW_merged <= 0.05] <- "*"
asterisks_top_N_MDR_WW_merged[is.na(pvalue_top_N_MDR_WW_merged)] <- ""





# order and rename columns for heatmap
columns_order <- c("mpal", "mpol", "physco")
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[, columns_order]
asterisks_top_N_MDR_WW_merged <- asterisks_top_N_MDR_WW_merged[, columns_order]
colnames(LFC_top_N_MDR_WW_merged) <- c("M. paleacea", "M. polymorpha", "P. patens")
colnames(asterisks_top_N_MDR_WW_merged) <- c("M. paleacea", "M. polymorpha", "P. patens")


# order rows for heatmap
# rows should be in the same order as the ASV phylogeny (hclust object)
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[fun_phylogeny_top_N$labels,]
asterisks_top_N_MDR_WW_merged <- asterisks_top_N_MDR_WW_merged[fun_phylogeny_top_N$labels,]




# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- -3
mid1 <- 0
max1 <- 3
myBreaks <- c(seq(min1, mid1, length.out=ceiling(paletteLength/2) + 1),
              seq(mid1+0.1, max1, length.out=floor(paletteLength/2)))




# re-order species for consistency (first P. patens, then M. paleacea, then M. polymorpha)
LFC_top_N_MDR_WW_merged <- LFC_top_N_MDR_WW_merged[, c(3, 1, 2)]
asterisks_top_N_MDR_WW_merged <- asterisks_top_N_MDR_WW_merged[, c(3, 1, 2)]
colnames(LFC_top_N_MDR_WW_merged) == colnames(asterisks_top_N_MDR_WW_merged)

# make column labels with italicized species names
hm_column_labels_italicized <- expression(italic("P. patens"), italic("M. paleacea"), italic("M. polymorpha"))

# make heatmap
hm <- pheatmap(LFC_top_N_MDR_WW_merged, cluster_cols=FALSE, cluster_rows=fun_phylogeny_top_N, labels_row=genus_labels_fun_top_N_italicized, labels_col=hm_column_labels_italicized, display_numbers=asterisks_top_N_MDR_WW_merged, color = myColor, breaks = myBreaks, na_col = "grey")

# export plot for Supp. Fig. 4
ggsave(plot=hm, filename="R_plots/fun_top_100_LFC_mild_drought.png", width=1300, height=5000, unit="px")





## 9 - Rhizobiales phylogenetic heatmap (Fig. 2C) ----

## 9.1 - Subset top N Rhizobiales ASVs ====

# subset Rhizobiales
rhizobiales_ps <- subset_taxa(filtered_relative_bac_ps, Order == "o__Rhizobiales")

# how many ASVs?
top_N <- 50

# list top N ASVs
rhizobiales_top_N_ASVs <- names(sort(taxa_sums(rhizobiales_ps), decreasing=TRUE)[1:top_N])

# filter top N ASVs
rhizobiales_ps_top_N <- prune_taxa(rhizobiales_top_N_ASVs, filtered_relative_bac_ps)



## 9.2 - Prepare phylogeny ====

# load bacterial phylogeny
bac_phylogeny <- read.tree("./data_bryophytes/16S/tree_Bac.nwk")
# plot.phylo(bac_phylogeny) # plot tree

# remove extra apostrophes from tip labels
bac_phylogeny$tip.label <- gsub("'", "", bac_phylogeny$tip.label)

# prune phylogeny to ASVs present in table
rhizobiales_phylogeny_top_N <- keep.tip(bac_phylogeny, rhizobiales_top_N_ASVs)
plot.phylo(rhizobiales_phylogeny_top_N) # plot tree

# align tree tips to make dendogram later
rhizobiales_phylogeny_top_N <- chronos(rhizobiales_phylogeny_top_N)

# turn the tree to a hclust object
# rhizobiales_phylogeny_top_N <- RRphylo::fix.poly(rhizobiales_phylogeny_top_N, type="resolve") # resolve polytomies to non-zero branch lengths, as the function as.hclust() only takes binary trees
rhizobiales_phylogeny_top_N <- as.hclust(rhizobiales_phylogeny_top_N) # convert to hclust object
plot(rhizobiales_phylogeny_top_N)



## prepare Genus names for rownames

# first let's load the tax table (output from qiime2-sklearn)
raw_bac_taxtab <- read.table(file = "C:/bryophytes_rhizosphere/data_bryophytes/16S/taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
raw_bac_taxtab <- separate(data = raw_bac_taxtab, col = Taxon, into =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

# get the genus names of the ASVs of our data frames (via rownames)
genus_labels_rhizobiales_top_N <- raw_bac_taxtab[rhizobiales_phylogeny_top_N$labels, "Genus"]

# remove "g__" before genus names
genus_labels_rhizobiales_top_N <- gsub(" g__", "", genus_labels_rhizobiales_top_N)


# change NA in genus labels to "NA" as text
genus_labels_rhizobiales_top_N[is.na(genus_labels_rhizobiales_top_N)] <- "NA"

# find which labels should be italicized: all genus labels except "uncultured" and "NA"
label_numbers_to_italicize <- which(!(genus_labels_rhizobiales_top_N %in% c("uncultured", "NA")))

# copy genus labels list to italicized genus labels list
genus_labels_rhizobiales_top_N_italicized <- genus_labels_rhizobiales_top_N

# italicize genus labels except "uncultured" and "NA"
for(i in label_numbers_to_italicize){
  genus_labels_rhizobiales_top_N_italicized[i] <- as.expression(bquote(italic(.(genus_labels_rhizobiales_top_N[i]))))
}



## 9.3 - Heatmap of relative abundance in each species/treatment ====

# use psmelt to combine the ASV table with metadata and taxonomy
rhizobiales_top_N_melted <- psmelt(rhizobiales_ps_top_N)

# calculate average relative abundance per species/treatment
rhizobiales_top_N_species_treatment <- rhizobiales_top_N_melted %>% 
  group_by(OTU, plant_species, treatment) %>%
  summarise(Abundance = mean(Abundance))



# make concatenated species/treatment column
rhizobiales_top_N_species_treatment$species_treatment <- paste0(rhizobiales_top_N_species_treatment$plant_species, "_", rhizobiales_top_N_species_treatment$treatment)

# remove columns plant_species and treatment
rhizobiales_top_N_species_treatment <- rhizobiales_top_N_species_treatment[, which(colnames(rhizobiales_top_N_species_treatment) %in% c("OTU", "species_treatment", "Abundance"))]

# pivot to wide format: ASVs as rows, species/treatment as columns
rhizobiales_top_N_species_treatment <- rhizobiales_top_N_species_treatment %>% pivot_wider(names_from='species_treatment',
                                                                           values_from='Abundance')
# change to dataframe
rhizobiales_top_N_species_treatment <- as.data.frame(rhizobiales_top_N_species_treatment)

# move ASV names to rownames
rownames(rhizobiales_top_N_species_treatment) <- rhizobiales_top_N_species_treatment$OTU

# remove column "OTU"
rhizobiales_top_N_species_treatment <- rhizobiales_top_N_species_treatment[,which(colnames(rhizobiales_top_N_species_treatment) != "OTU")]

# order and rename columns for heatmap
columns_order <- c("mpal_20", "mpal_13", "mpol_20", "mpol_13", "mpol_10", "physco_20", "physco_13")
rhizobiales_top_N_species_treatment <- rhizobiales_top_N_species_treatment[, columns_order]
colnames(rhizobiales_top_N_species_treatment) <- c("M. paleacea - well-watered", "M. paleacea - mild drought", "M. polymorpha - well-watered", "M. polymorpha - mild drought", "M. polymorpha - severe drought", "P. patens - well-watered", "P. patens - mild drought")





# order rows for heatmap
# rows should be in the same order as the ASV phylogeny (hclust object)
rhizobiales_top_N_species_treatment <- rhizobiales_top_N_species_treatment[rhizobiales_phylogeny_top_N$labels,]



# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 0.01
myBreaks <- c(seq(min1, max1, length=paletteLength))




# make heatmap
pheatmap(rhizobiales_top_N_species_treatment, cluster_cols=FALSE, cluster_rows=rhizobiales_phylogeny_top_N, labels_row=genus_labels_rhizobiales_top_N, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")



## 9.4 - Heatmap of normalized relative abundance in each species/treatment ====

# normalized by ASV-average, to see the host species and treatment preference of each ASV

# normalize rel. abundance by average of each ASV across plant species
rhizobiales_top_N_species_treatment_normalized <- rhizobiales_top_N_species_treatment*ncol(rhizobiales_top_N_species_treatment)/rowSums(rhizobiales_top_N_species_treatment)



## define colours and breaks for heatmaps - for normalized relative abundance

# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 3
myBreaks <- c(seq(min1, max1, length=paletteLength))


# re-order species/treatments for consistency (first P. patens, then M. paleacea, then M. polymorpha)
rhizobiales_top_N_species_treatment_normalized <- rhizobiales_top_N_species_treatment_normalized[, c(6, 7, 1:5)]

# make column labels with italicized species names
hm_column_labels_italicized <- expression(italic("P. patens") ~ "- well-watered", italic("P. patens") ~ "- mild drought", italic("M. paleacea") ~ "- well-watered", italic("M. paleacea") ~ "- mild drought", italic("M. polymorpha") ~ "- well-watered", italic("M. polymorpha") ~ "- mild drought", italic("M. polymorpha") ~ "- severe drought")

# make heatmap (not transposed, vertical)
hm <- pheatmap(rhizobiales_top_N_species_treatment_normalized, cluster_cols=FALSE, cluster_rows=rhizobiales_phylogeny_top_N, labels_row=genus_labels_rhizobiales_top_N_italicized, labels_col=hm_column_labels_italicized, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")
# export plot for Fig. 2C
ggsave(plot=hm, filename="R_plots/Rhizobiales_host_treatment_preferences.png", width=2000, height=2300, unit="px")



## 10 - Fusarium phylogenetic heatmap (Fig. 3B) ----

## 10.1 - Subset top N Fusarium ASVs ====

# subset Fusarium
fusarium_ps <- subset_taxa(filtered_relative_fun_ps, Genus == "g__Fusarium")

# how many ASVs?
top_N <- 5

# list top N ASVs
fusarium_top_N_ASVs <- names(sort(taxa_sums(fusarium_ps), decreasing=TRUE)[1:top_N])

# filter top N ASVs
fusarium_ps_top_N <- prune_taxa(fusarium_top_N_ASVs, filtered_relative_fun_ps)



## 10.2 - Prepare phylogeny ====

# load fungal phylogeny
fun_phylogeny <- read.tree("./data_bryophytes/ITS/tree_Fun.nwk")
# plot.phylo(fun_phylogeny) # plot tree

# remove extra apostrophes from tip labels
fun_phylogeny$tip.label <- gsub("'", "", fun_phylogeny$tip.label)

# prune phylogeny to ASVs present in table
fusarium_phylogeny_top_N <- keep.tip(fun_phylogeny, fusarium_top_N_ASVs)
plot.phylo(fusarium_phylogeny_top_N) # plot tree

# align tree tips to make dendogram later
fusarium_phylogeny_top_N <- chronos(fusarium_phylogeny_top_N)

# turn the tree to a hclust object
# fusarium_phylogeny_top_N <- RRphylo::fix.poly(fusarium_phylogeny_top_N, type="resolve") # resolve polytomies to non-zero branch lengths, as the function as.hclust() only takes binary trees
fusarium_phylogeny_top_N <- as.hclust(fusarium_phylogeny_top_N) # convert to hclust object
plot(fusarium_phylogeny_top_N)





## 10.3 - Heatmap of relative abundance in each species/treatment ====

# use psmelt to combine the ASV table with metadata and taxonomy
fusarium_top_N_melted <- psmelt(fusarium_ps_top_N)

# calculate average relative abundance per species/treatment
fusarium_top_N_species_treatment <- fusarium_top_N_melted %>% 
  group_by(OTU, plant_species, treatment) %>%
  summarise(Abundance = mean(Abundance))



# make concatenated species/treatment column
fusarium_top_N_species_treatment$species_treatment <- paste0(fusarium_top_N_species_treatment$plant_species, "_", fusarium_top_N_species_treatment$treatment)

# remove columns plant_species and treatment
fusarium_top_N_species_treatment <- fusarium_top_N_species_treatment[, which(colnames(fusarium_top_N_species_treatment) %in% c("OTU", "species_treatment", "Abundance"))]

# pivot to wide format: ASVs as rows, species/treatment as columns
fusarium_top_N_species_treatment <- fusarium_top_N_species_treatment %>% pivot_wider(names_from='species_treatment',
                                                                                             values_from='Abundance')
# change to dataframe
fusarium_top_N_species_treatment <- as.data.frame(fusarium_top_N_species_treatment)

# move ASV names to rownames
rownames(fusarium_top_N_species_treatment) <- fusarium_top_N_species_treatment$OTU

# remove column "OTU"
fusarium_top_N_species_treatment <- fusarium_top_N_species_treatment[,which(colnames(fusarium_top_N_species_treatment) != "OTU")]

# order and rename columns for heatmap
columns_order <- c("mpal_20", "mpal_13", "mpol_20", "mpol_13", "mpol_10", "physco_20", "physco_13")
fusarium_top_N_species_treatment <- fusarium_top_N_species_treatment[, columns_order]
colnames(fusarium_top_N_species_treatment) <- c("M. paleacea - well-watered", "M. paleacea - mild drought", "M. polymorpha - well-watered", "M. polymorpha - mild drought", "M. polymorpha - severe drought", "P. patens - well-watered", "P. patens - mild drought")



# order rows for heatmap
# rows should be in the same order as the ASV phylogeny (hclust object)
fusarium_top_N_species_treatment <- fusarium_top_N_species_treatment[fusarium_phylogeny_top_N$labels,]



# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 0.01
myBreaks <- c(seq(min1, max1, length=paletteLength))


# make heatmap
pheatmap(fusarium_top_N_species_treatment, cluster_cols=FALSE, cluster_rows=fusarium_phylogeny_top_N, labels_row=rownames(fusarium_top_N_species_treatment), display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")




## 10.4 - Heatmap of normalized relative abundance in each species/treatment ====

# normalized by ASV-average, to see the host species and treatment preference of each ASV

# normalize rel. abundance by average of each ASV across plant species
fusarium_top_N_species_treatment_normalized <- fusarium_top_N_species_treatment*ncol(fusarium_top_N_species_treatment)/rowSums(fusarium_top_N_species_treatment)



## define colours and breaks for heatmaps - for normalized relative abundance

# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 3
myBreaks <- c(seq(min1, max1, length=paletteLength))



# re-order species/treatments for consistency (first P. patens, then M. paleacea, then M. polymorpha)
fusarium_top_N_species_treatment_normalized <- fusarium_top_N_species_treatment_normalized[, c(6, 7, 1:5)]

# make column labels with italicized species names
hm_column_labels_italicized <- expression(italic("P. patens") ~ "- well-watered", italic("P. patens") ~ "- mild drought", italic("M. paleacea") ~ "- well-watered", italic("M. paleacea") ~ "- mild drought", italic("M. polymorpha") ~ "- well-watered", italic("M. polymorpha") ~ "- mild drought", italic("M. polymorpha") ~ "- severe drought")

# make heatmap (not transposed, vertical)
hm <- pheatmap(fusarium_top_N_species_treatment_normalized, cluster_cols=FALSE, cluster_rows=fusarium_phylogeny_top_N, labels_row=rownames(fusarium_top_N_species_treatment_normalized), labels_col=hm_column_labels_italicized, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")

# export plot for Fig. 3B
ggsave(plot=hm, filename="R_plots/Fusarium_host_treatment_preferences.png", width=1000, height=1000, unit="px")







## 11 - Glomeromycota phylogenetic heatmap (Supp. Fig. 8) ----


## 11.1 - Subset top N Glomeromycota ASVs ====

# subset Glomeromycota
glomeromycota_ps <- subset_taxa(filtered_relative_fun_ps, Phylum == "p__Glomeromycota")

glomeromycota_ASVs <- rownames(otu_table(glomeromycota_ps))


## 11.2 - Prepare phylogeny ====

# load fungal phylogeny
fun_phylogeny <- read.tree("./data_bryophytes/ITS/tree_Fun.nwk")
# plot.phylo(fun_phylogeny) # plot tree

# remove extra apostrophes from tip labels
fun_phylogeny$tip.label <- gsub("'", "", fun_phylogeny$tip.label)

# prune phylogeny to ASVs present in table
glomeromycota_phylogeny <- keep.tip(fun_phylogeny, glomeromycota_ASVs)
plot.phylo(glomeromycota_phylogeny) # plot tree

# align tree tips to make dendogram later
glomeromycota_phylogeny <- chronos(glomeromycota_phylogeny)

# turn the tree to a hclust object
# glomeromycota_phylogeny <- RRphylo::fix.poly(glomeromycota_phylogeny, type="resolve") # resolve polytomies to non-zero branch lengths, as the function as.hclust() only takes binary trees
glomeromycota_phylogeny <- as.hclust(glomeromycota_phylogeny) # convert to hclust object
plot(glomeromycota_phylogeny)




## prepare Genus names for rownames

# first let's load the tax table (output from qiime2-sklearn)
raw_fun_taxtab <- read.table(file = "C:/bryophytes_rhizosphere/data_bryophytes/ITS/taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
raw_fun_taxtab <- separate(data = raw_fun_taxtab, col = Taxon, into =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

# add "f_" prefix to fungal ASV names, because the original taxonomy table didn't have that prefix
rownames(raw_fun_taxtab) <- paste0("f", rownames(raw_fun_taxtab))

# get the genus names of the ASVs of our data frames (via rownames)
genus_labels_glomeromycota <- raw_fun_taxtab[glomeromycota_phylogeny$labels, "Genus"]

# remove "g__" before genus names
genus_labels_glomeromycota <- gsub("g__", "", genus_labels_glomeromycota)

# find which labels should be italicized: all genus labels except "unidentified"
label_numbers_to_italicize <- which(!(genus_labels_glomeromycota == "unidentified"))

# copy genus labels list to italicized genus labels list
genus_labels_glomeromycota_italicized <- genus_labels_glomeromycota

# italicize genus labels except "uncultured" and "NA"
for(i in label_numbers_to_italicize){
  genus_labels_glomeromycota_italicized[i] <- as.expression(bquote(italic(.(genus_labels_glomeromycota[i]))))
}



## 11.3 - Heatmap of relative abundance in each species/treatment ====

# use psmelt to combine the ASV table with metadata and taxonomy
glomeromycota_melted <- psmelt(glomeromycota_ps)

# calculate average relative abundance per species/treatment
glomeromycota_species_treatment <- glomeromycota_melted %>% 
  group_by(OTU, plant_species, treatment) %>%
  summarise(Abundance = mean(Abundance))



# make concatenated species/treatment column
glomeromycota_species_treatment$species_treatment <- paste0(glomeromycota_species_treatment$plant_species, "_", glomeromycota_species_treatment$treatment)

# remove columns plant_species and treatment
glomeromycota_species_treatment <- glomeromycota_species_treatment[, which(colnames(glomeromycota_species_treatment) %in% c("OTU", "species_treatment", "Abundance"))]

# pivot to wide format: ASVs as rows, species/treatment as columns
glomeromycota_species_treatment <- glomeromycota_species_treatment %>% pivot_wider(names_from='species_treatment',
                                                                                     values_from='Abundance')
# change to dataframe
glomeromycota_species_treatment <- as.data.frame(glomeromycota_species_treatment)

# move ASV names to rownames
rownames(glomeromycota_species_treatment) <- glomeromycota_species_treatment$OTU

# remove column "OTU"
glomeromycota_species_treatment <- glomeromycota_species_treatment[,which(colnames(glomeromycota_species_treatment) != "OTU")]

# order and rename columns for heatmap
columns_order <- c("mpal_20", "mpal_13", "mpol_20", "mpol_13", "mpol_10", "physco_20", "physco_13")
glomeromycota_species_treatment <- glomeromycota_species_treatment[, columns_order]
colnames(glomeromycota_species_treatment) <- c("M. paleacea - well-watered", "M. paleacea - mild drought", "M. polymorpha - well-watered", "M. polymorpha - mild drought", "M. polymorpha - severe drought", "P. patens - well-watered", "P. patens - mild drought")



# order rows for heatmap
# rows should be in the same order as the ASV phylogeny (hclust object)
glomeromycota_species_treatment <- glomeromycota_species_treatment[glomeromycota_phylogeny$labels,]



# define colour settings for heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)

# use floor and ceiling to deal with even/odd length palettelengths
min1 <- 0
max1 <- 0.002
myBreaks <- c(seq(min1, max1, length=paletteLength))



# re-order species/treatments for consistency (first P. patens, then M. paleacea, then M. polymorpha)
glomeromycota_species_treatment <- glomeromycota_species_treatment[, c(6, 7, 1:5)]

# make column labels with italicized species names
hm_column_labels_italicized <- expression(italic("P. patens") ~ "- well-watered", italic("P. patens") ~ "- mild drought", italic("M. paleacea") ~ "- well-watered", italic("M. paleacea") ~ "- mild drought", italic("M. polymorpha") ~ "- well-watered", italic("M. polymorpha") ~ "- mild drought", italic("M. polymorpha") ~ "- severe drought")

# make heatmap (not transposed, vertical)
hm <- pheatmap(glomeromycota_species_treatment, cluster_cols=FALSE, cluster_rows=glomeromycota_phylogeny, labels_row=genus_labels_glomeromycota_italicized, labels_col=hm_column_labels_italicized, display_numbers=FALSE, color = myColor, breaks = myBreaks, na_col = "grey")

# export plot for Supp. Fig. 8
ggsave(plot=hm, filename="R_plots/Glomeromycota_ASVs_relative_abundance.png", width=1000, height=1300, unit="px")







## 12 - LEfSe analysis ----



# we conduct LEfSe analysis using the lefser package
# first we need to convert our phyloseq object to a SummarizedExperiment object

# load filtered physeq objects with counts
load("./data_bryophytes/16S/phyloseq_objects/filtered_bac_ps.RData")
load("./data_bryophytes/ITS/phyloseq_objects/filtered_fun_ps.RData")


# extract ASV table and sample data from phyloseq
counts_bac <- unclass(otu_table(filtered_bac_ps))
colData_bac <- as(sample_data(filtered_bac_ps), "data.frame")

counts_fun <- unclass(otu_table(filtered_fun_ps))
colData_fun <- as(sample_data(filtered_fun_ps), "data.frame")

# create a SummarizedExperiment objects
filtered_bac_SE <- SummarizedExperiment(assays = list(counts = counts_bac), colData = colData_bac)
filtered_fun_SE <- SummarizedExperiment(assays = list(counts = counts_fun), colData = colData_fun)


## to be moved to top of script
library(lefser)


res <- lefser(filtered_fun_SE, groupCol = "plant_species", blockCol = "treatment")
