

# load libraries
library(ggplot2)
# library(dplyr)
library(car) # for leveneTest()
library(dunn.test)
library(multcompView) # for significance letters

setwd("C:/bryophytes_rhizosphere/for_publication/code_for_github/angiosperms_experiment/angiosperms_biomass")


## data loading and preparation ----

# load csv
data <- read.csv("angiosperms_biomass.csv", header=TRUE, sep = ";") 


# note on treatment codes:
# "15" = over-watered = 13.0% (w/w) soil moisture
# "12" = well-watered = 10.7% (w/w) soil moisture
# "9" = mild drought = 8.3% (w/w) soil moisture
# "6" = severe drought = 5.7% (w/w) soil moisture


# change weights to numeric
data$shoot_fresh_weight <- as.numeric(data$shoot_fresh_weight)
data$shoot_dry_weight <- as.numeric(data$shoot_dry_weight)
data$root_dry_weight <- as.numeric(data$root_dry_weight)

# change treatment to factor
# data$treatment <- as.factor(data$treatment)

# calculate root/shoot ratio (dry weight)
data$root_shoot_ratio_dry <- data$root_dry_weight / data$shoot_dry_weight

# # calculate total dry weight
# data$total_dry_weight <- data$shoot_dry_weight + data$root_dry_weight


# make new column with the soil moisture content
data$soil_moisture <- ""
data$soil_moisture[which(data$treatment=="15")] <- "13.0%"
data$soil_moisture[which(data$treatment=="12")] <- "10.7%"
data$soil_moisture[which(data$treatment=="9")] <- "8.3%"
data$soil_moisture[which(data$treatment=="6")] <- "5.7%"

# make new column with the treatment name
data$treatment_name <- ""
data$treatment_name[which(data$treatment=="12")] <- "well-watered"
data$treatment_name[which(data$treatment=="9")] <- "mild drought"
data$treatment_name[which(data$treatment=="6")] <- "severe drought"

# order treatments
data$treatment <- factor(data$treatment, levels = c("15", "12", "9", "6"))
data$soil_moisture <- factor(data$soil_moisture, levels = c("13.0%", "10.7%", "8.3%", "5.7%"))
data$treatment_name <- factor(data$treatment_name, levels = c("well-watered", "mild drought", "severe drought"))

# change species names to latin names
data$plant_species <- gsub("arabidopsis", "A. thaliana", data$plant_species)
data$plant_species <- gsub("tomato", "S. lycopersicum", data$plant_species)


# make column "species_treatment"
data$species_treatment <- paste0(data$plant_species, " / ", data$soil_moisture)




##  shoot dry weight ----

# test assumption for ANOVA: equal variances (homoscedasticity) - Levene's test
leveneTest((shoot_dry_weight) ~ treatment * plant_species, data = data) # significant heteroscedasticity

# significant heteroscedasticity, we will use a non-parametric test

# Kruskal-Wallis test
kw <- kruskal.test(shoot_dry_weight ~ species_treatment, data = data)
# check p-value of Kruskal-Wallis
kw$p.value # p < 0.05

# Dunn's test
dunn <- dunn.test(data$shoot_dry_weight, data$species_treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)

# put p-values in a named vector for multcompLetters()
pvec <- dunn$altP.adjusted
names(pvec) <- dunn$comparisons
# remove spaces around dash because it messes up multcompLetters()
names(pvec) <- gsub(" - ", "-", names(pvec))


# define group (species_treatment) order for significance letter computation
group_order <- c("A. thaliana / 13.0%", "A. thaliana / 10.7%", "A. thaliana / 8.3%", "A. thaliana / 5.7%", "S. lycopersicum / 13.0%", "S. lycopersicum / 10.7%", "S. lycopersicum / 8.3%", "S. lycopersicum / 5.7%")
data$species_treatment <- factor(data$species_treatment, levels = group_order)

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

# make a dataframe with significance letters
letters <- data.frame(cld)
# make column "group"
letters$group <- rownames(letters)
# re-order the table in the desired group order (needed to match the order on the boxplot)
letters <- letters[group_order,]
# add columns plant_species and treatment
letters$plant_species <- factor(c("A. thaliana", "A. thaliana", "A. thaliana", "A. thaliana", "S. lycopersicum", "S. lycopersicum", "S. lycopersicum", "S. lycopersicum"))
letters$treatment <- factor(c("13.0%", "10.7%", "8.3%", "5.7%", "13.0%", "10.7%", "8.3%", "5.7%"))


# boxplot
p <- ggplot(data, aes(x=soil_moisture, y=shoot_dry_weight)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2) +
  facet_wrap(~plant_species) +
  labs(title=NULL, x =NULL, y = "shoot dry weight [g]") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(plot.title = element_text(hjust = 0.5)) + # center main title
  expand_limits(y = 0) + # force y-axis to show 0
  theme(text = element_text(size = 10), strip.text = element_text(size = 13)) +  # increase text size for all text elements, and facet title size even bigger
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic

# define y-positions for labels: extract y-value of the top of the whisker of each boxplot from the ggplot object
letters$ypos <- layer_data(p)$ymax_final
# get y-axis range
yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
# add 5% of the y axis range to have the label slightly higher than the top of the whisker
letters$ypos <- letters$ypos + 0.1 * yrange

# show minimal boxplots with significance letters
p + geom_text(data = letters, aes(x = treatment, y = ypos, label = cld))


# export plot for Supp. Fig. 10A
ggsave("R_plots/angiosperms_shoot_dry_weight.png", width=1700, height=1000, unit="px")



## root dry weight ----

# for the root dry weight we only have measurements for the well-watered, mild drought and severe drought treatments
# exclude treatment "15" (overwatered, officially 13.0% soil moisture content) for root dry weight and root/shoot ratio analysis
data <- subset(data, treatment != "15")


# test assumption for ANOVA: equal variances (homoscedasticity) - Levene's test
leveneTest((root_dry_weight) ~ treatment * plant_species, data = data) # no significant heteroscedasticity

# no significant heteroscedasticity, but we will use a non-parametric test for consistency (as we did for all the others)

# Kruskal-Wallis test
kw <- kruskal.test(root_dry_weight ~ species_treatment, data = data)
# check p-value of Kruskal-Wallis
kw$p.value # p < 0.05

# Dunn's test
dunn <- dunn.test(data$root_dry_weight, data$species_treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)

# put p-values in a named vector for multcompLetters()
pvec <- dunn$altP.adjusted
names(pvec) <- dunn$comparisons
# remove spaces around dash because it messes up multcompLetters()
names(pvec) <- gsub(" - ", "-", names(pvec))


# define group (species_treatment) order for significance letter computation
group_order <- c("A. thaliana / 10.7%", "A. thaliana / 8.3%", "A. thaliana / 5.7%", "S. lycopersicum / 10.7%", "S. lycopersicum / 8.3%", "S. lycopersicum / 5.7%")
data$species_treatment <- factor(data$species_treatment, levels = group_order)

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

# make a dataframe with significance letters
letters <- data.frame(cld)
# make column "group"
letters$group <- rownames(letters)
# re-order the table in the desired group order (needed to match the order on the boxplot)
letters <- letters[group_order,]
# add columns plant_species and treatment
letters$plant_species <- factor(c("A. thaliana", "A. thaliana", "A. thaliana", "S. lycopersicum", "S. lycopersicum", "S. lycopersicum"))
letters$treatment <- factor(c("10.7%", "8.3%", "5.7%", "10.7%", "8.3%", "5.7%"))



# boxplot
p <- ggplot(data, aes(x=soil_moisture, y=root_dry_weight)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2) +
  facet_wrap(~plant_species) +
  labs(title=NULL, x =NULL, y = "root dry weight [g]") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(plot.title = element_text(hjust = 0.5)) + # center main title
  expand_limits(y = 0) + # force y-axis to show 0
  theme(text = element_text(size = 10), strip.text = element_text(size = 13)) +  # increase text size for all text elements, and facet title size even bigger
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic

# define y-positions for labels: extract y-value of the top of the whisker of each boxplot from the ggplot object
letters$ypos <- layer_data(p)$ymax_final
# get y-axis range
yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
# add 5% of the y axis range to have the label slightly higher than the top of the whisker
letters$ypos <- letters$ypos + 0.1 * yrange

# show minimal boxplots with significance letters
p + geom_text(data = letters, aes(x = treatment, y = ypos, label = cld))


# export plot for Supp. Fig. 10B
ggsave("R_plots/angiosperms_root_dry_weight.png", width=1700, height=1000, unit="px")






## root/shoot ratio (with dry weights) ----

# for the root dry weight we only have measurements for the well-watered, mild drought and severe drought treatments
# exclude treatment "15" (overwatered, officially 13.0% soil moisture content) for root dry weight and root/shoot ratio analysis
data <- subset(data, treatment != "15")


# test assumption for ANOVA: equal variances (homoscedasticity) - Levene's test
leveneTest((root_shoot_ratio_dry) ~ treatment * plant_species, data = data) # significant heteroscedasticity

# we will use a non-parametric test for the root/shoot ratio

# Kruskal-Wallis test
kw <- kruskal.test(root_shoot_ratio_dry ~ species_treatment, data = data)
# check p-value of Kruskal-Wallis
kw$p.value


# Dunn's test
dunn <- dunn.test(data$root_shoot_ratio_dry, data$species_treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)

# put p-values in a named vector for multcompLetters()
pvec <- dunn$altP.adjusted
names(pvec) <- dunn$comparisons
# remove spaces around dash because it messes up multcompLetters()
names(pvec) <- gsub(" - ", "-", names(pvec))



# define group (species_treatment) order for significance letter computation
group_order <- c("A. thaliana / 10.7%", "A. thaliana / 8.3%", "A. thaliana / 5.7%", "S. lycopersicum / 10.7%", "S. lycopersicum / 8.3%", "S. lycopersicum / 5.7%")
data$species_treatment <- factor(data$species_treatment, levels = group_order)

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

# make a dataframe with significance letters
letters <- data.frame(cld)
# make column "group"
letters$group <- rownames(letters)
# re-order the table in the desired group order (needed to match the order on the boxplot)
letters <- letters[group_order,]
# add columns plant_species and treatment
letters$plant_species <- factor(c("A. thaliana", "A. thaliana", "A. thaliana", "S. lycopersicum", "S. lycopersicum", "S. lycopersicum"))
letters$treatment_name <- factor(c("well-watered", "mild drought", "severe drought", "well-watered", "mild drought", "severe drought"))



# boxplot
p <- ggplot(data, aes(x=treatment_name, y=root_shoot_ratio_dry)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2) +
  facet_wrap(~plant_species) +
  labs(title=NULL, x =NULL, y = "root/shoot ratio [g/g]") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(plot.title = element_text(hjust = 0.5)) + # center main title
  expand_limits(y = 0) + # force y-axis to show 0
  theme(text = element_text(size = 12), strip.text = element_text(size = 13)) +  # increase text size for all text elements, and facet title size even bigger
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # angle x-axis labels
  theme(strip.background = element_blank()) + # no box around facet label
  theme(strip.text = element_text(face = "italic")) # facet label in italic

# define y-positions for labels: extract y-value of the top of the whisker of each boxplot from the ggplot object
letters$ypos <- layer_data(p)$ymax_final
# get y-axis range
yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
# add 5% of the y axis range to have the label slightly higher than the top of the whisker
letters$ypos <- letters$ypos + 0.1 * yrange

# show minimal boxplots with significance letters
p + geom_text(data = letters, aes(x = treatment_name, y = ypos, label = cld))


# export plot for Fig. 4C
ggsave("R_plots/angiosperms_root_shoot_ratio.png", width=1400, height=1000, unit="px")


