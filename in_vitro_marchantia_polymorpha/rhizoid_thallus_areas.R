


# load libraries
library(ggplot2)
library(car) # for leveneTest()
library(dunn.test)
library(multcompView) # for significance letters

setwd("C:/bryophytes_rhizosphere/for_publication/code_for_github/in_vitro_marchantia_polymorpha")

# read csv
data <- read.csv("rhizoid_thallus_areas_14days.csv", sep=";")

# calculate rhizoid/thallus area ratio
data$rhizoid_thallus_ratio <- data$rhizoid_area / data$thallus_area

# replace uM with μM for plot
data$treatment <- gsub("uM", "μM", data$treatment)

# order treatments for boxplots
data$treatment <- factor(data$treatment, levels=c("control", "100 μM sorbitol", "150 μM sorbitol"))


# define group order
group_order <- c("control", "100 μM sorbitol", "150 μM sorbitol")

# order the groups for the significance letters and plot
data <- data[order(factor(data$treatment, levels = group_order)),]
data$treatment <- factor(data$treatment, levels=group_order)




## thallus area ----

# test assumption for ANOVA: equal variances (homoscedasticity) - Levene's test
leveneTest((thallus_area) ~ treatment, data = data)

# significant heteroscedasticity, we will use a non-parametric test

# Kruskal-Wallis test
kw <- kruskal.test(thallus_area ~ treatment, data = data)
# check p-value of Kruskal-Wallis
kw$p.value # p < 0.05

# Dunn's test
dunn <- dunn.test(data$thallus_area, data$treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)

# put p-values in a named vector for multcompLetters()
pvec <- dunn$altP.adjusted
names(pvec) <- dunn$comparisons
# remove spaces around dash because it messes up multcompLetters()
names(pvec) <- gsub(" - ", "-", names(pvec))


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


# boxplot
p <- ggplot(data = data, aes(x = treatment, y = thallus_area)) +
  labs(x = "treatment", y = "thallus area [mm²]") +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 0)  +
  xlab(NULL) +
  theme(text = element_text(size = 10))  # increase text size for all text elements

# define y-positions for labels: extract y-value of the top of the whisker of each boxplot from the ggplot object
letters$ypos <- layer_data(p)$ymax_final
# get y-axis range
yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
# add 5% of the y axis range to have the label slightly higher than the top of the whisker
letters$ypos <- letters$ypos + 0.1 * yrange

# show minimal boxplots with significance letters
p + geom_text(data = letters, aes(x = group, y = ypos, label = cld))


# export plot for Supp. Fig. 11A
ggsave("R_plots/polymorpha_thallus_area.png", width=1000, height=1000, unit="px")




## rhizoid area ----


# test assumption for ANOVA: equal variances (homoscedasticity) - Levene's test
leveneTest((rhizoid_area) ~ treatment, data = data) # significant heteroscedasticity

# significant heteroscedasticity, we will use a non-parametric test

# Kruskal-Wallis test
kw <- kruskal.test(rhizoid_area ~ treatment, data = data)
# check p-value of Kruskal-Wallis
kw$p.value # p < 0.05

# Dunn's test
dunn <- dunn.test(data$rhizoid_area, data$treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)

# put p-values in a named vector for multcompLetters()
pvec <- dunn$altP.adjusted
names(pvec) <- dunn$comparisons
# remove spaces around dash because it messes up multcompLetters()
names(pvec) <- gsub(" - ", "-", names(pvec))


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


# boxplot
p <- ggplot(data = data, aes(x = treatment, y = rhizoid_area)) +
  labs(x = "treatment", y = "rhizoid area [mm²]") +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 0) +
  xlab(NULL) +
  theme(text = element_text(size = 10))  # increase text size for all text elements

# define y-positions for labels: extract y-value of the top of the whisker of each boxplot from the ggplot object
letters$ypos <- layer_data(p)$ymax_final
# get y-axis range
yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
# add 5% of the y axis range to have the label slightly higher than the top of the whisker
letters$ypos <- letters$ypos + 0.1 * yrange

# show minimal boxplots with significance letters
p + geom_text(data = letters, aes(x = group, y = ypos, label = cld))


# export plot for Supp. Fig. 11B
ggsave("R_plots/polymorpha_rhizoid_area.png", width=1000, height=1000, unit="px")





## rhizoid/thallus ratio ----

# test assumption for ANOVA: equal variances (homoscedasticity) - Levene's test
leveneTest((rhizoid_thallus_ratio) ~ treatment, data = data) # significant heteroscedasticity

# we will use a non-parametric test for the rhizoid/thallus ratio

# Kruskal-Wallis test
kw <- kruskal.test(rhizoid_thallus_ratio ~ treatment, data = data)
# check p-value of Kruskal-Wallis
kw$p.value # p < 0.05

# Dunn's test
dunn <- dunn.test(data$rhizoid_thallus_ratio, data$treatment, method="bh", altp=TRUE) # Benjamini-Hochberg FDR correction. altp=TRUE makes sure we get p-values for test rejection at p<alpha and not p<alpha/2 (see documentation)

# put p-values in a named vector for multcompLetters()
pvec <- dunn$altP.adjusted
names(pvec) <- dunn$comparisons
# remove spaces around dash because it messes up multcompLetters()
names(pvec) <- gsub(" - ", "-", names(pvec))


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


# boxplot
p <- ggplot(data = data, aes(x = treatment, y = rhizoid_thallus_ratio)) +
  labs(x = "treatment", y = "rhizoid/thallus area ratio [mm²/mm²]") +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.2) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove gridlines
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 0) +
  xlab(NULL) +
  theme(text = element_text(size = 12)) +  # increase text size for all text elements
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # angle x-axis labels
  theme(axis.title.y = element_text(size=10))  # change y-axis label size

# define y-positions for labels: extract y-value of the top of the whisker of each boxplot from the ggplot object
letters$ypos <- layer_data(p)$ymax_final
# get y-axis range
yrange <- layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]
# add 5% of the y axis range to have the label slightly higher than the top of the whisker
letters$ypos <- letters$ypos + 0.1 * yrange

# show minimal boxplots with significance letters
p + geom_text(data = letters, aes(x = group, y = ypos, label = cld))


# export plot for Fig. 4E
ggsave("R_plots/polymorpha_rhizoid_thallus_ratio.png", width=900, height=900, unit="px")

