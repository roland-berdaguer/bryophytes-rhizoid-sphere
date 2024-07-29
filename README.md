This repository contains all scripts related to the microbiome analysis for the manuscript "The bryophyte rhizoid-sphere microbiome responds to water deficit".

This includes two separate experiments:
- the "bryophytes experiment": 16S and ITS sequencing of the rhizoid-sphere of Marchantia polymorpha, Marchantia paleacea and Physcomitrium patens
- the "angiosperms experiment": 16S sequencing of the rhizosphere of Solanum lycopersicum and Arabidopsis thaliana, and unplanted soil

For both experiments, raw reads are available on NCBI under Bioprojects PRJNA1049656 (bryophytes 16S), PRJNA1049759 (bryophytes ITS) and PRJNA1049762 (angiosperms 16S).
Raw reads are first processed using DADA2, leading to ASV tables, which are then filtered and analyzed in R.

Note: for all scripts, file paths and working directories need to be adapted to your own folder structure.
