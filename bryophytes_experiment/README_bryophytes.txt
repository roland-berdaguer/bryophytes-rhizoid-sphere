For the bryophytes experiment, NovaSeq reads are first processed with DADA2 using a modified "Ernakovich lab pipeline" (https://github.com/ErnakovichLab/dada2_ernakovichlab).
This results in the files "seqtab_final.txt", "taxonomy.tsv" and "repset.fasta".
These files are then loaded in the script "2_bryophytes_processing.R", which filters the data, and then exports the ASV sequences as "refSeq_Bac.fna" and "refSeq_Fun.fna".
Trees are then built from the exported sequences using the script "3_bryophytes_processing_qiime.sh". In order to run this script, the qiime2 conda environment needs to first be activated using "conda activate qiime2-2022.2". The resultilng trees are included in the "data_bryophytes" folder for your convenience.
The processed data is then analyzed using the script "4_bryophytes_analysis.R".
