# Qiime version: qiime2-2022.2


# IMPORT DATA

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path angiosperms_manifest.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path angiosperms_preDada2_withPrimers.qza


# create visualization of demultiplexed files
qiime demux summarize \
--i-data angiosperms_preDada2_withPrimers.qza \
--o-visualization angiosperms_preDada2_withPrimers.qzv


# REMOVE PRIMERS

# run cutadapt plugin of qiime, removing primer sequences and discarding reads without primers
qiime cutadapt trim-paired \
--p-cores 10 \
--i-demultiplexed-sequences angiosperms_preDada2_withPrimers.qza \
--p-front-f CCTACGGGNGGCWGCAG \
--p-front-r GGACTACHVGGGTATCTAATCC \
--p-error-rate 0 \
--p-no-indels \
--p-discard-untrimmed \
--o-trimmed-sequences angiosperms_preDada2_TrimmedPrimers.qza \
--verbose

# create visualization of artefact
qiime demux summarize \
--i-data  angiosperms_preDada2_TrimmedPrimers.qza \
--o-visualization  angiosperms_preDada2_TrimmedPrimers.qzv

# use dada2 to trim low-quality ends, remove chimeric sequences, remove singletons, join denoised paired-end reads, and then dereplicate into ASVs. trim size based on quality plot.
qiime dada2 denoise-paired \
--i-demultiplexed-seqs angiosperms_preDada2_TrimmedPrimers.qza \
--p-trunc-len-f 280 \
--p-trunc-len-r 200 \
--p-max-ee-f 1 \
--p-max-ee-r 1 \
--p-n-threads 0 \
--o-table feature-table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza

# visualize artefact, feature table
qiime feature-table summarize \
--i-table feature-table.qza \
--o-visualization feature-table.qzv \
--m-sample-metadata-file angiosperms_mapping.txt

# visualize artefact, representative sequences table
qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

# visualize artefact, denoising stats
qiime metadata tabulate \
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv

# assign taxonomy ranks on representative sequences, based on the trained SILVA reference
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138.1-ssu-nr99-341f-805r-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# create a vizualization of the taxonomy
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv


# make phylogeny
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza




# export feature table in BIOM format
qiime tools export \
--input-path feature-table.qza \
--output-path angiosperms_qiime_output

# convert feature table from BIOM to tsv
biom convert -i angiosperms_qiime_output/feature-table.biom -o angiosperms_qiime_output/feature-table.tsv --to-tsv
# MANUALLY EDIT feature-table.tsv to remove the first line ("# Constructed from BIOM file") and beginning of second line: remove "#OTU_ID" and the tab following that

# export taxonomy
qiime tools export \
--input-path taxonomy.qza \
--output-path angiosperms_qiime_output

# export unrooted phylogenetic tree in Newick format
qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path angiosperms_qiime_output

# export representative sequences in FASTA format
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path angiosperms_qiime_output



