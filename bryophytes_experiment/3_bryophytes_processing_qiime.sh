## treebuilding with MAFFT in qiime2 on server, then export in Newick format:
#
# first: copy refseq files to server


## 16S

conda activate qiime2-2022.2

qiime tools import \
--input-path refseq_Bac.fna \
--output-path refseq_Bac.qza \
--type 'FeatureData[Sequence]'

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences refseq_Bac.qza \
--output-dir mafft-fasttree-output-16S

qiime tools export \
--input-path mafft-fasttree-output-16S/rooted_tree.qza \
--output-path exported_rooted_tree_Bac

conda deactivate

mv exported_rooted_tree_Bac/tree.nwk ./tree_Bac.nwk
rm -r exported_rooted_tree_Bac

# copy newick tree to ./data_bryophytes/16S/





## ITS

conda activate qiime2-2022.2

qiime tools import \
--input-path refseq_Fun.fna \
--output-path refseq_Fun.qza \
--type 'FeatureData[Sequence]'

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences refseq_Fun.qza \
--output-dir mafft-fasttree-output-ITS

qiime tools export \
--input-path mafft-fasttree-output-ITS/rooted_tree.qza \
--output-path exported_rooted_tree_Fun

conda deactivate

mv exported_rooted_tree_Fun/tree.nwk ./tree_Fun.nwk
rm -r exported_rooted_tree_Fun

# copy newick tree to ./data_bryophytes/ITS/




