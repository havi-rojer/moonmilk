Hello welcome to our analysis folder for the 16S Amplicon Sequencing Project. This work was done by Havi Rojer and Sophia Aredas. Some of the files have "SA" in the name but please ignore that, that is a previous naming system but the work reflects the collaboration between Havi and Sophia. Here is an explanation of what the different folders include:

- 01_DADA2_MoonMilk_SA_files: 
Working with DADA2 to process our sample reads""
  1. Uploading our raw fastq files to generate quality plots and assess the quality of our reads
  2. Filtering and trimming bad sequences and bases from our sequencing files.
  3. Our dataset required that we needed to use the function in DADA2 "tryRC= TRUE" because reverse compliments (RC) were found in the dataset which ultimately skewed our results leading to errors. Once tryRC=TRUE was implemented the dataset produced a lot more normal data
  4. Infer errors on forward and reverse reads individually
  5. Identified ASVs on forward and reverse reads separately using the error model
  6. Merge ASVs into contiguous sequences
  7. Generate the ASV count table which will be used for downstream analysis with phyloseq
  8. Remove chimeras
  9. Assign taxonomy
  10. Track read counts
  
- 02_PreProcessing_SA: 
We processed the data into a phyloseq object (raw_preprocessed_physeq.RData)

We processed the data by:
1. removing chloroplasts, mitochondria, lowly abundant sample
2. writing the data file of the phyloseq output and saving it as (raw_preprocessed_physeq.RData)


- 03_Phylogenetics is broken down into 2 RMarkdown files:

03a_Phylogenetic_Tree
The main goal here is to make a phylogenetic tree with FastTree (ML tree)
1. Using our preprocessed phyloseq object we created ASV fasta file from phyloseq object
2. We then aligned the 16S sequences from the fasta file with MAFT
3. Made phylogenetic tree with FastTree

03b_Phylogenetic_Tree
1. We visualized the tree produced by FastTree
2. Added fasttree to phyloseq object
3. We created both unrooted and midrooted trees
4. We then pruned lowly abundant and NA taxa

04_Biodiversity
1. Calculate the Hill Diversity of the samples. 
2. Evaluate the rarefaction curves. 
3. Evaluate the Diversity values. 
4. Makes notes of specific samples and their seq depth. 

05_CommunityAnalysis
1. We made PCoA ordination plots for sorenson, bray-curtis, and weighted unifrac
2. We also made NMDS ordination plots for sorenson, bray-curtis, and weighted unifrac
3. We also looked at phylogenetic abundance for the moonmilk samples 

