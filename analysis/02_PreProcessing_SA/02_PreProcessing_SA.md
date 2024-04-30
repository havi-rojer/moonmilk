---
title: "02_PreProcessing"
author: "sophiaaredas"
date: "2024-04-30"
output:
  html_document:
    code_folding: show
    highlight: default
    keep_md: yes
    theme: journal
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 3
  pdf_document:
    toc: yes
editor_options:
  chunk_output_type: console
---



## Load libraries 

```r
#devtools::install_github("joey711/phyloseq")
pacman::p_load(devtools, phyloseq, tidyverse, dada2, install = FALSE)
```

## Goals

Here, we will process the data into a phyloseq object. 

- ASV table 
- Taxonomy Table 
- Track Reads (metadata)

Then, we will remove the following: 

1. Remove Chloroplasts
2. Remove Mitochondria. 
3. Remove samples without "enough" reads. 

Finally, write data file of phyloseq output.


## Load Data 

### ASV Table


```r
# First, load asv table
load("/local/workdir/sna49/moon_milk/moonmilk/data/01_DADA2/ASV_counts.RData")

# Inspect asv_tab
head(asv_tab)[,1:5]
```

```
##       ERR11588428_R1_filtered.fastq.gz ERR11588429_R1_filtered.fastq.gz
## ASV_1                                0                                0
## ASV_2                                0                                0
## ASV_3                                0                                0
## ASV_4                                0                                0
## ASV_5                                0                                0
## ASV_6                                0                                0
##       ERR11588430_R1_filtered.fastq.gz ERR11588431_R1_filtered.fastq.gz
## ASV_1                                0                             2261
## ASV_2                                0                             1807
## ASV_3                                0                             1593
## ASV_4                                0                             1396
## ASV_5                                0                             1191
## ASV_6                                0                             1232
##       ERR11588432_R1_filtered.fastq.gz
## ASV_1                             2077
## ASV_2                             1698
## ASV_3                             1527
## ASV_4                             1293
## ASV_5                             1140
## ASV_6                             1100
```

```r
# Fix names 
sample_titles <- colnames(asv_tab)
samples_fixed <- sapply(strsplit(basename(sample_titles), "_"), `[`,1) 
head(samples_fixed)
```

```
## [1] "ERR11588428" "ERR11588429" "ERR11588430" "ERR11588431" "ERR11588432"
## [6] "ERR11588433"
```

```r
# re-write the ASV count file to fix names 
colnames(asv_tab) <- samples_fixed
str(asv_tab)
```

```
##  int [1:4534, 1:10] 0 0 0 0 0 0 0 0 0 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:4534] "ASV_1" "ASV_2" "ASV_3" "ASV_4" ...
##   ..$ : chr [1:10] "ERR11588428" "ERR11588429" "ERR11588430" "ERR11588431" ...
```

###  Taxonomy Table

```r
tax_df <- read.table("/local/workdir/sna49/moon_milk/moonmilk/data/01_DADA2/ASV_taxonomy.tsv", sep = "\t", skip = 1)
head(tax_df)
```

```
##      V1       V2    V3   V4   V5   V6   V7   V8    V9
## 1 ASV_1 Bacteria GAL15 <NA> <NA> <NA> <NA> <NA> ASV_1
## 2 ASV_2 Bacteria GAL15 <NA> <NA> <NA> <NA> <NA> ASV_2
## 3 ASV_3 Bacteria GAL15 <NA> <NA> <NA> <NA> <NA> ASV_3
## 4 ASV_4 Bacteria GAL15 <NA> <NA> <NA> <NA> <NA> ASV_4
## 5 ASV_5 Bacteria GAL15 <NA> <NA> <NA> <NA> <NA> ASV_5
## 6 ASV_6 Bacteria GAL15 <NA> <NA> <NA> <NA> <NA> ASV_6
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                 V10
## 1 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCAGTAGTC
## 2 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCGAGTAGTC
## 3 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTC
## 4 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCTAGTAGTC
## 5 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCTGTAGTC
## 6 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCGGGTAGTC
```

```r
# fix column names 
colnames(tax_df) <- c("asv_names", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV", "ASVseq")

head(tax_df)
```

```
##   asv_names  Kingdom Phylum Class Order Family Genus Species   ASV
## 1     ASV_1 Bacteria  GAL15  <NA>  <NA>   <NA>  <NA>    <NA> ASV_1
## 2     ASV_2 Bacteria  GAL15  <NA>  <NA>   <NA>  <NA>    <NA> ASV_2
## 3     ASV_3 Bacteria  GAL15  <NA>  <NA>   <NA>  <NA>    <NA> ASV_3
## 4     ASV_4 Bacteria  GAL15  <NA>  <NA>   <NA>  <NA>    <NA> ASV_4
## 5     ASV_5 Bacteria  GAL15  <NA>  <NA>   <NA>  <NA>    <NA> ASV_5
## 6     ASV_6 Bacteria  GAL15  <NA>  <NA>   <NA>  <NA>    <NA> ASV_6
##                                                                                                                                                                                                                                                                                                                                                                                                                                                              ASVseq
## 1 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCAGTAGTC
## 2 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCGAGTAGTC
## 3 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTC
## 4 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCTAGTAGTC
## 5 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCTGTAGTC
## 6 CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCGGGTAGTC
```

```r
# Taxonomy Table Matrix
tax_mat <- 
  tax_df %>%
  tibble::column_to_rownames(., var = "asv_names") %>%
  as.matrix()
```

###  Track Reads Data

```r
load("/local/workdir/sna49/moon_milk/moonmilk/data/01_DADA2/track_read_counts.RData")

# Take a look at the data
head(track_counts_df)
```

```
##         names  input filtered denoisedF denoisedR merged nochim
## 1 ERR11588428  96845    57420     53646     52928  39698  21236
## 2 ERR11588429 101519    57834     52927     52528  37117  20768
## 3 ERR11588430  95676    56552     51951     51872  37096  20697
## 4 ERR11588431  90880    51923     51671     51312  50857  29177
## 5 ERR11588432  85608    50841     50594     50229  49720  28614
## 6 ERR11588433  73835    42323     42137     41784  41173  24959
##   perc_reads_retained
## 1            21.92782
## 2            20.45725
## 3            21.63238
## 4            32.10497
## 5            33.42445
## 6            33.80375
```

```r
dim(track_counts_df)
```

```
## [1] 10  8
```

```r
# Load in metadata
metadata_df <- read.csv("/local/workdir/sna49/moon_milk/moonmilk/data/02_PreProcessing/metadata.csv")
dim(metadata_df)
```

```
## [1] 10  2
```

```r
colnames(metadata_df)
```

```
## [1] "run_accession" "sample_title"
```

```r
# Replace the name of the column
colnames(metadata_df)[1] <- "names"
print(metadata_df)
```

```
##          names          sample_title
## 1  ERR11588428 Fata Apei Cave - FA_1
## 2  ERR11588429 Fata Apei Cave - FA_2
## 3  ERR11588430 Fata Apei Cave - FA_3
## 4  ERR11588431   Ferice Cave - PF5_1
## 5  ERR11588432   Ferice Cave - PF5_2
## 6  ERR11588433   Ferice Cave - PF5_3
## 7  ERR11588434   Ferice Cave - PF2_1
## 8  ERR11588435   Ferice Cave - PF2_2
## 9  ERR11588436      Nestor Cave - NE
## 10 ERR11588437    Tausoare Cave - TS
```

```r
# Add column with cave_name
metadata_df <- separate(metadata_df, sample_title, into = c("cave", "sample"), sep = "-")  

print(metadata_df)
```

```
##          names            cave sample
## 1  ERR11588428 Fata Apei Cave    FA_1
## 2  ERR11588429 Fata Apei Cave    FA_2
## 3  ERR11588430 Fata Apei Cave    FA_3
## 4  ERR11588431    Ferice Cave   PF5_1
## 5  ERR11588432    Ferice Cave   PF5_2
## 6  ERR11588433    Ferice Cave   PF5_3
## 7  ERR11588434    Ferice Cave   PF2_1
## 8  ERR11588435    Ferice Cave   PF2_2
## 9  ERR11588436    Nestor Cave      NE
## 10 ERR11588437  Tausoare Cave      TS
```

```r
head(track_counts_df)
```

```
##         names  input filtered denoisedF denoisedR merged nochim
## 1 ERR11588428  96845    57420     53646     52928  39698  21236
## 2 ERR11588429 101519    57834     52927     52528  37117  20768
## 3 ERR11588430  95676    56552     51951     51872  37096  20697
## 4 ERR11588431  90880    51923     51671     51312  50857  29177
## 5 ERR11588432  85608    50841     50594     50229  49720  28614
## 6 ERR11588433  73835    42323     42137     41784  41173  24959
##   perc_reads_retained
## 1            21.92782
## 2            20.45725
## 3            21.63238
## 4            32.10497
## 5            33.42445
## 6            33.80375
```

```r
metadata_track_reads_df <- 
  track_counts_df %>%
  left_join(., metadata_df, by = "names")

# Intuition check 
head(metadata_track_reads_df)
```

```
##         names  input filtered denoisedF denoisedR merged nochim
## 1 ERR11588428  96845    57420     53646     52928  39698  21236
## 2 ERR11588429 101519    57834     52927     52528  37117  20768
## 3 ERR11588430  95676    56552     51951     51872  37096  20697
## 4 ERR11588431  90880    51923     51671     51312  50857  29177
## 5 ERR11588432  85608    50841     50594     50229  49720  28614
## 6 ERR11588433  73835    42323     42137     41784  41173  24959
##   perc_reads_retained            cave sample
## 1            21.92782 Fata Apei Cave    FA_1
## 2            20.45725 Fata Apei Cave    FA_2
## 3            21.63238 Fata Apei Cave    FA_3
## 4            32.10497    Ferice Cave   PF5_1
## 5            33.42445    Ferice Cave   PF5_2
## 6            33.80375    Ferice Cave   PF5_3
```

```r
# Update row.names to be sample names 
## Before 
row.names(metadata_track_reads_df)
```

```
##  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"
```

```r
# Rewrite 
row.names(metadata_track_reads_df) <- metadata_track_reads_df$names
# Check afterwards that it worked 
row.names(metadata_track_reads_df)
```

```
##  [1] "ERR11588428" "ERR11588429" "ERR11588430" "ERR11588431" "ERR11588432"
##  [6] "ERR11588433" "ERR11588434" "ERR11588435" "ERR11588436" "ERR11588437"
```

```r
# intuition check
head(metadata_track_reads_df)
```

```
##                   names  input filtered denoisedF denoisedR merged nochim
## ERR11588428 ERR11588428  96845    57420     53646     52928  39698  21236
## ERR11588429 ERR11588429 101519    57834     52927     52528  37117  20768
## ERR11588430 ERR11588430  95676    56552     51951     51872  37096  20697
## ERR11588431 ERR11588431  90880    51923     51671     51312  50857  29177
## ERR11588432 ERR11588432  85608    50841     50594     50229  49720  28614
## ERR11588433 ERR11588433  73835    42323     42137     41784  41173  24959
##             perc_reads_retained            cave sample
## ERR11588428            21.92782 Fata Apei Cave    FA_1
## ERR11588429            20.45725 Fata Apei Cave    FA_2
## ERR11588430            21.63238 Fata Apei Cave    FA_3
## ERR11588431            32.10497    Ferice Cave   PF5_1
## ERR11588432            33.42445    Ferice Cave   PF5_2
## ERR11588433            33.80375    Ferice Cave   PF5_3
```

## Handoff to phyloseq

```r
# double check it's all good 
dim(asv_tab)
```

```
## [1] 4534   10
```

```r
dim(tax_mat)
```

```
## [1] 4534    9
```

```r
# Intuition check 
stopifnot(row.names(asv_tab) == row.names(tax_mat))

# Construct the phyloseq object 
raw_physeq <- phyloseq(otu_table(asv_tab, taxa_are_rows = TRUE),
                       sample_data(metadata_track_reads_df), tax_table(tax_mat))
                       
raw_physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4534 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 4534 taxa by 9 taxonomic ranks ]
```

```r
# Save this raw phyloseq object 
save(raw_physeq, file = "/local/workdir/sna49/moon_milk/moonmilk/data/02_PreProcessing/raw_physeq.RData")
```

## Clean up the data

Remove: 

1. Chloroplasts
2. mitochondria  


```r
# Remind myself of tax table 
#View(tax_mat)

# Make new physeq without chloroplasts
noChloros_physeq <- 
  raw_physeq %>% 
  # rm chloroplasts
  subset_taxa(Order != "Chloroplast" | is.na(Order))
  
# How many taxa were chloroplasts? 
num_chloro_ASVs <- ntaxa(raw_physeq) - ntaxa(noChloros_physeq)
num_chloro_ASVs
```

```
## [1] 3
```

```r
# Intuition chek 
#noChloros_physeq %>%
#  tax_table() %>%
#  data.frame() %>%
#  View()

# remove mitochondria 
noChlorosMitos_physeq <- 
  noChloros_physeq %>%
  subset_taxa(Family != "Mitochondria" | is.na(Family))

# How many mitochondrial ASVs? 
num_mito_ASVs <- ntaxa(noChloros_physeq) - ntaxa(noChlorosMitos_physeq)
num_mito_ASVs
```

```
## [1] 0
```

```r
noChlorosMitos_physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4531 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 4531 taxa by 9 taxonomic ranks ]
```

```r
# How many total asvs were removed from chloros and mitos 
ntaxa(raw_physeq) - ntaxa(noChlorosMitos_physeq)
```

```
## [1] 3
```

```r
# proportion of asvs kept? 
ntaxa(noChlorosMitos_physeq)/ntaxa(raw_physeq)
```

```
## [1] 0.9993383
```

# Evaulate and remove the control samples 

We do not have any control samples or a mock community.

# Evaluate the Sequencing Depth 


```r
# The current data object
noChlorosMitos_physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4531 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 4531 taxa by 9 taxonomic ranks ]
```

```r
# What is the library size/sequencing depth for each sample? 
seqSums_df <- 
  noChlorosMitos_physeq %>%
  otu_table() %>%
  # Sum each sample column 
  colSums() %>%
  data.frame() %>%
  rownames_to_column(var = "names") %>%
  left_join(., metadata_track_reads_df, by = "names") 

# Rename second column 
colnames(seqSums_df)[2] <- "TotalSeqs"

# check
dim(seqSums_df)
```

```
## [1] 10 11
```

```r
head(seqSums_df)
```

```
##         names TotalSeqs  input filtered denoisedF denoisedR merged nochim
## 1 ERR11588428     21236  96845    57420     53646     52928  39698  21236
## 2 ERR11588429     20768 101519    57834     52927     52528  37117  20768
## 3 ERR11588430     20697  95676    56552     51951     51872  37096  20697
## 4 ERR11588431     29177  90880    51923     51671     51312  50857  29177
## 5 ERR11588432     28614  85608    50841     50594     50229  49720  28614
## 6 ERR11588433     24959  73835    42323     42137     41784  41173  24959
##   perc_reads_retained            cave sample
## 1            21.92782 Fata Apei Cave    FA_1
## 2            20.45725 Fata Apei Cave    FA_2
## 3            21.63238 Fata Apei Cave    FA_3
## 4            32.10497    Ferice Cave   PF5_1
## 5            33.42445    Ferice Cave   PF5_2
## 6            33.80375    Ferice Cave   PF5_3
```

```r
# Show the depth of samples 
seqSums_df %>%
  dplyr::select(names, TotalSeqs) %>%
  arrange(TotalSeqs) %>%
  head()
```

```
##         names TotalSeqs
## 1 ERR11588436     13720
## 2 ERR11588434     16019
## 3 ERR11588435     17222
## 4 ERR11588437     18630
## 5 ERR11588430     20697
## 6 ERR11588429     20768
```

```r
# plot it! 
seqSums_df %>%
  ggplot(aes(x=reorder(names, TotalSeqs), y = TotalSeqs, fill = cave)) + 
  geom_bar(stat = "identity") 
```

<img src="/local/workdir/sna49/moon_milk/moonmilk/figures/02_PreProcessing_SAseq-depth-1.png" style="display: block; margin: auto;" />

```r
# Density plot

seqSums_df %>%
  ggplot(aes(TotalSeqs, fill = cave)) +
  geom_density(alpha = 0.5)
```

```
## Warning: Groups with fewer than two data points have been dropped.
## Groups with fewer than two data points have been dropped.
```

```
## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
## -Inf

## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
## -Inf
```

<img src="/local/workdir/sna49/moon_milk/moonmilk/figures/02_PreProcessing_SAseq-depth-2.png" style="display: block; margin: auto;" />

# Remove samples with few reads 


```r
# What's the min seq depth? 
min(sample_sums(noChlorosMitos_physeq))
```

```
## [1] 13720
```

```r
# We don't have any samples with zero seq depth so we're good

raw_preprocessed_physeq <- noChlorosMitos_physeq
```

# Save Preprocessed Phyloseq Object

```r
#set your WD to wherever you wanna put this file

save(raw_preprocessed_physeq, file = "./raw_preprocessed_physeq.RData")
```

## Session Information 

```r
# Ensure reproducibility 
devtools::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.3.2 (2023-10-31)
##  os       Rocky Linux 9.0 (Blue Onyx)
##  system   x86_64, linux-gnu
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/New_York
##  date     2024-04-30
##  pandoc   3.1.1 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package              * version   date (UTC) lib source
##  abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.2)
##  ade4                   1.7-22    2023-02-06 [1] CRAN (R 4.3.2)
##  ape                    5.7-1     2023-03-13 [2] CRAN (R 4.3.2)
##  Biobase                2.62.0    2023-10-24 [2] Bioconductor
##  BiocGenerics           0.48.1    2023-11-01 [2] Bioconductor
##  BiocParallel           1.36.0    2023-10-24 [2] Bioconductor
##  biomformat             1.30.0    2023-10-24 [1] Bioconductor
##  Biostrings             2.70.1    2023-10-25 [2] Bioconductor
##  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.2)
##  bslib                  0.5.1     2023-08-11 [2] CRAN (R 4.3.2)
##  cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.3.2)
##  callr                  3.7.3     2022-11-02 [2] CRAN (R 4.3.2)
##  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.2)
##  cluster                2.1.4     2022-08-22 [2] CRAN (R 4.3.2)
##  codetools              0.2-19    2023-02-01 [2] CRAN (R 4.3.2)
##  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.2)
##  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.2)
##  dada2                * 1.30.0    2023-10-24 [1] Bioconductor
##  data.table             1.14.8    2023-02-17 [2] CRAN (R 4.3.2)
##  DelayedArray           0.28.0    2023-10-24 [2] Bioconductor
##  deldir                 1.0-9     2023-05-17 [2] CRAN (R 4.3.2)
##  devtools             * 2.4.4     2022-07-20 [2] CRAN (R 4.2.1)
##  digest                 0.6.33    2023-07-07 [2] CRAN (R 4.3.2)
##  dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.2)
##  ellipsis               0.3.2     2021-04-29 [2] CRAN (R 4.3.2)
##  evaluate               0.23      2023-11-01 [2] CRAN (R 4.3.2)
##  fansi                  1.0.5     2023-10-08 [2] CRAN (R 4.3.2)
##  farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.2)
##  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.2)
##  forcats              * 1.0.0     2023-01-29 [1] CRAN (R 4.3.2)
##  foreach                1.5.2     2022-02-02 [2] CRAN (R 4.3.2)
##  fs                     1.6.3     2023-07-20 [2] CRAN (R 4.3.2)
##  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.2)
##  GenomeInfoDb           1.38.0    2023-10-24 [2] Bioconductor
##  GenomeInfoDbData       1.2.11    2023-11-07 [2] Bioconductor
##  GenomicAlignments      1.38.0    2023-10-24 [2] Bioconductor
##  GenomicRanges          1.54.1    2023-10-29 [2] Bioconductor
##  ggplot2              * 3.5.0     2024-02-23 [2] CRAN (R 4.3.2)
##  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.2)
##  gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.2)
##  highr                  0.10      2022-12-22 [2] CRAN (R 4.3.2)
##  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.2)
##  htmltools              0.5.7     2023-11-03 [2] CRAN (R 4.3.2)
##  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.3.2)
##  httpuv                 1.6.12    2023-10-23 [2] CRAN (R 4.3.2)
##  hwriter                1.3.2.1   2022-04-08 [1] CRAN (R 4.3.2)
##  igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.2)
##  interp                 1.1-6     2024-01-26 [1] CRAN (R 4.3.2)
##  IRanges                2.36.0    2023-10-24 [2] Bioconductor
##  iterators              1.0.14    2022-02-05 [2] CRAN (R 4.3.2)
##  jpeg                   0.1-10    2022-11-29 [1] CRAN (R 4.3.2)
##  jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.3.2)
##  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.2)
##  knitr                  1.45      2023-10-30 [2] CRAN (R 4.3.2)
##  labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.2)
##  later                  1.3.1     2023-05-02 [2] CRAN (R 4.3.2)
##  lattice                0.21-9    2023-10-01 [2] CRAN (R 4.3.2)
##  latticeExtra           0.6-30    2022-07-04 [1] CRAN (R 4.3.2)
##  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.2)
##  lubridate            * 1.9.3     2023-09-27 [1] CRAN (R 4.3.2)
##  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.2)
##  MASS                   7.3-60    2023-05-04 [2] CRAN (R 4.3.2)
##  Matrix                 1.6-1.1   2023-09-18 [2] CRAN (R 4.3.2)
##  MatrixGenerics         1.14.0    2023-10-24 [2] Bioconductor
##  matrixStats            1.1.0     2023-11-07 [2] CRAN (R 4.3.2)
##  memoise                2.0.1     2021-11-26 [2] CRAN (R 4.3.2)
##  mgcv                   1.9-0     2023-07-11 [2] CRAN (R 4.3.2)
##  mime                   0.12      2021-09-28 [2] CRAN (R 4.3.2)
##  miniUI                 0.1.1.1   2018-05-18 [2] CRAN (R 4.3.2)
##  multtest               2.58.0    2023-10-24 [1] Bioconductor
##  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.2)
##  nlme                   3.1-163   2023-08-09 [2] CRAN (R 4.3.2)
##  pacman                 0.5.1     2019-03-11 [1] CRAN (R 4.3.2)
##  permute                0.9-7     2022-01-27 [1] CRAN (R 4.3.2)
##  phyloseq             * 1.46.0    2023-10-24 [1] Bioconductor
##  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.2)
##  pkgbuild               1.4.2     2023-06-26 [2] CRAN (R 4.3.2)
##  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.2)
##  pkgload                1.3.3     2023-09-22 [2] CRAN (R 4.3.2)
##  plyr                   1.8.9     2023-10-02 [2] CRAN (R 4.3.2)
##  png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.2)
##  prettyunits            1.2.0     2023-09-24 [2] CRAN (R 4.3.2)
##  processx               3.8.2     2023-06-30 [2] CRAN (R 4.3.2)
##  profvis                0.3.8     2023-05-02 [2] CRAN (R 4.3.2)
##  promises               1.2.1     2023-08-10 [2] CRAN (R 4.3.2)
##  ps                     1.7.5     2023-04-18 [2] CRAN (R 4.3.2)
##  purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.3.2)
##  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.2)
##  RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.2)
##  Rcpp                 * 1.0.11    2023-07-06 [2] CRAN (R 4.3.2)
##  RcppParallel           5.1.7     2023-02-27 [2] CRAN (R 4.3.2)
##  RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.2)
##  readr                * 2.1.5     2024-01-10 [1] CRAN (R 4.3.2)
##  remotes                2.4.2.1   2023-07-18 [2] CRAN (R 4.3.2)
##  reshape2               1.4.4     2020-04-09 [2] CRAN (R 4.3.2)
##  rhdf5                  2.46.1    2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
##  rhdf5filters           1.14.1    2023-11-06 [1] Bioconductor
##  Rhdf5lib               1.24.2    2024-02-07 [1] Bioconductor 3.18 (R 4.3.2)
##  rlang                  1.1.2     2023-11-04 [2] CRAN (R 4.3.2)
##  rmarkdown              2.25      2023-09-18 [2] CRAN (R 4.3.2)
##  Rsamtools              2.18.0    2023-10-24 [2] Bioconductor
##  rstudioapi             0.15.0    2023-07-07 [2] CRAN (R 4.3.2)
##  S4Arrays               1.2.0     2023-10-24 [2] Bioconductor
##  S4Vectors              0.40.1    2023-10-26 [2] Bioconductor
##  sass                   0.4.7     2023-07-15 [2] CRAN (R 4.3.2)
##  scales                 1.3.0     2023-11-28 [2] CRAN (R 4.3.2)
##  sessioninfo            1.2.2     2021-12-06 [2] CRAN (R 4.3.2)
##  shiny                  1.7.5.1   2023-10-14 [2] CRAN (R 4.3.2)
##  ShortRead              1.60.0    2023-10-24 [1] Bioconductor
##  SparseArray            1.2.1     2023-11-05 [2] Bioconductor
##  stringi                1.7.12    2023-01-11 [2] CRAN (R 4.3.2)
##  stringr              * 1.5.0     2022-12-02 [2] CRAN (R 4.3.2)
##  SummarizedExperiment   1.32.0    2023-10-24 [2] Bioconductor
##  survival               3.5-7     2023-08-14 [2] CRAN (R 4.3.2)
##  tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.3.2)
##  tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.3.2)
##  tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.3.2)
##  tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.3.2)
##  timechange             0.3.0     2024-01-18 [1] CRAN (R 4.3.2)
##  tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.2)
##  urlchecker             1.0.1     2021-11-30 [2] CRAN (R 4.3.2)
##  usethis              * 2.2.2     2023-07-06 [2] CRAN (R 4.3.2)
##  utf8                   1.2.4     2023-10-22 [2] CRAN (R 4.3.2)
##  vctrs                  0.6.4     2023-10-12 [2] CRAN (R 4.3.2)
##  vegan                  2.6-4     2022-10-11 [1] CRAN (R 4.3.2)
##  withr                  2.5.2     2023-10-30 [2] CRAN (R 4.3.2)
##  xfun                   0.41      2023-11-01 [2] CRAN (R 4.3.2)
##  xtable                 1.8-4     2019-04-21 [2] CRAN (R 4.3.2)
##  XVector                0.42.0    2023-10-24 [2] Bioconductor
##  yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.2)
##  zlibbioc               1.48.0    2023-10-24 [2] Bioconductor
## 
##  [1] /home/sna49/R/x86_64-pc-linux-gnu-library/4.3
##  [2] /programs/R-4.3.2/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```
