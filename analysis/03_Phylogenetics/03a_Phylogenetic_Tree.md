---
title: "03a_Phylogenetics"
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



---
# Goals of this file 

The goal is to create a phylogenetic tree! 

1. Load in preprocessed phyloseq object.  
2. Create ASV fasta file from the phyloseq object. 
3. Align the 16S sequences from fasta file with MAFFT. 
4. Create a tree with FastTree2.  

# Set the seed 

```r
set.seed(0909199)
```


## Load Packages & Phyloseq Object

```r
# phytools, ggtree, RColorBrewer
pacman::p_load(phytools, ggtree, RColorBrewer, install = FALSE)

# Load physeq 
load("/local/workdir/sna49/moon_milk/moonmilk/data/02_PreProcessing/raw_preprocessed_physeq.RData")
raw_preprocessed_physeq
```

```
## Loading required package: phyloseq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4531 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 4531 taxa by 9 taxonomic ranks ]
```

# Create Fasta File of ASV and their Sequences

This fasta file will be used to create our alignment in MAFFT


```r
# pull out ASV seqs and ASV names 
asv_seq_df <- 
  raw_preprocessed_physeq@tax_table %>%
  data.frame() %>%
  dplyr::select(ASV, ASVseq)

#View(asv_seq_df)

# Add the > to make fasta header 
asv_seq_df$ASV <- paste0(">",asv_seq_df$ASV)
#View(asv_seq_df)

# Create the fasta object 
asv_seq_fasta <- c(rbind(asv_seq_df$ASV, asv_seq_df$ASVseq))
head(asv_seq_fasta)
```

```
## [1] ">ASV_1"                                                                                                                                                                                                                                                                                                                                                                                                                                                           
## [2] "CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCAGTAGTC"
## [3] ">ASV_2"                                                                                                                                                                                                                                                                                                                                                                                                                                                           
## [4] "CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCGAGTAGTC"
## [5] ">ASV_3"                                                                                                                                                                                                                                                                                                                                                                                                                                                           
## [6] "CCTACGGGAGGCAGCAGTTGGGAATTTTGGGCAATGGGCGCAAGCCTGACCCAGCGACGCCGCGTGGGGGATGAAGGCCTTCGGGTCGTAAACCCCTGTCAGGAGGGAAGAAGCCCGCAAGGGTGACGGTACCTCCAGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGGGCGAGCGTTGTCCGGATTCACTGGGCGTAAAGCGCGAGTAGGTGGCCTGTTAAGTGGCGTGTGAAAGCCTGCGGCTCAACCGTAGGAGGTCGCGCCAGACTGGCGGGCTTGAGGGTGGCAGAGGGTGATGGAATTCCCGGTGTAGCGGTGAAATGCGCAGATATCGGGAGGAACGCCAGTGGCGAAGGCGGTCACCTGGGCCATCCCTAACGCTGAGTCGCGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTC"
```

```r
# Write to a file 
write(asv_seq_fasta, 
      file = "data/03_Phylogenetic_Tree/preprocessed_ASVs.fasta")
```

# Align the 16S sequences from fasta file with MAFFT

- `engine.opts = '-l'`: give us original terminal base environment.


```bash
# Write bash code to run mafft 
# First provide the path to mafft 
export PATH=/programs/mafft/bin:$PATH
# change directories to provide the fasta file we made above 
cd data/03_Phylogenetic_Tree/
pwd

# Set a seed  - using same seed as before for consistency 
RANDOM=1234

# Run Mafft 
# To test in the shell directly from Rmd 
# mac: command + option + enter 
# Windows: control + alt + enter 
# For now, use default options, note the version 
# MAFFT automatically knows that it's a nucleotide alignment 
/programs/mafft/bin/mafft --auto preprocessed_ASVs.fasta > MAFFT_aligned_ASVs.fasta

# Change back to the project directory 
cd ../../
pwd
```

```
## /local/workdir/sna49/moon_milk/moonmilk/data/03_Phylogenetic_Tree
## nthread = 0
## nthreadpair = 0
## nthreadtb = 0
## ppenalty_ex = 0
## stacksize: 8192 kb
## generating a scoring matrix for nucleotide (dist=200) ... done
## Gap Penalty = -1.53, +0.00, +0.00
## 
## 
## 
## Making a distance matrix ..
## 
    1 / 4531
  101 / 4531
  201 / 4531
  301 / 4531
  401 / 4531
  501 / 4531
  601 / 4531
  701 / 4531
  801 / 4531
  901 / 4531
 1001 / 4531
 1101 / 4531
 1201 / 4531
 1301 / 4531
 1401 / 4531
 1501 / 4531
 1601 / 4531
 1701 / 4531
 1801 / 4531
 1901 / 4531
 2001 / 4531
 2101 / 4531
 2201 / 4531
 2301 / 4531
 2401 / 4531
 2501 / 4531
 2601 / 4531
 2701 / 4531
 2801 / 4531
 2901 / 4531
 3001 / 4531
 3101 / 4531
 3201 / 4531
 3301 / 4531
 3401 / 4531
 3501 / 4531
 3601 / 4531
 3701 / 4531
 3801 / 4531
 3901 / 4531
 4001 / 4531
 4101 / 4531
 4201 / 4531
 4301 / 4531
 4401 / 4531
 4501 / 4531
## done.
## 
## Constructing a UPGMA tree (efffree=0) ... 
## 
    0 / 4531
   10 / 4531
   20 / 4531
   30 / 4531
   40 / 4531
   50 / 4531
   60 / 4531
   70 / 4531
   80 / 4531
   90 / 4531
  100 / 4531
  110 / 4531
  120 / 4531
  130 / 4531
  140 / 4531
  150 / 4531
  160 / 4531
  170 / 4531
  180 / 4531
  190 / 4531
  200 / 4531
  210 / 4531
  220 / 4531
  230 / 4531
  240 / 4531
  250 / 4531
  260 / 4531
  270 / 4531
  280 / 4531
  290 / 4531
  300 / 4531
  310 / 4531
  320 / 4531
  330 / 4531
  340 / 4531
  350 / 4531
  360 / 4531
  370 / 4531
  380 / 4531
  390 / 4531
  400 / 4531
  410 / 4531
  420 / 4531
  430 / 4531
  440 / 4531
  450 / 4531
  460 / 4531
  470 / 4531
  480 / 4531
  490 / 4531
  500 / 4531
  510 / 4531
  520 / 4531
  530 / 4531
  540 / 4531
  550 / 4531
  560 / 4531
  570 / 4531
  580 / 4531
  590 / 4531
  600 / 4531
  610 / 4531
  620 / 4531
  630 / 4531
  640 / 4531
  650 / 4531
  660 / 4531
  670 / 4531
  680 / 4531
  690 / 4531
  700 / 4531
  710 / 4531
  720 / 4531
  730 / 4531
  740 / 4531
  750 / 4531
  760 / 4531
  770 / 4531
  780 / 4531
  790 / 4531
  800 / 4531
  810 / 4531
  820 / 4531
  830 / 4531
  840 / 4531
  850 / 4531
  860 / 4531
  870 / 4531
  880 / 4531
  890 / 4531
  900 / 4531
  910 / 4531
  920 / 4531
  930 / 4531
  940 / 4531
  950 / 4531
  960 / 4531
  970 / 4531
  980 / 4531
  990 / 4531
 1000 / 4531
 1010 / 4531
 1020 / 4531
 1030 / 4531
 1040 / 4531
 1050 / 4531
 1060 / 4531
 1070 / 4531
 1080 / 4531
 1090 / 4531
 1100 / 4531
 1110 / 4531
 1120 / 4531
 1130 / 4531
 1140 / 4531
 1150 / 4531
 1160 / 4531
 1170 / 4531
 1180 / 4531
 1190 / 4531
 1200 / 4531
 1210 / 4531
 1220 / 4531
 1230 / 4531
 1240 / 4531
 1250 / 4531
 1260 / 4531
 1270 / 4531
 1280 / 4531
 1290 / 4531
 1300 / 4531
 1310 / 4531
 1320 / 4531
 1330 / 4531
 1340 / 4531
 1350 / 4531
 1360 / 4531
 1370 / 4531
 1380 / 4531
 1390 / 4531
 1400 / 4531
 1410 / 4531
 1420 / 4531
 1430 / 4531
 1440 / 4531
 1450 / 4531
 1460 / 4531
 1470 / 4531
 1480 / 4531
 1490 / 4531
 1500 / 4531
 1510 / 4531
 1520 / 4531
 1530 / 4531
 1540 / 4531
 1550 / 4531
 1560 / 4531
 1570 / 4531
 1580 / 4531
 1590 / 4531
 1600 / 4531
 1610 / 4531
 1620 / 4531
 1630 / 4531
 1640 / 4531
 1650 / 4531
 1660 / 4531
 1670 / 4531
 1680 / 4531
 1690 / 4531
 1700 / 4531
 1710 / 4531
 1720 / 4531
 1730 / 4531
 1740 / 4531
 1750 / 4531
 1760 / 4531
 1770 / 4531
 1780 / 4531
 1790 / 4531
 1800 / 4531
 1810 / 4531
 1820 / 4531
 1830 / 4531
 1840 / 4531
 1850 / 4531
 1860 / 4531
 1870 / 4531
 1880 / 4531
 1890 / 4531
 1900 / 4531
 1910 / 4531
 1920 / 4531
 1930 / 4531
 1940 / 4531
 1950 / 4531
 1960 / 4531
 1970 / 4531
 1980 / 4531
 1990 / 4531
 2000 / 4531
 2010 / 4531
 2020 / 4531
 2030 / 4531
 2040 / 4531
 2050 / 4531
 2060 / 4531
 2070 / 4531
 2080 / 4531
 2090 / 4531
 2100 / 4531
 2110 / 4531
 2120 / 4531
 2130 / 4531
 2140 / 4531
 2150 / 4531
 2160 / 4531
 2170 / 4531
 2180 / 4531
 2190 / 4531
 2200 / 4531
 2210 / 4531
 2220 / 4531
 2230 / 4531
 2240 / 4531
 2250 / 4531
 2260 / 4531
 2270 / 4531
 2280 / 4531
 2290 / 4531
 2300 / 4531
 2310 / 4531
 2320 / 4531
 2330 / 4531
 2340 / 4531
 2350 / 4531
 2360 / 4531
 2370 / 4531
 2380 / 4531
 2390 / 4531
 2400 / 4531
 2410 / 4531
 2420 / 4531
 2430 / 4531
 2440 / 4531
 2450 / 4531
 2460 / 4531
 2470 / 4531
 2480 / 4531
 2490 / 4531
 2500 / 4531
 2510 / 4531
 2520 / 4531
 2530 / 4531
 2540 / 4531
 2550 / 4531
 2560 / 4531
 2570 / 4531
 2580 / 4531
 2590 / 4531
 2600 / 4531
 2610 / 4531
 2620 / 4531
 2630 / 4531
 2640 / 4531
 2650 / 4531
 2660 / 4531
 2670 / 4531
 2680 / 4531
 2690 / 4531
 2700 / 4531
 2710 / 4531
 2720 / 4531
 2730 / 4531
 2740 / 4531
 2750 / 4531
 2760 / 4531
 2770 / 4531
 2780 / 4531
 2790 / 4531
 2800 / 4531
 2810 / 4531
 2820 / 4531
 2830 / 4531
 2840 / 4531
 2850 / 4531
 2860 / 4531
 2870 / 4531
 2880 / 4531
 2890 / 4531
 2900 / 4531
 2910 / 4531
 2920 / 4531
 2930 / 4531
 2940 / 4531
 2950 / 4531
 2960 / 4531
 2970 / 4531
 2980 / 4531
 2990 / 4531
 3000 / 4531
 3010 / 4531
 3020 / 4531
 3030 / 4531
 3040 / 4531
 3050 / 4531
 3060 / 4531
 3070 / 4531
 3080 / 4531
 3090 / 4531
 3100 / 4531
 3110 / 4531
 3120 / 4531
 3130 / 4531
 3140 / 4531
 3150 / 4531
 3160 / 4531
 3170 / 4531
 3180 / 4531
 3190 / 4531
 3200 / 4531
 3210 / 4531
 3220 / 4531
 3230 / 4531
 3240 / 4531
 3250 / 4531
 3260 / 4531
 3270 / 4531
 3280 / 4531
 3290 / 4531
 3300 / 4531
 3310 / 4531
 3320 / 4531
 3330 / 4531
 3340 / 4531
 3350 / 4531
 3360 / 4531
 3370 / 4531
 3380 / 4531
 3390 / 4531
 3400 / 4531
 3410 / 4531
 3420 / 4531
 3430 / 4531
 3440 / 4531
 3450 / 4531
 3460 / 4531
 3470 / 4531
 3480 / 4531
 3490 / 4531
 3500 / 4531
 3510 / 4531
 3520 / 4531
 3530 / 4531
 3540 / 4531
 3550 / 4531
 3560 / 4531
 3570 / 4531
 3580 / 4531
 3590 / 4531
 3600 / 4531
 3610 / 4531
 3620 / 4531
 3630 / 4531
 3640 / 4531
 3650 / 4531
 3660 / 4531
 3670 / 4531
 3680 / 4531
 3690 / 4531
 3700 / 4531
 3710 / 4531
 3720 / 4531
 3730 / 4531
 3740 / 4531
 3750 / 4531
 3760 / 4531
 3770 / 4531
 3780 / 4531
 3790 / 4531
 3800 / 4531
 3810 / 4531
 3820 / 4531
 3830 / 4531
 3840 / 4531
 3850 / 4531
 3860 / 4531
 3870 / 4531
 3880 / 4531
 3890 / 4531
 3900 / 4531
 3910 / 4531
 3920 / 4531
 3930 / 4531
 3940 / 4531
 3950 / 4531
 3960 / 4531
 3970 / 4531
 3980 / 4531
 3990 / 4531
 4000 / 4531
 4010 / 4531
 4020 / 4531
 4030 / 4531
 4040 / 4531
 4050 / 4531
 4060 / 4531
 4070 / 4531
 4080 / 4531
 4090 / 4531
 4100 / 4531
 4110 / 4531
 4120 / 4531
 4130 / 4531
 4140 / 4531
 4150 / 4531
 4160 / 4531
 4170 / 4531
 4180 / 4531
 4190 / 4531
 4200 / 4531
 4210 / 4531
 4220 / 4531
 4230 / 4531
 4240 / 4531
 4250 / 4531
 4260 / 4531
 4270 / 4531
 4280 / 4531
 4290 / 4531
 4300 / 4531
 4310 / 4531
 4320 / 4531
 4330 / 4531
 4340 / 4531
 4350 / 4531
 4360 / 4531
 4370 / 4531
 4380 / 4531
 4390 / 4531
 4400 / 4531
 4410 / 4531
 4420 / 4531
 4430 / 4531
 4440 / 4531
 4450 / 4531
 4460 / 4531
 4470 / 4531
 4480 / 4531
 4490 / 4531
 4500 / 4531
 4510 / 4531
 4520 / 4531
## done.
## 
## Progressive alignment 1/2... 
## 
STEP     1 / 4530  f
STEP     2 / 4530  f
STEP     3 / 4530  f
STEP     4 / 4530  f
STEP     5 / 4530  f
STEP     6 / 4530  f
STEP     7 / 4530  f
STEP     8 / 4530  f
STEP     9 / 4530  f
STEP    10 / 4530  f
STEP    11 / 4530  f
STEP    12 / 4530  f
STEP    13 / 4530  f
STEP    14 / 4530  f
STEP    15 / 4530  f
STEP    16 / 4530  f
STEP    17 / 4530  f
STEP    18 / 4530  f
STEP    19 / 4530  f
STEP    20 / 4530  f
STEP    21 / 4530  f
STEP    22 / 4530  f
STEP    23 / 4530  f
STEP    24 / 4530  f
STEP    25 / 4530  f
STEP    26 / 4530  f
STEP    27 / 4530  f
STEP    28 / 4530  f
STEP    29 / 4530  f
STEP    30 / 4530  f
STEP    31 / 4530  f
STEP    32 / 4530  f
STEP    33 / 4530  f
STEP    34 / 4530  f
STEP    35 / 4530  f
STEP    36 / 4530  f
STEP    37 / 4530  f
STEP    38 / 4530  f
STEP    39 / 4530  f
STEP    40 / 4530  f
STEP    41 / 4530  f
STEP    42 / 4530  f
STEP    43 / 4530  f
STEP    44 / 4530  f
STEP    45 / 4530  f
STEP    46 / 4530  f
STEP    47 / 4530  f
STEP    48 / 4530  f
STEP    49 / 4530  f
STEP    50 / 4530  f
STEP    51 / 4530  f
STEP    52 / 4530  f
STEP    53 / 4530  f
STEP    54 / 4530  f
STEP    55 / 4530  f
STEP    56 / 4530  f
STEP    57 / 4530  f
STEP    58 / 4530  f
STEP    59 / 4530  f
STEP    60 / 4530  f
STEP    61 / 4530  f
STEP    62 / 4530  f
STEP    63 / 4530  f
STEP    64 / 4530  f
STEP    65 / 4530  f
STEP    66 / 4530  f
STEP    67 / 4530  f
STEP    68 / 4530  f
STEP    69 / 4530  f
STEP    70 / 4530  f
STEP    71 / 4530  f
STEP    72 / 4530  f
STEP    73 / 4530  f
STEP    74 / 4530  f
STEP    75 / 4530  f
STEP    76 / 4530  f
STEP    77 / 4530  f
STEP    78 / 4530  f
STEP    79 / 4530  f
STEP    80 / 4530  f
STEP    81 / 4530  f
STEP    82 / 4530  f
STEP    83 / 4530  f
STEP    84 / 4530  f
STEP    85 / 4530  f
STEP    86 / 4530  f
STEP    87 / 4530  f
STEP    88 / 4530  f
STEP    89 / 4530  f
STEP    90 / 4530  f
STEP    91 / 4530  f
STEP    92 / 4530  f
STEP    93 / 4530  f
STEP    94 / 4530  f
STEP    95 / 4530  f
STEP    96 / 4530  f
STEP    97 / 4530  f
STEP    98 / 4530  f
STEP    99 / 4530  f
STEP   100 / 4530  f
STEP   101 / 4530  f
STEP   102 / 4530  f
STEP   103 / 4530  f
STEP   104 / 4530  f
STEP   105 / 4530  f
STEP   106 / 4530  f
STEP   107 / 4530  f
STEP   108 / 4530  f
STEP   109 / 4530  f
STEP   110 / 4530  f
STEP   111 / 4530  f
STEP   112 / 4530  f
STEP   113 / 4530  f
STEP   114 / 4530  f
STEP   115 / 4530  f
STEP   116 / 4530  f
STEP   117 / 4530  f
STEP   118 / 4530  f
STEP   119 / 4530  f
STEP   120 / 4530  f
STEP   121 / 4530  f
STEP   122 / 4530  f
STEP   123 / 4530  f
STEP   124 / 4530  f
STEP   125 / 4530  f
STEP   126 / 4530  f
STEP   127 / 4530  f
STEP   128 / 4530  f
STEP   129 / 4530  f
STEP   130 / 4530  f
STEP   131 / 4530  f
STEP   132 / 4530  f
STEP   133 / 4530  f
STEP   134 / 4530  f
STEP   135 / 4530  f
STEP   136 / 4530  f
STEP   137 / 4530  f
STEP   138 / 4530  f
STEP   139 / 4530  f
STEP   140 / 4530  f
STEP   141 / 4530  f
STEP   142 / 4530  f
STEP   143 / 4530  f
STEP   144 / 4530  f
STEP   145 / 4530  f
STEP   146 / 4530  f
STEP   147 / 4530  f
STEP   148 / 4530  f
STEP   149 / 4530  f
STEP   150 / 4530  f
STEP   151 / 4530  f
STEP   152 / 4530  f
STEP   153 / 4530  f
STEP   154 / 4530  f
STEP   155 / 4530  f
STEP   156 / 4530  f
STEP   157 / 4530  f
STEP   158 / 4530  f
STEP   159 / 4530  f
STEP   160 / 4530  f
STEP   161 / 4530  f
STEP   162 / 4530  f
STEP   163 / 4530  f
STEP   164 / 4530  f
STEP   165 / 4530  f
STEP   166 / 4530  f
STEP   167 / 4530  f
STEP   168 / 4530  f
STEP   169 / 4530  f
STEP   170 / 4530  f
STEP   171 / 4530  f
STEP   172 / 4530  f
STEP   173 / 4530  f
STEP   174 / 4530  f
STEP   175 / 4530  f
STEP   176 / 4530  f
STEP   177 / 4530  f
STEP   178 / 4530  f
STEP   179 / 4530  f
STEP   180 / 4530  f
STEP   181 / 4530  f
STEP   182 / 4530  f
STEP   183 / 4530  f
STEP   184 / 4530  f
STEP   185 / 4530  f
STEP   186 / 4530  f
STEP   187 / 4530  f
STEP   188 / 4530  f
STEP   189 / 4530  f
STEP   190 / 4530  f
STEP   191 / 4530  f
STEP   192 / 4530  f
STEP   193 / 4530  f
STEP   194 / 4530  f
STEP   195 / 4530  f
STEP   196 / 4530  f
STEP   197 / 4530  f
STEP   198 / 4530  f
STEP   199 / 4530  f
STEP   200 / 4530  f
STEP   201 / 4530  f
STEP   202 / 4530  f
STEP   203 / 4530  f
STEP   204 / 4530  f
STEP   205 / 4530  f
STEP   206 / 4530  f
STEP   207 / 4530  f
STEP   208 / 4530  f
STEP   209 / 4530  f
STEP   210 / 4530  f
STEP   211 / 4530  f
STEP   212 / 4530  f
STEP   213 / 4530  f
STEP   214 / 4530  f
STEP   215 / 4530  f
STEP   216 / 4530  f
STEP   217 / 4530  f
STEP   218 / 4530  f
STEP   219 / 4530  f
STEP   220 / 4530  f
STEP   221 / 4530  f
STEP   222 / 4530  f
STEP   223 / 4530  f
STEP   224 / 4530  f
STEP   225 / 4530  f
STEP   226 / 4530  f
STEP   227 / 4530  f
STEP   228 / 4530  f
STEP   229 / 4530  f
STEP   230 / 4530  f
STEP   231 / 4530  f
STEP   232 / 4530  f
STEP   233 / 4530  f
STEP   234 / 4530  f
STEP   235 / 4530  f
STEP   236 / 4530  f
STEP   237 / 4530  f
STEP   238 / 4530  f
STEP   239 / 4530  f
STEP   240 / 4530  f
STEP   241 / 4530  f
STEP   242 / 4530  f
STEP   243 / 4530  f
STEP   244 / 4530  f
STEP   245 / 4530  f
STEP   246 / 4530  f
STEP   247 / 4530  f
STEP   248 / 4530  f
STEP   249 / 4530  f
STEP   250 / 4530  f
STEP   251 / 4530  f
STEP   252 / 4530  f
STEP   253 / 4530  f
STEP   254 / 4530  f
STEP   255 / 4530  f
STEP   256 / 4530  f
STEP   257 / 4530  f
STEP   258 / 4530  f
STEP   259 / 4530  f
STEP   260 / 4530  f
STEP   261 / 4530  f
STEP   262 / 4530  f
STEP   263 / 4530  f
STEP   264 / 4530  f
STEP   265 / 4530  f
STEP   266 / 4530  f
STEP   267 / 4530  f
STEP   268 / 4530  f
STEP   269 / 4530  f
STEP   270 / 4530  f
STEP   271 / 4530  f
STEP   272 / 4530  f
STEP   273 / 4530  f
STEP   274 / 4530  f
STEP   275 / 4530  f
STEP   276 / 4530  f
STEP   277 / 4530  f
STEP   278 / 4530  f
STEP   279 / 4530  f
STEP   280 / 4530  f
STEP   281 / 4530  f
STEP   282 / 4530  f
STEP   283 / 4530  f
STEP   284 / 4530  f
STEP   285 / 4530  f
STEP   286 / 4530  f
STEP   287 / 4530  f
STEP   288 / 4530  f
STEP   289 / 4530  f
STEP   290 / 4530  f
STEP   291 / 4530  f
STEP   292 / 4530  f
STEP   293 / 4530  f
STEP   294 / 4530  f
STEP   295 / 4530  f
STEP   296 / 4530  f
STEP   297 / 4530  f
STEP   298 / 4530  f
STEP   299 / 4530  f
STEP   300 / 4530  f
STEP   301 / 4530  f
STEP   302 / 4530  f
STEP   303 / 4530  f
STEP   304 / 4530  f
STEP   305 / 4530  f
STEP   306 / 4530  f
STEP   307 / 4530  f
STEP   308 / 4530  f
STEP   309 / 4530  f
STEP   310 / 4530  f
STEP   311 / 4530  f
STEP   312 / 4530  f
STEP   313 / 4530  f
STEP   314 / 4530  f
STEP   315 / 4530  f
STEP   316 / 4530  f
STEP   317 / 4530  f
STEP   318 / 4530  f
STEP   319 / 4530  f
STEP   320 / 4530  f
STEP   321 / 4530  f
STEP   322 / 4530  f
STEP   323 / 4530  f
STEP   324 / 4530  f
STEP   325 / 4530  f
STEP   326 / 4530  f
STEP   327 / 4530  f
STEP   328 / 4530  f
STEP   329 / 4530  f
STEP   330 / 4530  f
STEP   331 / 4530  f
STEP   332 / 4530  f
STEP   333 / 4530  f
STEP   334 / 4530  f
STEP   335 / 4530  f
STEP   336 / 4530  f
STEP   337 / 4530  f
STEP   338 / 4530  f
STEP   339 / 4530  f
STEP   340 / 4530  f
STEP   341 / 4530  f
STEP   342 / 4530  f
STEP   343 / 4530  f
STEP   344 / 4530  f
STEP   345 / 4530  f
STEP   346 / 4530  f
STEP   347 / 4530  f
STEP   348 / 4530  f
STEP   349 / 4530  f
STEP   350 / 4530  f
STEP   351 / 4530  f
STEP   352 / 4530  f
STEP   353 / 4530  f
STEP   354 / 4530  f
STEP   355 / 4530  f
STEP   356 / 4530  f
STEP   357 / 4530  f
STEP   358 / 4530  f
STEP   359 / 4530  f
STEP   360 / 4530  f
STEP   361 / 4530  f
STEP   362 / 4530  f
STEP   363 / 4530  f
STEP   364 / 4530  f
STEP   365 / 4530  f
STEP   366 / 4530  f
STEP   367 / 4530  f
STEP   368 / 4530  f
STEP   369 / 4530  f
STEP   370 / 4530  f
STEP   371 / 4530  f
STEP   372 / 4530  f
STEP   373 / 4530  f
STEP   374 / 4530  f
STEP   375 / 4530  f
STEP   376 / 4530  f
STEP   377 / 4530  f
STEP   378 / 4530  f
STEP   379 / 4530  f
STEP   380 / 4530  f
STEP   381 / 4530  f
STEP   382 / 4530  f
STEP   383 / 4530  f
STEP   384 / 4530  f
STEP   385 / 4530  f
STEP   386 / 4530  f
STEP   387 / 4530  f
STEP   388 / 4530  f
STEP   389 / 4530  f
STEP   390 / 4530  f
STEP   391 / 4530  f
STEP   392 / 4530  f
STEP   393 / 4530  f
STEP   394 / 4530  f
STEP   395 / 4530  f
STEP   396 / 4530  f
STEP   397 / 4530  f
STEP   398 / 4530  f
STEP   399 / 4530  f
STEP   400 / 4530  f
STEP   401 / 4530  f
STEP   402 / 4530  f
STEP   403 / 4530  f
STEP   404 / 4530  f
STEP   405 / 4530  f
STEP   406 / 4530  f
STEP   407 / 4530  f
STEP   408 / 4530  f
STEP   409 / 4530  f
STEP   410 / 4530  f
STEP   411 / 4530  f
STEP   412 / 4530  f
STEP   413 / 4530  f
STEP   414 / 4530  f
STEP   415 / 4530  f
STEP   416 / 4530  f
STEP   417 / 4530  f
STEP   418 / 4530  f
STEP   419 / 4530  f
STEP   420 / 4530  f
STEP   421 / 4530  f
STEP   422 / 4530  f
STEP   423 / 4530  f
STEP   424 / 4530  f
STEP   425 / 4530  f
STEP   426 / 4530  f
STEP   427 / 4530  f
STEP   428 / 4530  f
STEP   429 / 4530  f
STEP   430 / 4530  f
STEP   431 / 4530  f
STEP   432 / 4530  f
STEP   433 / 4530  f
STEP   434 / 4530  f
STEP   435 / 4530  f
STEP   436 / 4530  f
STEP   437 / 4530  f
STEP   438 / 4530  f
STEP   439 / 4530  f
STEP   440 / 4530  f
STEP   441 / 4530  f
STEP   442 / 4530  f
STEP   443 / 4530  f
STEP   444 / 4530  f
STEP   445 / 4530  f
STEP   446 / 4530  f
STEP   447 / 4530  f
STEP   448 / 4530  f
STEP   449 / 4530  f
STEP   450 / 4530  f
STEP   451 / 4530  f
STEP   452 / 4530  f
STEP   453 / 4530  f
STEP   454 / 4530  f
STEP   455 / 4530  f
STEP   456 / 4530  f
STEP   457 / 4530  f
STEP   458 / 4530  f
STEP   459 / 4530  f
STEP   460 / 4530  f
STEP   461 / 4530  f
STEP   462 / 4530  f
STEP   463 / 4530  f
STEP   464 / 4530  f
STEP   465 / 4530  f
STEP   466 / 4530  f
STEP   467 / 4530  f
STEP   468 / 4530  f
STEP   469 / 4530  f
STEP   470 / 4530  f
STEP   471 / 4530  f
STEP   472 / 4530  f
STEP   473 / 4530  f
STEP   474 / 4530  f
STEP   475 / 4530  f
STEP   476 / 4530  f
STEP   477 / 4530  f
STEP   478 / 4530  f
STEP   479 / 4530  f
STEP   480 / 4530  f
STEP   481 / 4530  f
STEP   482 / 4530  f
STEP   483 / 4530  f
STEP   484 / 4530  f
STEP   485 / 4530  f
STEP   486 / 4530  f
STEP   487 / 4530  f
STEP   488 / 4530  f
STEP   489 / 4530  f
STEP   490 / 4530  f
STEP   491 / 4530  f
STEP   492 / 4530  f
STEP   493 / 4530  f
STEP   494 / 4530  f
STEP   495 / 4530  f
STEP   496 / 4530  f
STEP   497 / 4530  f
STEP   498 / 4530  f
STEP   499 / 4530  f
STEP   500 / 4530  f
STEP   501 / 4530  f
STEP   601 / 4530  f
STEP   701 / 4530  f
STEP   801 / 4530  f
STEP   901 / 4530  f
STEP  1001 / 4530  f
STEP  1101 / 4530  f
STEP  1201 / 4530  f
STEP  1301 / 4530  f
STEP  1401 / 4530  f
STEP  1501 / 4530  f
STEP  1601 / 4530  f
STEP  1701 / 4530  f
STEP  1801 / 4530  f
STEP  1901 / 4530  f
STEP  2001 / 4530  f
STEP  2101 / 4530  f
STEP  2201 / 4530  f
STEP  2301 / 4530  f
STEP  2401 / 4530  f
STEP  2501 / 4530  f
STEP  2601 / 4530  f
STEP  2701 / 4530  f
STEP  2801 / 4530  f
STEP  2901 / 4530  f
STEP  3001 / 4530  f
STEP  3101 / 4530  f
STEP  3201 / 4530  f
STEP  3301 / 4530  f
## Reallocating..done. *alloclen = 1935
## 
STEP  3401 / 4530  f
STEP  3501 / 4530  f
STEP  3601 / 4530  f
STEP  3701 / 4530  f
STEP  3801 / 4530  f
STEP  3901 / 4530  f
STEP  4001 / 4530  f
STEP  4101 / 4530  f
STEP  4201 / 4530  f
STEP  4301 / 4530  f
STEP  4401 / 4530  f
STEP  4501 / 4530  d h
## done.
## 
## Making a distance matrix from msa.. 
## 
    0 / 4531
  100 / 4531
  200 / 4531
  300 / 4531
  400 / 4531
  500 / 4531
  600 / 4531
  700 / 4531
  800 / 4531
  900 / 4531
 1000 / 4531
 1100 / 4531
 1200 / 4531
 1300 / 4531
 1400 / 4531
 1500 / 4531
 1600 / 4531
 1700 / 4531
 1800 / 4531
 1900 / 4531
 2000 / 4531
 2100 / 4531
 2200 / 4531
 2300 / 4531
 2400 / 4531
 2500 / 4531
 2600 / 4531
 2700 / 4531
 2800 / 4531
 2900 / 4531
 3000 / 4531
 3100 / 4531
 3200 / 4531
 3300 / 4531
 3400 / 4531
 3500 / 4531
 3600 / 4531
 3700 / 4531
 3800 / 4531
 3900 / 4531
 4000 / 4531
 4100 / 4531
 4200 / 4531
 4300 / 4531
 4400 / 4531
 4500 / 4531
## done.
## 
## Constructing a UPGMA tree (efffree=1) ... 
## 
    0 / 4531
   10 / 4531
   20 / 4531
   30 / 4531
   40 / 4531
   50 / 4531
   60 / 4531
   70 / 4531
   80 / 4531
   90 / 4531
  100 / 4531
  110 / 4531
  120 / 4531
  130 / 4531
  140 / 4531
  150 / 4531
  160 / 4531
  170 / 4531
  180 / 4531
  190 / 4531
  200 / 4531
  210 / 4531
  220 / 4531
  230 / 4531
  240 / 4531
  250 / 4531
  260 / 4531
  270 / 4531
  280 / 4531
  290 / 4531
  300 / 4531
  310 / 4531
  320 / 4531
  330 / 4531
  340 / 4531
  350 / 4531
  360 / 4531
  370 / 4531
  380 / 4531
  390 / 4531
  400 / 4531
  410 / 4531
  420 / 4531
  430 / 4531
  440 / 4531
  450 / 4531
  460 / 4531
  470 / 4531
  480 / 4531
  490 / 4531
  500 / 4531
  510 / 4531
  520 / 4531
  530 / 4531
  540 / 4531
  550 / 4531
  560 / 4531
  570 / 4531
  580 / 4531
  590 / 4531
  600 / 4531
  610 / 4531
  620 / 4531
  630 / 4531
  640 / 4531
  650 / 4531
  660 / 4531
  670 / 4531
  680 / 4531
  690 / 4531
  700 / 4531
  710 / 4531
  720 / 4531
  730 / 4531
  740 / 4531
  750 / 4531
  760 / 4531
  770 / 4531
  780 / 4531
  790 / 4531
  800 / 4531
  810 / 4531
  820 / 4531
  830 / 4531
  840 / 4531
  850 / 4531
  860 / 4531
  870 / 4531
  880 / 4531
  890 / 4531
  900 / 4531
  910 / 4531
  920 / 4531
  930 / 4531
  940 / 4531
  950 / 4531
  960 / 4531
  970 / 4531
  980 / 4531
  990 / 4531
 1000 / 4531
 1010 / 4531
 1020 / 4531
 1030 / 4531
 1040 / 4531
 1050 / 4531
 1060 / 4531
 1070 / 4531
 1080 / 4531
 1090 / 4531
 1100 / 4531
 1110 / 4531
 1120 / 4531
 1130 / 4531
 1140 / 4531
 1150 / 4531
 1160 / 4531
 1170 / 4531
 1180 / 4531
 1190 / 4531
 1200 / 4531
 1210 / 4531
 1220 / 4531
 1230 / 4531
 1240 / 4531
 1250 / 4531
 1260 / 4531
 1270 / 4531
 1280 / 4531
 1290 / 4531
 1300 / 4531
 1310 / 4531
 1320 / 4531
 1330 / 4531
 1340 / 4531
 1350 / 4531
 1360 / 4531
 1370 / 4531
 1380 / 4531
 1390 / 4531
 1400 / 4531
 1410 / 4531
 1420 / 4531
 1430 / 4531
 1440 / 4531
 1450 / 4531
 1460 / 4531
 1470 / 4531
 1480 / 4531
 1490 / 4531
 1500 / 4531
 1510 / 4531
 1520 / 4531
 1530 / 4531
 1540 / 4531
 1550 / 4531
 1560 / 4531
 1570 / 4531
 1580 / 4531
 1590 / 4531
 1600 / 4531
 1610 / 4531
 1620 / 4531
 1630 / 4531
 1640 / 4531
 1650 / 4531
 1660 / 4531
 1670 / 4531
 1680 / 4531
 1690 / 4531
 1700 / 4531
 1710 / 4531
 1720 / 4531
 1730 / 4531
 1740 / 4531
 1750 / 4531
 1760 / 4531
 1770 / 4531
 1780 / 4531
 1790 / 4531
 1800 / 4531
 1810 / 4531
 1820 / 4531
 1830 / 4531
 1840 / 4531
 1850 / 4531
 1860 / 4531
 1870 / 4531
 1880 / 4531
 1890 / 4531
 1900 / 4531
 1910 / 4531
 1920 / 4531
 1930 / 4531
 1940 / 4531
 1950 / 4531
 1960 / 4531
 1970 / 4531
 1980 / 4531
 1990 / 4531
 2000 / 4531
 2010 / 4531
 2020 / 4531
 2030 / 4531
 2040 / 4531
 2050 / 4531
 2060 / 4531
 2070 / 4531
 2080 / 4531
 2090 / 4531
 2100 / 4531
 2110 / 4531
 2120 / 4531
 2130 / 4531
 2140 / 4531
 2150 / 4531
 2160 / 4531
 2170 / 4531
 2180 / 4531
 2190 / 4531
 2200 / 4531
 2210 / 4531
 2220 / 4531
 2230 / 4531
 2240 / 4531
 2250 / 4531
 2260 / 4531
 2270 / 4531
 2280 / 4531
 2290 / 4531
 2300 / 4531
 2310 / 4531
 2320 / 4531
 2330 / 4531
 2340 / 4531
 2350 / 4531
 2360 / 4531
 2370 / 4531
 2380 / 4531
 2390 / 4531
 2400 / 4531
 2410 / 4531
 2420 / 4531
 2430 / 4531
 2440 / 4531
 2450 / 4531
 2460 / 4531
 2470 / 4531
 2480 / 4531
 2490 / 4531
 2500 / 4531
 2510 / 4531
 2520 / 4531
 2530 / 4531
 2540 / 4531
 2550 / 4531
 2560 / 4531
 2570 / 4531
 2580 / 4531
 2590 / 4531
 2600 / 4531
 2610 / 4531
 2620 / 4531
 2630 / 4531
 2640 / 4531
 2650 / 4531
 2660 / 4531
 2670 / 4531
 2680 / 4531
 2690 / 4531
 2700 / 4531
 2710 / 4531
 2720 / 4531
 2730 / 4531
 2740 / 4531
 2750 / 4531
 2760 / 4531
 2770 / 4531
 2780 / 4531
 2790 / 4531
 2800 / 4531
 2810 / 4531
 2820 / 4531
 2830 / 4531
 2840 / 4531
 2850 / 4531
 2860 / 4531
 2870 / 4531
 2880 / 4531
 2890 / 4531
 2900 / 4531
 2910 / 4531
 2920 / 4531
 2930 / 4531
 2940 / 4531
 2950 / 4531
 2960 / 4531
 2970 / 4531
 2980 / 4531
 2990 / 4531
 3000 / 4531
 3010 / 4531
 3020 / 4531
 3030 / 4531
 3040 / 4531
 3050 / 4531
 3060 / 4531
 3070 / 4531
 3080 / 4531
 3090 / 4531
 3100 / 4531
 3110 / 4531
 3120 / 4531
 3130 / 4531
 3140 / 4531
 3150 / 4531
 3160 / 4531
 3170 / 4531
 3180 / 4531
 3190 / 4531
 3200 / 4531
 3210 / 4531
 3220 / 4531
 3230 / 4531
 3240 / 4531
 3250 / 4531
 3260 / 4531
 3270 / 4531
 3280 / 4531
 3290 / 4531
 3300 / 4531
 3310 / 4531
 3320 / 4531
 3330 / 4531
 3340 / 4531
 3350 / 4531
 3360 / 4531
 3370 / 4531
 3380 / 4531
 3390 / 4531
 3400 / 4531
 3410 / 4531
 3420 / 4531
 3430 / 4531
 3440 / 4531
 3450 / 4531
 3460 / 4531
 3470 / 4531
 3480 / 4531
 3490 / 4531
 3500 / 4531
 3510 / 4531
 3520 / 4531
 3530 / 4531
 3540 / 4531
 3550 / 4531
 3560 / 4531
 3570 / 4531
 3580 / 4531
 3590 / 4531
 3600 / 4531
 3610 / 4531
 3620 / 4531
 3630 / 4531
 3640 / 4531
 3650 / 4531
 3660 / 4531
 3670 / 4531
 3680 / 4531
 3690 / 4531
 3700 / 4531
 3710 / 4531
 3720 / 4531
 3730 / 4531
 3740 / 4531
 3750 / 4531
 3760 / 4531
 3770 / 4531
 3780 / 4531
 3790 / 4531
 3800 / 4531
 3810 / 4531
 3820 / 4531
 3830 / 4531
 3840 / 4531
 3850 / 4531
 3860 / 4531
 3870 / 4531
 3880 / 4531
 3890 / 4531
 3900 / 4531
 3910 / 4531
 3920 / 4531
 3930 / 4531
 3940 / 4531
 3950 / 4531
 3960 / 4531
 3970 / 4531
 3980 / 4531
 3990 / 4531
 4000 / 4531
 4010 / 4531
 4020 / 4531
 4030 / 4531
 4040 / 4531
 4050 / 4531
 4060 / 4531
 4070 / 4531
 4080 / 4531
 4090 / 4531
 4100 / 4531
 4110 / 4531
 4120 / 4531
 4130 / 4531
 4140 / 4531
 4150 / 4531
 4160 / 4531
 4170 / 4531
 4180 / 4531
 4190 / 4531
 4200 / 4531
 4210 / 4531
 4220 / 4531
 4230 / 4531
 4240 / 4531
 4250 / 4531
 4260 / 4531
 4270 / 4531
 4280 / 4531
 4290 / 4531
 4300 / 4531
 4310 / 4531
 4320 / 4531
 4330 / 4531
 4340 / 4531
 4350 / 4531
 4360 / 4531
 4370 / 4531
 4380 / 4531
 4390 / 4531
 4400 / 4531
 4410 / 4531
 4420 / 4531
 4430 / 4531
 4440 / 4531
 4450 / 4531
 4460 / 4531
 4470 / 4531
 4480 / 4531
 4490 / 4531
 4500 / 4531
 4510 / 4531
 4520 / 4531
## done.
## 
## Progressive alignment 2/2... 
## 
STEP     1 / 4530  f
STEP     2 / 4530  f
STEP     3 / 4530  f
STEP     4 / 4530  f
STEP     5 / 4530  f
STEP     6 / 4530  f
STEP     7 / 4530  f
STEP     8 / 4530  f
STEP     9 / 4530  f
STEP    10 / 4530  f
STEP    11 / 4530  f
STEP    12 / 4530  f
STEP    13 / 4530  f
STEP    14 / 4530  f
STEP    15 / 4530  f
STEP    16 / 4530  f
STEP    17 / 4530  f
STEP    18 / 4530  f
STEP    19 / 4530  f
STEP    20 / 4530  f
STEP    21 / 4530  f
STEP    22 / 4530  f
STEP    23 / 4530  f
STEP    24 / 4530  f
STEP    25 / 4530  f
STEP    26 / 4530  f
STEP    27 / 4530  f
STEP    28 / 4530  f
STEP    29 / 4530  f
STEP    30 / 4530  f
STEP    31 / 4530  f
STEP    32 / 4530  f
STEP    33 / 4530  f
STEP    34 / 4530  f
STEP    35 / 4530  f
STEP    36 / 4530  f
STEP    37 / 4530  f
STEP    38 / 4530  f
STEP    39 / 4530  f
STEP    40 / 4530  f
STEP    41 / 4530  f
STEP    42 / 4530  f
STEP    43 / 4530  f
STEP    44 / 4530  f
STEP    45 / 4530  f
STEP    46 / 4530  f
STEP    47 / 4530  f
STEP    48 / 4530  f
STEP    49 / 4530  f
STEP    50 / 4530  f
STEP    51 / 4530  f
STEP    52 / 4530  f
STEP    53 / 4530  f
STEP    54 / 4530  f
STEP    55 / 4530  f
STEP    56 / 4530  f
STEP    57 / 4530  f
STEP    58 / 4530  f
STEP    59 / 4530  f
STEP    60 / 4530  f
STEP    61 / 4530  f
STEP    62 / 4530  f
STEP    63 / 4530  f
STEP    64 / 4530  f
STEP    65 / 4530  f
STEP    66 / 4530  f
STEP    67 / 4530  f
STEP    68 / 4530  f
STEP    69 / 4530  f
STEP    70 / 4530  f
STEP    71 / 4530  f
STEP    72 / 4530  f
STEP    73 / 4530  f
STEP    74 / 4530  f
STEP    75 / 4530  f
STEP    76 / 4530  f
STEP    77 / 4530  f
STEP    78 / 4530  f
STEP    79 / 4530  f
STEP    80 / 4530  f
STEP    81 / 4530  f
STEP    82 / 4530  f
STEP    83 / 4530  f
STEP    84 / 4530  f
STEP    85 / 4530  f
STEP    86 / 4530  f
STEP    87 / 4530  f
STEP    88 / 4530  f
STEP    89 / 4530  f
STEP    90 / 4530  f
STEP    91 / 4530  f
STEP    92 / 4530  f
STEP    93 / 4530  f
STEP    94 / 4530  f
STEP    95 / 4530  f
STEP    96 / 4530  f
STEP    97 / 4530  f
STEP    98 / 4530  f
STEP    99 / 4530  f
STEP   100 / 4530  f
STEP   101 / 4530  f
STEP   102 / 4530  f
STEP   103 / 4530  f
STEP   104 / 4530  f
STEP   105 / 4530  f
STEP   106 / 4530  f
STEP   107 / 4530  f
STEP   108 / 4530  f
STEP   109 / 4530  f
STEP   110 / 4530  f
STEP   111 / 4530  f
STEP   112 / 4530  f
STEP   113 / 4530  f
STEP   114 / 4530  f
STEP   115 / 4530  f
STEP   116 / 4530  f
STEP   117 / 4530  f
STEP   118 / 4530  f
STEP   119 / 4530  f
STEP   120 / 4530  f
STEP   121 / 4530  f
STEP   122 / 4530  f
STEP   123 / 4530  f
STEP   124 / 4530  f
STEP   125 / 4530  f
STEP   126 / 4530  f
STEP   127 / 4530  f
STEP   128 / 4530  f
STEP   129 / 4530  f
STEP   130 / 4530  f
STEP   131 / 4530  f
STEP   132 / 4530  f
STEP   133 / 4530  f
STEP   134 / 4530  f
STEP   135 / 4530  f
STEP   136 / 4530  f
STEP   137 / 4530  f
STEP   138 / 4530  f
STEP   139 / 4530  f
STEP   140 / 4530  f
STEP   141 / 4530  f
STEP   142 / 4530  f
STEP   143 / 4530  f
STEP   144 / 4530  f
STEP   145 / 4530  f
STEP   146 / 4530  f
STEP   147 / 4530  f
STEP   148 / 4530  f
STEP   149 / 4530  f
STEP   150 / 4530  f
STEP   151 / 4530  f
STEP   152 / 4530  f
STEP   153 / 4530  f
STEP   154 / 4530  f
STEP   155 / 4530  f
STEP   156 / 4530  f
STEP   157 / 4530  f
STEP   158 / 4530  f
STEP   159 / 4530  f
STEP   160 / 4530  f
STEP   161 / 4530  f
STEP   162 / 4530  f
STEP   163 / 4530  f
STEP   164 / 4530  f
STEP   165 / 4530  f
STEP   166 / 4530  f
STEP   167 / 4530  f
STEP   168 / 4530  f
STEP   169 / 4530  f
STEP   170 / 4530  f
STEP   171 / 4530  f
STEP   172 / 4530  f
STEP   173 / 4530  f
STEP   174 / 4530  f
STEP   175 / 4530  f
STEP   176 / 4530  f
STEP   177 / 4530  f
STEP   178 / 4530  f
STEP   179 / 4530  f
STEP   180 / 4530  f
STEP   181 / 4530  f
STEP   182 / 4530  f
STEP   183 / 4530  f
STEP   184 / 4530  f
STEP   185 / 4530  f
STEP   186 / 4530  f
STEP   187 / 4530  f
STEP   188 / 4530  f
STEP   189 / 4530  f
STEP   190 / 4530  f
STEP   191 / 4530  f
STEP   192 / 4530  f
STEP   193 / 4530  f
STEP   194 / 4530  f
STEP   195 / 4530  f
STEP   196 / 4530  f
STEP   197 / 4530  f
STEP   198 / 4530  f
STEP   199 / 4530  f
STEP   200 / 4530  f
STEP   201 / 4530  f
STEP   202 / 4530  f
STEP   203 / 4530  f
STEP   204 / 4530  f
STEP   205 / 4530  f
STEP   206 / 4530  f
STEP   207 / 4530  f
STEP   208 / 4530  f
STEP   209 / 4530  f
STEP   210 / 4530  f
STEP   211 / 4530  f
STEP   212 / 4530  f
STEP   213 / 4530  f
STEP   214 / 4530  f
STEP   215 / 4530  f
STEP   216 / 4530  f
STEP   217 / 4530  f
STEP   218 / 4530  f
STEP   219 / 4530  f
STEP   220 / 4530  f
STEP   221 / 4530  f
STEP   222 / 4530  f
STEP   223 / 4530  f
STEP   224 / 4530  f
STEP   225 / 4530  f
STEP   226 / 4530  f
STEP   227 / 4530  f
STEP   228 / 4530  f
STEP   229 / 4530  f
STEP   230 / 4530  f
STEP   231 / 4530  f
STEP   232 / 4530  f
STEP   233 / 4530  f
STEP   234 / 4530  f
STEP   235 / 4530  f
STEP   236 / 4530  f
STEP   237 / 4530  f
STEP   238 / 4530  f
STEP   239 / 4530  f
STEP   240 / 4530  f
STEP   241 / 4530  f
STEP   242 / 4530  f
STEP   243 / 4530  f
STEP   244 / 4530  f
STEP   245 / 4530  f
STEP   246 / 4530  f
STEP   247 / 4530  f
STEP   248 / 4530  f
STEP   249 / 4530  f
STEP   250 / 4530  f
STEP   251 / 4530  f
STEP   252 / 4530  f
STEP   253 / 4530  f
STEP   254 / 4530  f
STEP   255 / 4530  f
STEP   256 / 4530  f
STEP   257 / 4530  f
STEP   258 / 4530  f
STEP   259 / 4530  f
STEP   260 / 4530  f
STEP   261 / 4530  f
STEP   262 / 4530  f
STEP   263 / 4530  f
STEP   264 / 4530  f
STEP   265 / 4530  f
STEP   266 / 4530  f
STEP   267 / 4530  f
STEP   268 / 4530  f
STEP   269 / 4530  f
STEP   270 / 4530  f
STEP   271 / 4530  f
STEP   272 / 4530  f
STEP   273 / 4530  f
STEP   274 / 4530  f
STEP   275 / 4530  f
STEP   276 / 4530  f
STEP   277 / 4530  f
STEP   278 / 4530  f
STEP   279 / 4530  f
STEP   280 / 4530  f
STEP   281 / 4530  f
STEP   282 / 4530  f
STEP   283 / 4530  f
STEP   284 / 4530  f
STEP   285 / 4530  f
STEP   286 / 4530  f
STEP   287 / 4530  f
STEP   288 / 4530  f
STEP   289 / 4530  f
STEP   290 / 4530  f
STEP   291 / 4530  f
STEP   292 / 4530  f
STEP   293 / 4530  f
STEP   294 / 4530  f
STEP   295 / 4530  f
STEP   296 / 4530  f
STEP   297 / 4530  f
STEP   298 / 4530  f
STEP   299 / 4530  f
STEP   300 / 4530  f
STEP   301 / 4530  f
STEP   302 / 4530  f
STEP   303 / 4530  f
STEP   304 / 4530  f
STEP   305 / 4530  f
STEP   306 / 4530  f
STEP   307 / 4530  f
STEP   308 / 4530  f
STEP   309 / 4530  f
STEP   310 / 4530  f
STEP   311 / 4530  f
STEP   312 / 4530  f
STEP   313 / 4530  f
STEP   314 / 4530  f
STEP   315 / 4530  f
STEP   316 / 4530  f
STEP   317 / 4530  f
STEP   318 / 4530  f
STEP   319 / 4530  f
STEP   320 / 4530  f
STEP   321 / 4530  f
STEP   322 / 4530  f
STEP   323 / 4530  f
STEP   324 / 4530  f
STEP   325 / 4530  f
STEP   326 / 4530  f
STEP   327 / 4530  f
STEP   328 / 4530  f
STEP   329 / 4530  f
STEP   330 / 4530  f
STEP   331 / 4530  f
STEP   332 / 4530  f
STEP   333 / 4530  f
STEP   334 / 4530  f
STEP   335 / 4530  f
STEP   336 / 4530  f
STEP   337 / 4530  f
STEP   338 / 4530  f
STEP   339 / 4530  f
STEP   340 / 4530  f
STEP   341 / 4530  f
STEP   342 / 4530  f
STEP   343 / 4530  f
STEP   344 / 4530  f
STEP   345 / 4530  f
STEP   346 / 4530  f
STEP   347 / 4530  f
STEP   348 / 4530  f
STEP   349 / 4530  f
STEP   350 / 4530  f
STEP   351 / 4530  f
STEP   352 / 4530  f
STEP   353 / 4530  f
STEP   354 / 4530  f
STEP   355 / 4530  f
STEP   356 / 4530  f
STEP   357 / 4530  f
STEP   358 / 4530  f
STEP   359 / 4530  f
STEP   360 / 4530  f
STEP   361 / 4530  f
STEP   362 / 4530  f
STEP   363 / 4530  f
STEP   364 / 4530  f
STEP   365 / 4530  f
STEP   366 / 4530  f
STEP   367 / 4530  f
STEP   368 / 4530  f
STEP   369 / 4530  f
STEP   370 / 4530  f
STEP   371 / 4530  f
STEP   372 / 4530  f
STEP   373 / 4530  f
STEP   374 / 4530  f
STEP   375 / 4530  f
STEP   376 / 4530  f
STEP   377 / 4530  f
STEP   378 / 4530  f
STEP   379 / 4530  f
STEP   380 / 4530  f
STEP   381 / 4530  f
STEP   382 / 4530  f
STEP   383 / 4530  f
STEP   384 / 4530  f
STEP   385 / 4530  f
STEP   386 / 4530  f
STEP   387 / 4530  f
STEP   388 / 4530  f
STEP   389 / 4530  f
STEP   390 / 4530  f
STEP   391 / 4530  f
STEP   392 / 4530  f
STEP   393 / 4530  f
STEP   394 / 4530  f
STEP   395 / 4530  f
STEP   396 / 4530  f
STEP   397 / 4530  f
STEP   398 / 4530  f
STEP   399 / 4530  f
STEP   400 / 4530  f
STEP   401 / 4530  f
STEP   402 / 4530  f
STEP   403 / 4530  f
STEP   404 / 4530  f
STEP   405 / 4530  f
STEP   406 / 4530  f
STEP   407 / 4530  f
STEP   408 / 4530  f
STEP   409 / 4530  f
STEP   410 / 4530  f
STEP   411 / 4530  f
STEP   412 / 4530  f
STEP   413 / 4530  f
STEP   414 / 4530  f
STEP   415 / 4530  f
STEP   416 / 4530  f
STEP   417 / 4530  f
STEP   418 / 4530  f
STEP   419 / 4530  f
STEP   420 / 4530  f
STEP   421 / 4530  f
STEP   422 / 4530  f
STEP   423 / 4530  f
STEP   424 / 4530  f
STEP   425 / 4530  f
STEP   426 / 4530  f
STEP   427 / 4530  f
STEP   428 / 4530  f
STEP   429 / 4530  f
STEP   430 / 4530  f
STEP   431 / 4530  f
STEP   432 / 4530  f
STEP   433 / 4530  f
STEP   434 / 4530  f
STEP   435 / 4530  f
STEP   436 / 4530  f
STEP   437 / 4530  f
STEP   438 / 4530  f
STEP   439 / 4530  f
STEP   440 / 4530  f
STEP   441 / 4530  f
STEP   442 / 4530  f
STEP   443 / 4530  f
STEP   444 / 4530  f
STEP   445 / 4530  f
STEP   446 / 4530  f
STEP   447 / 4530  f
STEP   448 / 4530  f
STEP   449 / 4530  f
STEP   450 / 4530  f
STEP   451 / 4530  f
STEP   452 / 4530  f
STEP   453 / 4530  f
STEP   454 / 4530  f
STEP   455 / 4530  f
STEP   456 / 4530  f
STEP   457 / 4530  f
STEP   458 / 4530  f
STEP   459 / 4530  f
STEP   460 / 4530  f
STEP   461 / 4530  f
STEP   462 / 4530  f
STEP   463 / 4530  f
STEP   464 / 4530  f
STEP   465 / 4530  f
STEP   466 / 4530  f
STEP   467 / 4530  f
STEP   468 / 4530  f
STEP   469 / 4530  f
STEP   470 / 4530  f
STEP   471 / 4530  f
STEP   472 / 4530  f
STEP   473 / 4530  f
STEP   474 / 4530  f
STEP   475 / 4530  f
STEP   476 / 4530  f
STEP   477 / 4530  f
STEP   478 / 4530  f
STEP   479 / 4530  f
STEP   480 / 4530  f
STEP   481 / 4530  f
STEP   482 / 4530  f
STEP   483 / 4530  f
STEP   484 / 4530  f
STEP   485 / 4530  f
STEP   486 / 4530  f
STEP   487 / 4530  f
STEP   488 / 4530  f
STEP   489 / 4530  f
STEP   490 / 4530  f
STEP   491 / 4530  f
STEP   492 / 4530  f
STEP   493 / 4530  f
STEP   494 / 4530  f
STEP   495 / 4530  f
STEP   496 / 4530  f
STEP   497 / 4530  f
STEP   498 / 4530  f
STEP   499 / 4530  f
STEP   500 / 4530  f
STEP   501 / 4530  f
STEP   601 / 4530  f
STEP   701 / 4530  f
STEP   801 / 4530  f
STEP   901 / 4530  f
STEP  1001 / 4530  f
STEP  1101 / 4530  f
STEP  1201 / 4530  f
STEP  1301 / 4530  f
STEP  1401 / 4530  f
STEP  1501 / 4530  f
STEP  1601 / 4530  f
STEP  1701 / 4530  f
STEP  1801 / 4530  f
STEP  1901 / 4530  f
STEP  2001 / 4530  f
STEP  2101 / 4530  f
STEP  2201 / 4530  f
STEP  2301 / 4530  f
STEP  2401 / 4530  f
STEP  2501 / 4530  f
STEP  2601 / 4530  f
STEP  2701 / 4530  f
STEP  2801 / 4530  f
STEP  2901 / 4530  f
STEP  3001 / 4530  f
STEP  3101 / 4530  f
STEP  3201 / 4530  f
STEP  3301 / 4530  f
STEP  3401 / 4530  f
STEP  3501 / 4530  f
STEP  3601 / 4530  f
STEP  3701 / 4530  f
STEP  3801 / 4530  f
STEP  3901 / 4530  f
## Reallocating..done. *alloclen = 1934
## 
STEP  4001 / 4530  f
STEP  4101 / 4530  f
STEP  4201 / 4530  f
STEP  4301 / 4530  f
STEP  4401 / 4530  f
STEP  4501 / 4530  f
## done.
## 
## disttbfast (nuc) Version 7.522
## alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
## 0 thread(s)
## 
## 
## Strategy:
##  FFT-NS-2 (Fast but rough)
##  Progressive method (guide trees were built 2 times.)
## 
## If unsure which option to use, try 'mafft --auto input > output'.
## For more information, see 'mafft --help', 'mafft --man' and the mafft page.
## 
## The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
## It tends to insert more gaps into gap-rich regions than previous versions.
## To disable this change, add the --leavegappyregion option.
## 
## /local/workdir/sna49/moon_milk/moonmilk
```



# FastTree2 

```bash
# Provide export path to fasttree 
export PATH=/programs/FastTree-2.1.11:$PATH

# cd into the alignment file folder 
cd data/03_Phylogenetic_Tree/
pwd

# Run Fasttree to generate phylogenetic tree 
# parameters: 
    # -nt = indicates it's a nucleotide alignment
    # -gtr = generalized time reversible substitution model 
    # -fastest speed up the model, reduce memory usage (recommended for datasets that have >50,000)
    # -log = output a log file 
    # input alignment file 
    # specify the output tree file 
FastTree -nt -gtr -fastest -log FastTree.log MAFFT_aligned_ASVs.fasta > ASVs_unrooted.tree

# Change back to the main directory 
cd ../../
pwd 
echo "The working directory is $PWD"
```

```
## /local/workdir/sna49/moon_milk/moonmilk/data/03_Phylogenetic_Tree
## FastTree Version 2.1.11 SSE3
## Alignment: MAFFT_aligned_ASVs.fasta
## Nucleotide distances: Jukes-Cantor Joins: balanced Support: SH-like 1000
## Search: Fastest+2nd +NNI +SPR (2 rounds range 10) +ML-NNI opt-each=1
## TopHits: 1.00*sqrtN close=default refresh=0.50
## ML Model: Generalized Time-Reversible, CAT approximation with 20 rate categories
##       0.10 seconds: Top hits for    370 of   4531 seqs (at seed    100)
##       0.22 seconds: Top hits for   1390 of   4531 seqs (at seed    500)
##       0.33 seconds: Top hits for   2161 of   4531 seqs (at seed   1000)
##       0.44 seconds: Top hits for   2628 of   4531 seqs (at seed   1600)
##       0.54 seconds: Top hits for   2847 of   4531 seqs (at seed   2100)
##       0.64 seconds: Top hits for   3517 of   4531 seqs (at seed   2400)
##       0.74 seconds: Top hits for   3997 of   4531 seqs (at seed   3100)
##       0.86 seconds: Top hits for   4412 of   4531 seqs (at seed   3700)
##       0.97 seconds: Joined    100 of   4528
##       1.11 seconds: Joined    300 of   4528
##       1.30 seconds: Joined    500 of   4528
##       1.48 seconds: Joined    600 of   4528
##       1.63 seconds: Joined    700 of   4528
##       1.74 seconds: Joined    800 of   4528
##       1.89 seconds: Joined   1000 of   4528
##       2.00 seconds: Joined   1100 of   4528
##       2.12 seconds: Joined   1200 of   4528
##       2.28 seconds: Joined   1400 of   4528
##       2.42 seconds: Joined   1500 of   4528
##       2.55 seconds: Joined   1700 of   4528
##       2.78 seconds: Joined   1900 of   4528
##       2.88 seconds: Joined   2000 of   4528
##       3.13 seconds: Joined   2200 of   4528
##       3.32 seconds: Joined   2400 of   4528
##       3.54 seconds: Joined   2600 of   4528
##       3.66 seconds: Joined   2700 of   4528
##       3.86 seconds: Joined   2800 of   4528
##       4.06 seconds: Joined   3000 of   4528
##       4.17 seconds: Joined   3100 of   4528
##       4.31 seconds: Joined   3200 of   4528
##       4.43 seconds: Joined   3300 of   4528
##       4.59 seconds: Joined   3500 of   4528
##       4.78 seconds: Joined   3700 of   4528
##       4.95 seconds: Joined   3800 of   4528
##       5.07 seconds: Joined   3900 of   4528
##       5.18 seconds: Joined   4000 of   4528
##       5.30 seconds: Joined   4100 of   4528
##       5.48 seconds: Joined   4200 of   4528
##       5.61 seconds: Joined   4400 of   4528
## Initial topology in 5.71 seconds
## Refining topology: 49 rounds ME-NNIs, 2 rounds ME-SPRs, 24 rounds ML-NNIs
##       5.72 seconds: ME NNI round 1 of 49, 301 of 4529 splits, 44 changes (max delta 0.031)
##       5.82 seconds: ME NNI round 1 of 49, 2701 of 4529 splits, 541 changes (max delta 0.071)
##       5.92 seconds: ME NNI round 2 of 49, 501 of 4529 splits, 59 changes (max delta 0.020)
##       6.02 seconds: ME NNI round 2 of 49, 2901 of 4529 splits, 294 changes (max delta 0.036)
##       6.12 seconds: ME NNI round 3 of 49, 501 of 4529 splits, 34 changes (max delta 0.012)
##       6.23 seconds: ME NNI round 3 of 49, 2801 of 4529 splits, 202 changes (max delta 0.018)
##       6.33 seconds: ME NNI round 4 of 49, 1601 of 4529 splits, 97 changes (max delta 0.021)
##       6.43 seconds: ME NNI round 6 of 49, 101 of 4529 splits, 7 changes (max delta 0.001)
##       6.61 seconds: SPR round   1 of   2, 101 of 9060 nodes
##       6.76 seconds: SPR round   1 of   2, 301 of 9060 nodes
##       6.93 seconds: SPR round   1 of   2, 501 of 9060 nodes
##       7.06 seconds: SPR round   1 of   2, 701 of 9060 nodes
##       7.19 seconds: SPR round   1 of   2, 901 of 9060 nodes
##       7.32 seconds: SPR round   1 of   2, 1101 of 9060 nodes
##       7.47 seconds: SPR round   1 of   2, 1301 of 9060 nodes
##       7.61 seconds: SPR round   1 of   2, 1501 of 9060 nodes
##       7.77 seconds: SPR round   1 of   2, 1701 of 9060 nodes
##       7.93 seconds: SPR round   1 of   2, 1901 of 9060 nodes
##       8.08 seconds: SPR round   1 of   2, 2101 of 9060 nodes
##       8.25 seconds: SPR round   1 of   2, 2301 of 9060 nodes
##       8.42 seconds: SPR round   1 of   2, 2501 of 9060 nodes
##       8.55 seconds: SPR round   1 of   2, 2701 of 9060 nodes
##       8.69 seconds: SPR round   1 of   2, 2901 of 9060 nodes
##       8.87 seconds: SPR round   1 of   2, 3101 of 9060 nodes
##       8.99 seconds: SPR round   1 of   2, 3301 of 9060 nodes
##       9.13 seconds: SPR round   1 of   2, 3501 of 9060 nodes
##       9.30 seconds: SPR round   1 of   2, 3701 of 9060 nodes
##       9.46 seconds: SPR round   1 of   2, 3901 of 9060 nodes
##       9.61 seconds: SPR round   1 of   2, 4101 of 9060 nodes
##       9.78 seconds: SPR round   1 of   2, 4301 of 9060 nodes
##       9.95 seconds: SPR round   1 of   2, 4501 of 9060 nodes
##      10.08 seconds: SPR round   1 of   2, 4701 of 9060 nodes
##      10.24 seconds: SPR round   1 of   2, 4901 of 9060 nodes
##      10.36 seconds: SPR round   1 of   2, 5101 of 9060 nodes
##      10.52 seconds: SPR round   1 of   2, 5301 of 9060 nodes
##      10.70 seconds: SPR round   1 of   2, 5501 of 9060 nodes
##      10.87 seconds: SPR round   1 of   2, 5701 of 9060 nodes
##      10.99 seconds: SPR round   1 of   2, 5901 of 9060 nodes
##      11.12 seconds: SPR round   1 of   2, 6101 of 9060 nodes
##      11.27 seconds: SPR round   1 of   2, 6301 of 9060 nodes
##      11.41 seconds: SPR round   1 of   2, 6501 of 9060 nodes
##      11.55 seconds: SPR round   1 of   2, 6701 of 9060 nodes
##      11.70 seconds: SPR round   1 of   2, 6901 of 9060 nodes
##      11.86 seconds: SPR round   1 of   2, 7101 of 9060 nodes
##      11.99 seconds: SPR round   1 of   2, 7301 of 9060 nodes
##      12.13 seconds: SPR round   1 of   2, 7501 of 9060 nodes
##      12.25 seconds: SPR round   1 of   2, 7601 of 9060 nodes
##      12.41 seconds: SPR round   1 of   2, 7801 of 9060 nodes
##      12.55 seconds: SPR round   1 of   2, 8001 of 9060 nodes
##      12.69 seconds: SPR round   1 of   2, 8201 of 9060 nodes
##      12.83 seconds: SPR round   1 of   2, 8401 of 9060 nodes
##      13.00 seconds: SPR round   1 of   2, 8601 of 9060 nodes
##      13.14 seconds: SPR round   1 of   2, 8801 of 9060 nodes
##      13.30 seconds: SPR round   1 of   2, 9001 of 9060 nodes
##      13.40 seconds: ME NNI round 17 of 49, 1001 of 4529 splits, 7 changes (max delta 0.004)
##      13.50 seconds: ME NNI round 17 of 49, 3501 of 4529 splits, 32 changes (max delta 0.006)
##      13.60 seconds: ME NNI round 18 of 49, 1501 of 4529 splits, 4 changes (max delta 0.002)
##      13.70 seconds: ME NNI round 18 of 49, 3901 of 4529 splits, 17 changes (max delta 0.011)
##      13.87 seconds: SPR round   2 of   2, 101 of 9060 nodes
##      14.01 seconds: SPR round   2 of   2, 301 of 9060 nodes
##      14.17 seconds: SPR round   2 of   2, 501 of 9060 nodes
##      14.34 seconds: SPR round   2 of   2, 701 of 9060 nodes
##      14.50 seconds: SPR round   2 of   2, 901 of 9060 nodes
##      14.67 seconds: SPR round   2 of   2, 1101 of 9060 nodes
##      14.83 seconds: SPR round   2 of   2, 1301 of 9060 nodes
##      14.96 seconds: SPR round   2 of   2, 1501 of 9060 nodes
##      15.11 seconds: SPR round   2 of   2, 1701 of 9060 nodes
##      15.25 seconds: SPR round   2 of   2, 1901 of 9060 nodes
##      15.39 seconds: SPR round   2 of   2, 2101 of 9060 nodes
##      15.52 seconds: SPR round   2 of   2, 2301 of 9060 nodes
##      15.66 seconds: SPR round   2 of   2, 2501 of 9060 nodes
##      15.84 seconds: SPR round   2 of   2, 2701 of 9060 nodes
##      15.96 seconds: SPR round   2 of   2, 2901 of 9060 nodes
##      16.09 seconds: SPR round   2 of   2, 3101 of 9060 nodes
##      16.23 seconds: SPR round   2 of   2, 3301 of 9060 nodes
##      16.37 seconds: SPR round   2 of   2, 3501 of 9060 nodes
##      16.51 seconds: SPR round   2 of   2, 3701 of 9060 nodes
##      16.66 seconds: SPR round   2 of   2, 3901 of 9060 nodes
##      16.82 seconds: SPR round   2 of   2, 4101 of 9060 nodes
##      16.94 seconds: SPR round   2 of   2, 4301 of 9060 nodes
##      17.10 seconds: SPR round   2 of   2, 4501 of 9060 nodes
##      17.28 seconds: SPR round   2 of   2, 4701 of 9060 nodes
##      17.44 seconds: SPR round   2 of   2, 4901 of 9060 nodes
##      17.58 seconds: SPR round   2 of   2, 5101 of 9060 nodes
##      17.74 seconds: SPR round   2 of   2, 5301 of 9060 nodes
##      17.90 seconds: SPR round   2 of   2, 5501 of 9060 nodes
##      18.05 seconds: SPR round   2 of   2, 5701 of 9060 nodes
##      18.22 seconds: SPR round   2 of   2, 5901 of 9060 nodes
##      18.37 seconds: SPR round   2 of   2, 6101 of 9060 nodes
##      18.52 seconds: SPR round   2 of   2, 6301 of 9060 nodes
##      18.68 seconds: SPR round   2 of   2, 6501 of 9060 nodes
##      18.82 seconds: SPR round   2 of   2, 6701 of 9060 nodes
##      18.98 seconds: SPR round   2 of   2, 6901 of 9060 nodes
##      19.12 seconds: SPR round   2 of   2, 7101 of 9060 nodes
##      19.25 seconds: SPR round   2 of   2, 7301 of 9060 nodes
##      19.39 seconds: SPR round   2 of   2, 7501 of 9060 nodes
##      19.54 seconds: SPR round   2 of   2, 7701 of 9060 nodes
##      19.70 seconds: SPR round   2 of   2, 7901 of 9060 nodes
##      19.80 seconds: SPR round   2 of   2, 8101 of 9060 nodes
##      19.94 seconds: SPR round   2 of   2, 8301 of 9060 nodes
##      20.08 seconds: SPR round   2 of   2, 8501 of 9060 nodes
##      20.21 seconds: SPR round   2 of   2, 8701 of 9060 nodes
##      20.36 seconds: SPR round   2 of   2, 8901 of 9060 nodes
##      20.47 seconds: ME NNI round 33 of 49, 1 of 4529 splits
##      20.58 seconds: ME NNI round 33 of 49, 2601 of 4529 splits, 14 changes (max delta 0.074)
##      20.68 seconds: ME NNI round 34 of 49, 601 of 4529 splits, 0 changes
##      20.78 seconds: ME NNI round 34 of 49, 3101 of 4529 splits, 6 changes (max delta 0.017)
## Total branch-length 112.551 after 21.09 sec
##      21.11 seconds: ML Lengths 1 of 4529 splits
##      21.23 seconds: ML Lengths 301 of 4529 splits
##      21.34 seconds: ML Lengths 601 of 4529 splits
##      21.47 seconds: ML Lengths 901 of 4529 splits
##      21.59 seconds: ML Lengths 1201 of 4529 splits
##      21.70 seconds: ML Lengths 1501 of 4529 splits
##      21.83 seconds: ML Lengths 1801 of 4529 splits
##      21.94 seconds: ML Lengths 2101 of 4529 splits
##      22.06 seconds: ML Lengths 2401 of 4529 splits
##      22.18 seconds: ML Lengths 2701 of 4529 splits
##      22.29 seconds: ML Lengths 3001 of 4529 splits
##      22.40 seconds: ML Lengths 3301 of 4529 splits
##      22.51 seconds: ML Lengths 3601 of 4529 splits
##      22.62 seconds: ML Lengths 3901 of 4529 splits
##      22.73 seconds: ML Lengths 4201 of 4529 splits
##      22.84 seconds: ML Lengths 4501 of 4529 splits
##      22.98 seconds: ML NNI round 1 of 24, 101 of 4529 splits, 24 changes (max delta 3.666)
##      23.10 seconds: ML NNI round 1 of 24, 201 of 4529 splits, 41 changes (max delta 5.778)
##      23.23 seconds: ML NNI round 1 of 24, 301 of 4529 splits, 60 changes (max delta 9.178)
##      23.36 seconds: ML NNI round 1 of 24, 401 of 4529 splits, 67 changes (max delta 9.178)
##      23.50 seconds: ML NNI round 1 of 24, 501 of 4529 splits, 85 changes (max delta 10.374)
##      23.62 seconds: ML NNI round 1 of 24, 601 of 4529 splits, 102 changes (max delta 17.349)
##      23.76 seconds: ML NNI round 1 of 24, 701 of 4529 splits, 110 changes (max delta 17.349)
##      23.89 seconds: ML NNI round 1 of 24, 801 of 4529 splits, 129 changes (max delta 17.349)
##      24.03 seconds: ML NNI round 1 of 24, 901 of 4529 splits, 142 changes (max delta 17.349)
##      24.16 seconds: ML NNI round 1 of 24, 1001 of 4529 splits, 153 changes (max delta 17.349)
##      24.30 seconds: ML NNI round 1 of 24, 1101 of 4529 splits, 168 changes (max delta 17.349)
##      24.44 seconds: ML NNI round 1 of 24, 1201 of 4529 splits, 185 changes (max delta 17.349)
##      24.56 seconds: ML NNI round 1 of 24, 1301 of 4529 splits, 201 changes (max delta 17.349)
##      24.68 seconds: ML NNI round 1 of 24, 1401 of 4529 splits, 220 changes (max delta 17.349)
##      24.82 seconds: ML NNI round 1 of 24, 1501 of 4529 splits, 238 changes (max delta 17.349)
##      24.96 seconds: ML NNI round 1 of 24, 1601 of 4529 splits, 255 changes (max delta 17.349)
##      25.10 seconds: ML NNI round 1 of 24, 1701 of 4529 splits, 267 changes (max delta 17.349)
##      25.23 seconds: ML NNI round 1 of 24, 1801 of 4529 splits, 285 changes (max delta 17.349)
##      25.38 seconds: ML NNI round 1 of 24, 1901 of 4529 splits, 302 changes (max delta 17.349)
##      25.50 seconds: ML NNI round 1 of 24, 2001 of 4529 splits, 322 changes (max delta 17.349)
##      25.64 seconds: ML NNI round 1 of 24, 2101 of 4529 splits, 339 changes (max delta 17.349)
##      25.78 seconds: ML NNI round 1 of 24, 2201 of 4529 splits, 351 changes (max delta 17.349)
##      25.90 seconds: ML NNI round 1 of 24, 2301 of 4529 splits, 367 changes (max delta 17.349)
##      26.02 seconds: ML NNI round 1 of 24, 2401 of 4529 splits, 379 changes (max delta 17.349)
##      26.17 seconds: ML NNI round 1 of 24, 2501 of 4529 splits, 396 changes (max delta 17.349)
##      26.31 seconds: ML NNI round 1 of 24, 2601 of 4529 splits, 411 changes (max delta 17.349)
##      26.44 seconds: ML NNI round 1 of 24, 2701 of 4529 splits, 427 changes (max delta 17.349)
##      26.56 seconds: ML NNI round 1 of 24, 2801 of 4529 splits, 445 changes (max delta 17.349)
##      26.70 seconds: ML NNI round 1 of 24, 2901 of 4529 splits, 455 changes (max delta 17.349)
##      26.82 seconds: ML NNI round 1 of 24, 3001 of 4529 splits, 473 changes (max delta 17.349)
##      26.94 seconds: ML NNI round 1 of 24, 3101 of 4529 splits, 492 changes (max delta 17.349)
##      27.07 seconds: ML NNI round 1 of 24, 3201 of 4529 splits, 505 changes (max delta 17.349)
##      27.20 seconds: ML NNI round 1 of 24, 3301 of 4529 splits, 521 changes (max delta 17.875)
##      27.32 seconds: ML NNI round 1 of 24, 3401 of 4529 splits, 542 changes (max delta 17.875)
##      27.45 seconds: ML NNI round 1 of 24, 3501 of 4529 splits, 558 changes (max delta 17.875)
##      27.59 seconds: ML NNI round 1 of 24, 3601 of 4529 splits, 577 changes (max delta 17.875)
##      27.73 seconds: ML NNI round 1 of 24, 3701 of 4529 splits, 596 changes (max delta 17.875)
##      27.87 seconds: ML NNI round 1 of 24, 3801 of 4529 splits, 616 changes (max delta 17.875)
##      27.98 seconds: ML NNI round 1 of 24, 3901 of 4529 splits, 639 changes (max delta 17.875)
##      28.12 seconds: ML NNI round 1 of 24, 4001 of 4529 splits, 667 changes (max delta 17.875)
##      28.23 seconds: ML NNI round 1 of 24, 4101 of 4529 splits, 690 changes (max delta 17.875)
##      28.35 seconds: ML NNI round 1 of 24, 4201 of 4529 splits, 707 changes (max delta 17.875)
##      28.48 seconds: ML NNI round 1 of 24, 4301 of 4529 splits, 726 changes (max delta 17.875)
##      28.60 seconds: ML NNI round 1 of 24, 4401 of 4529 splits, 740 changes (max delta 17.875)
##      28.72 seconds: ML NNI round 1 of 24, 4501 of 4529 splits, 758 changes (max delta 17.875)
## ML-NNI round 1: LogLk = -321168.356 NNIs 765 max delta 17.88 Time 28.80
##      30.94 seconds: Optimizing GTR model, step 2 of 12
##      32.34 seconds: Optimizing GTR model, step 3 of 12
##      33.75 seconds: Optimizing GTR model, step 4 of 12
##      36.01 seconds: Optimizing GTR model, step 5 of 12
##      38.13 seconds: Optimizing GTR model, step 6 of 12
##      39.55 seconds: Optimizing GTR model, step 7 of 12
##      41.08 seconds: Optimizing GTR model, step 8 of 12
##      42.39 seconds: Optimizing GTR model, step 9 of 12
##      44.20 seconds: Optimizing GTR model, step 10 of 12
##      45.45 seconds: Optimizing GTR model, step 11 of 12
##      46.85 seconds: Optimizing GTR model, step 12 of 12
## GTR Frequencies: 0.2537 0.2145 0.3395 0.1924
## GTR rates(ac ag at cg ct gt) 0.6695 1.3983 1.1820 0.7451 2.9596 1.0000
##      48.40 seconds: ML Lengths 1 of 4529 splits
##      48.50 seconds: ML Lengths 201 of 4529 splits
##      48.61 seconds: ML Lengths 401 of 4529 splits
##      48.72 seconds: ML Lengths 601 of 4529 splits
##      48.83 seconds: ML Lengths 801 of 4529 splits
##      48.93 seconds: ML Lengths 1001 of 4529 splits
##      49.04 seconds: ML Lengths 1201 of 4529 splits
##      49.14 seconds: ML Lengths 1401 of 4529 splits
##      49.25 seconds: ML Lengths 1601 of 4529 splits
##      49.36 seconds: ML Lengths 1801 of 4529 splits
##      49.46 seconds: ML Lengths 2001 of 4529 splits
##      49.57 seconds: ML Lengths 2201 of 4529 splits
##      49.68 seconds: ML Lengths 2401 of 4529 splits
##      49.79 seconds: ML Lengths 2601 of 4529 splits
##      49.89 seconds: ML Lengths 2801 of 4529 splits
##      50.00 seconds: ML Lengths 3001 of 4529 splits
##      50.10 seconds: ML Lengths 3201 of 4529 splits
##      50.21 seconds: ML Lengths 3401 of 4529 splits
##      50.32 seconds: ML Lengths 3601 of 4529 splits
##      50.44 seconds: ML Lengths 3801 of 4529 splits
##      50.54 seconds: ML Lengths 4001 of 4529 splits
##      50.65 seconds: ML Lengths 4201 of 4529 splits
##      50.76 seconds: ML Lengths 4401 of 4529 splits
##      50.97 seconds: Site likelihoods with rate category 1 of 20
##      51.12 seconds: Site likelihoods with rate category 2 of 20
##      51.26 seconds: Site likelihoods with rate category 3 of 20
##      51.41 seconds: Site likelihoods with rate category 4 of 20
##      51.55 seconds: Site likelihoods with rate category 5 of 20
##      51.70 seconds: Site likelihoods with rate category 6 of 20
##      51.85 seconds: Site likelihoods with rate category 7 of 20
##      52.00 seconds: Site likelihoods with rate category 8 of 20
##      52.14 seconds: Site likelihoods with rate category 9 of 20
##      52.29 seconds: Site likelihoods with rate category 10 of 20
##      52.43 seconds: Site likelihoods with rate category 11 of 20
##      52.57 seconds: Site likelihoods with rate category 12 of 20
##      52.72 seconds: Site likelihoods with rate category 13 of 20
##      52.86 seconds: Site likelihoods with rate category 14 of 20
##      53.01 seconds: Site likelihoods with rate category 15 of 20
##      53.15 seconds: Site likelihoods with rate category 16 of 20
##      53.29 seconds: Site likelihoods with rate category 17 of 20
##      53.44 seconds: Site likelihoods with rate category 18 of 20
##      53.58 seconds: Site likelihoods with rate category 19 of 20
##      53.73 seconds: Site likelihoods with rate category 20 of 20
## Switched to using 20 rate categories (CAT approximation)
## Rate categories were divided by 1.145 so that average rate = 1.0
## CAT-based log-likelihoods may not be comparable across runs
## Use -gamma for approximate but comparable Gamma(20) log-likelihoods
##      53.96 seconds: ML NNI round 2 of 24, 1 of 4529 splits
##      54.14 seconds: ML NNI round 2 of 24, 101 of 4529 splits, 12 changes (max delta 3.472)
##      54.34 seconds: ML NNI round 2 of 24, 201 of 4529 splits, 31 changes (max delta 4.873)
##      54.51 seconds: ML NNI round 2 of 24, 301 of 4529 splits, 41 changes (max delta 5.772)
##      54.67 seconds: ML NNI round 2 of 24, 401 of 4529 splits, 51 changes (max delta 5.772)
##      54.83 seconds: ML NNI round 2 of 24, 501 of 4529 splits, 65 changes (max delta 6.603)
##      54.98 seconds: ML NNI round 2 of 24, 601 of 4529 splits, 73 changes (max delta 6.603)
##      55.13 seconds: ML NNI round 2 of 24, 701 of 4529 splits, 82 changes (max delta 6.603)
##      55.26 seconds: ML NNI round 2 of 24, 801 of 4529 splits, 93 changes (max delta 9.523)
##      55.44 seconds: ML NNI round 2 of 24, 901 of 4529 splits, 103 changes (max delta 9.523)
##      55.63 seconds: ML NNI round 2 of 24, 1001 of 4529 splits, 111 changes (max delta 9.523)
##      55.79 seconds: ML NNI round 2 of 24, 1101 of 4529 splits, 125 changes (max delta 9.523)
##      55.97 seconds: ML NNI round 2 of 24, 1201 of 4529 splits, 142 changes (max delta 9.523)
##      56.14 seconds: ML NNI round 2 of 24, 1301 of 4529 splits, 153 changes (max delta 9.523)
##      56.31 seconds: ML NNI round 2 of 24, 1401 of 4529 splits, 163 changes (max delta 9.523)
##      56.45 seconds: ML NNI round 2 of 24, 1501 of 4529 splits, 170 changes (max delta 9.523)
##      56.63 seconds: ML NNI round 2 of 24, 1601 of 4529 splits, 181 changes (max delta 9.523)
##      56.78 seconds: ML NNI round 2 of 24, 1701 of 4529 splits, 195 changes (max delta 9.523)
##      56.97 seconds: ML NNI round 2 of 24, 1801 of 4529 splits, 211 changes (max delta 9.523)
##      57.13 seconds: ML NNI round 2 of 24, 1901 of 4529 splits, 224 changes (max delta 9.523)
##      57.30 seconds: ML NNI round 2 of 24, 2001 of 4529 splits, 235 changes (max delta 9.523)
##      57.46 seconds: ML NNI round 2 of 24, 2101 of 4529 splits, 246 changes (max delta 9.523)
##      57.63 seconds: ML NNI round 2 of 24, 2201 of 4529 splits, 255 changes (max delta 9.523)
##      57.79 seconds: ML NNI round 2 of 24, 2301 of 4529 splits, 268 changes (max delta 9.523)
##      57.97 seconds: ML NNI round 2 of 24, 2401 of 4529 splits, 277 changes (max delta 9.523)
##      58.15 seconds: ML NNI round 2 of 24, 2501 of 4529 splits, 292 changes (max delta 9.523)
##      58.33 seconds: ML NNI round 2 of 24, 2601 of 4529 splits, 306 changes (max delta 9.523)
##      58.50 seconds: ML NNI round 2 of 24, 2701 of 4529 splits, 319 changes (max delta 9.523)
##      58.67 seconds: ML NNI round 2 of 24, 2801 of 4529 splits, 328 changes (max delta 9.523)
##      58.83 seconds: ML NNI round 2 of 24, 2901 of 4529 splits, 336 changes (max delta 9.523)
##      59.01 seconds: ML NNI round 2 of 24, 3001 of 4529 splits, 349 changes (max delta 9.523)
##      59.19 seconds: ML NNI round 2 of 24, 3101 of 4529 splits, 363 changes (max delta 9.523)
##      59.35 seconds: ML NNI round 2 of 24, 3201 of 4529 splits, 374 changes (max delta 9.523)
##      59.51 seconds: ML NNI round 2 of 24, 3301 of 4529 splits, 383 changes (max delta 9.523)
##      59.70 seconds: ML NNI round 2 of 24, 3401 of 4529 splits, 397 changes (max delta 9.523)
##      59.85 seconds: ML NNI round 2 of 24, 3501 of 4529 splits, 409 changes (max delta 9.523)
##      60.00 seconds: ML NNI round 2 of 24, 3601 of 4529 splits, 417 changes (max delta 9.523)
##      60.17 seconds: ML NNI round 2 of 24, 3701 of 4529 splits, 431 changes (max delta 9.523)
##      60.34 seconds: ML NNI round 2 of 24, 3801 of 4529 splits, 447 changes (max delta 9.523)
##      60.51 seconds: ML NNI round 2 of 24, 3901 of 4529 splits, 457 changes (max delta 9.523)
##      60.67 seconds: ML NNI round 2 of 24, 4001 of 4529 splits, 467 changes (max delta 9.523)
##      60.85 seconds: ML NNI round 2 of 24, 4101 of 4529 splits, 486 changes (max delta 9.523)
##      61.03 seconds: ML NNI round 2 of 24, 4201 of 4529 splits, 500 changes (max delta 9.523)
##      61.19 seconds: ML NNI round 2 of 24, 4301 of 4529 splits, 513 changes (max delta 9.523)
##      61.37 seconds: ML NNI round 2 of 24, 4401 of 4529 splits, 529 changes (max delta 9.523)
##      61.54 seconds: ML NNI round 2 of 24, 4501 of 4529 splits, 545 changes (max delta 9.523)
## ML-NNI round 2: LogLk = -259690.144 NNIs 550 max delta 9.52 Time 61.63
##      61.77 seconds: ML NNI round 3 of 24, 101 of 4529 splits, 7 changes (max delta 0.976)
##      61.94 seconds: ML NNI round 3 of 24, 201 of 4529 splits, 11 changes (max delta 0.976)
##      62.14 seconds: ML NNI round 3 of 24, 301 of 4529 splits, 20 changes (max delta 1.117)
##      62.32 seconds: ML NNI round 3 of 24, 401 of 4529 splits, 35 changes (max delta 6.231)
##      62.47 seconds: ML NNI round 3 of 24, 501 of 4529 splits, 44 changes (max delta 11.822)
##      62.64 seconds: ML NNI round 3 of 24, 601 of 4529 splits, 52 changes (max delta 11.822)
##      62.81 seconds: ML NNI round 3 of 24, 701 of 4529 splits, 60 changes (max delta 11.822)
##      62.99 seconds: ML NNI round 3 of 24, 801 of 4529 splits, 66 changes (max delta 11.822)
##      63.15 seconds: ML NNI round 3 of 24, 901 of 4529 splits, 75 changes (max delta 18.711)
##      63.29 seconds: ML NNI round 3 of 24, 1001 of 4529 splits, 76 changes (max delta 18.711)
##      63.46 seconds: ML NNI round 3 of 24, 1101 of 4529 splits, 82 changes (max delta 18.711)
##      63.61 seconds: ML NNI round 3 of 24, 1201 of 4529 splits, 84 changes (max delta 18.711)
##      63.78 seconds: ML NNI round 3 of 24, 1301 of 4529 splits, 92 changes (max delta 18.711)
##      63.94 seconds: ML NNI round 3 of 24, 1401 of 4529 splits, 102 changes (max delta 18.711)
##      64.11 seconds: ML NNI round 3 of 24, 1501 of 4529 splits, 111 changes (max delta 18.711)
##      64.27 seconds: ML NNI round 3 of 24, 1601 of 4529 splits, 114 changes (max delta 18.711)
##      64.42 seconds: ML NNI round 3 of 24, 1701 of 4529 splits, 118 changes (max delta 18.711)
##      64.58 seconds: ML NNI round 3 of 24, 1801 of 4529 splits, 121 changes (max delta 18.711)
##      64.73 seconds: ML NNI round 3 of 24, 1901 of 4529 splits, 124 changes (max delta 18.711)
##      64.90 seconds: ML NNI round 3 of 24, 2001 of 4529 splits, 133 changes (max delta 18.711)
##      65.04 seconds: ML NNI round 3 of 24, 2101 of 4529 splits, 135 changes (max delta 18.711)
##      65.21 seconds: ML NNI round 3 of 24, 2201 of 4529 splits, 138 changes (max delta 18.711)
##      65.37 seconds: ML NNI round 3 of 24, 2301 of 4529 splits, 145 changes (max delta 18.711)
##      65.54 seconds: ML NNI round 3 of 24, 2401 of 4529 splits, 153 changes (max delta 18.711)
##      65.71 seconds: ML NNI round 3 of 24, 2501 of 4529 splits, 161 changes (max delta 18.711)
##      65.88 seconds: ML NNI round 3 of 24, 2601 of 4529 splits, 173 changes (max delta 18.711)
##      66.04 seconds: ML NNI round 3 of 24, 2701 of 4529 splits, 178 changes (max delta 18.711)
##      66.22 seconds: ML NNI round 3 of 24, 2801 of 4529 splits, 184 changes (max delta 18.711)
##      66.39 seconds: ML NNI round 3 of 24, 2901 of 4529 splits, 190 changes (max delta 18.711)
##      66.55 seconds: ML NNI round 3 of 24, 3001 of 4529 splits, 195 changes (max delta 18.711)
##      66.68 seconds: ML NNI round 3 of 24, 3101 of 4529 splits, 196 changes (max delta 18.711)
##      66.82 seconds: ML NNI round 3 of 24, 3201 of 4529 splits, 201 changes (max delta 18.711)
## ML-NNI round 3: LogLk = -259414.801 NNIs 202 max delta 18.71 Time 66.90
##      67.05 seconds: ML NNI round 4 of 24, 101 of 4529 splits, 3 changes (max delta 0.533)
##      67.23 seconds: ML NNI round 4 of 24, 201 of 4529 splits, 12 changes (max delta 2.963)
##      67.38 seconds: ML NNI round 4 of 24, 301 of 4529 splits, 18 changes (max delta 4.361)
##      67.55 seconds: ML NNI round 4 of 24, 401 of 4529 splits, 21 changes (max delta 6.238)
##      67.72 seconds: ML NNI round 4 of 24, 501 of 4529 splits, 26 changes (max delta 6.238)
##      67.88 seconds: ML NNI round 4 of 24, 601 of 4529 splits, 34 changes (max delta 7.991)
##      68.03 seconds: ML NNI round 4 of 24, 701 of 4529 splits, 40 changes (max delta 7.991)
##      68.18 seconds: ML NNI round 4 of 24, 801 of 4529 splits, 44 changes (max delta 7.991)
##      68.33 seconds: ML NNI round 4 of 24, 901 of 4529 splits, 48 changes (max delta 7.991)
##      68.47 seconds: ML NNI round 4 of 24, 1001 of 4529 splits, 53 changes (max delta 7.991)
##      68.62 seconds: ML NNI round 4 of 24, 1101 of 4529 splits, 57 changes (max delta 7.991)
##      68.78 seconds: ML NNI round 4 of 24, 1201 of 4529 splits, 60 changes (max delta 7.991)
##      68.94 seconds: ML NNI round 4 of 24, 1301 of 4529 splits, 65 changes (max delta 7.991)
##      69.10 seconds: ML NNI round 4 of 24, 1401 of 4529 splits, 69 changes (max delta 7.991)
##      69.24 seconds: ML NNI round 4 of 24, 1501 of 4529 splits, 71 changes (max delta 7.991)
##      69.38 seconds: ML NNI round 4 of 24, 1601 of 4529 splits, 77 changes (max delta 7.991)
##      69.52 seconds: ML NNI round 4 of 24, 1701 of 4529 splits, 84 changes (max delta 7.991)
##      69.66 seconds: ML NNI round 4 of 24, 1801 of 4529 splits, 87 changes (max delta 7.991)
##      69.80 seconds: ML NNI round 4 of 24, 1901 of 4529 splits, 89 changes (max delta 7.991)
##      69.95 seconds: ML NNI round 4 of 24, 2001 of 4529 splits, 93 changes (max delta 7.991)
##      70.09 seconds: ML NNI round 4 of 24, 2101 of 4529 splits, 94 changes (max delta 7.991)
## ML-NNI round 4: LogLk = -259303.147 NNIs 97 max delta 7.99 Time 70.17
##      70.29 seconds: ML NNI round 5 of 24, 101 of 4529 splits, 3 changes (max delta 0.007)
##      70.43 seconds: ML NNI round 5 of 24, 201 of 4529 splits, 10 changes (max delta 6.690)
##      70.59 seconds: ML NNI round 5 of 24, 301 of 4529 splits, 12 changes (max delta 6.690)
##      70.74 seconds: ML NNI round 5 of 24, 401 of 4529 splits, 17 changes (max delta 6.690)
##      70.88 seconds: ML NNI round 5 of 24, 501 of 4529 splits, 24 changes (max delta 6.690)
##      71.04 seconds: ML NNI round 5 of 24, 601 of 4529 splits, 29 changes (max delta 6.690)
##      71.19 seconds: ML NNI round 5 of 24, 701 of 4529 splits, 33 changes (max delta 6.690)
##      71.36 seconds: ML NNI round 5 of 24, 801 of 4529 splits, 40 changes (max delta 6.690)
##      71.53 seconds: ML NNI round 5 of 24, 901 of 4529 splits, 48 changes (max delta 6.690)
## ML-NNI round 5: LogLk = -259252.340 NNIs 55 max delta 6.69 Time 71.69
##      71.69 seconds: ML NNI round 6 of 24, 1 of 4529 splits
##      71.81 seconds: ML NNI round 6 of 24, 101 of 4529 splits, 2 changes (max delta 3.025)
##      71.94 seconds: ML NNI round 6 of 24, 201 of 4529 splits, 5 changes (max delta 3.025)
##      72.09 seconds: ML NNI round 6 of 24, 301 of 4529 splits, 11 changes (max delta 3.025)
##      72.24 seconds: ML NNI round 6 of 24, 401 of 4529 splits, 12 changes (max delta 3.025)
##      72.39 seconds: ML NNI round 6 of 24, 501 of 4529 splits, 17 changes (max delta 6.444)
## ML-NNI round 6: LogLk = -259222.106 NNIs 19 max delta 6.44 Time 72.53
##      72.52 seconds: ML NNI round 7 of 24, 1 of 4529 splits
##      72.65 seconds: ML NNI round 7 of 24, 101 of 4529 splits, 4 changes (max delta 0.952)
##      72.80 seconds: ML NNI round 7 of 24, 201 of 4529 splits, 8 changes (max delta 0.952)
##      72.95 seconds: ML NNI round 7 of 24, 301 of 4529 splits, 10 changes (max delta 3.920)
## ML-NNI round 7: LogLk = -259214.213 NNIs 11 max delta 3.92 Time 73.04
##      73.17 seconds: ML NNI round 8 of 24, 101 of 4529 splits, 6 changes (max delta 6.079)
## ML-NNI round 8: LogLk = -259204.204 NNIs 9 max delta 6.08 Time 73.34
##      73.33 seconds: ML NNI round 9 of 24, 1 of 4529 splits
##      73.47 seconds: ML NNI round 9 of 24, 101 of 4529 splits, 4 changes (max delta 0.412)
## ML-NNI round 9: LogLk = -259203.411 NNIs 6 max delta 0.41 Time 73.56
## ML-NNI round 10: LogLk = -259203.400 NNIs 6 max delta 0.66 Time 73.72
## Turning off heuristics for final round of ML NNIs (converged)
##      73.72 seconds: ML NNI round 11 of 24, 1 of 4529 splits
##      73.94 seconds: ML NNI round 11 of 24, 101 of 4529 splits, 1 changes (max delta 0.879)
##      74.15 seconds: ML NNI round 11 of 24, 201 of 4529 splits, 5 changes (max delta 0.879)
##      74.38 seconds: ML NNI round 11 of 24, 301 of 4529 splits, 5 changes (max delta 0.879)
##      74.58 seconds: ML NNI round 11 of 24, 401 of 4529 splits, 9 changes (max delta 0.879)
##      74.80 seconds: ML NNI round 11 of 24, 501 of 4529 splits, 11 changes (max delta 5.384)
##      75.02 seconds: ML NNI round 11 of 24, 601 of 4529 splits, 14 changes (max delta 5.384)
##      75.24 seconds: ML NNI round 11 of 24, 701 of 4529 splits, 16 changes (max delta 5.384)
##      75.45 seconds: ML NNI round 11 of 24, 801 of 4529 splits, 17 changes (max delta 5.384)
##      75.66 seconds: ML NNI round 11 of 24, 901 of 4529 splits, 18 changes (max delta 5.384)
##      75.89 seconds: ML NNI round 11 of 24, 1001 of 4529 splits, 19 changes (max delta 5.384)
##      76.10 seconds: ML NNI round 11 of 24, 1101 of 4529 splits, 23 changes (max delta 5.384)
##      76.32 seconds: ML NNI round 11 of 24, 1201 of 4529 splits, 23 changes (max delta 5.384)
##      76.53 seconds: ML NNI round 11 of 24, 1301 of 4529 splits, 24 changes (max delta 5.384)
##      76.73 seconds: ML NNI round 11 of 24, 1401 of 4529 splits, 30 changes (max delta 5.384)
##      76.94 seconds: ML NNI round 11 of 24, 1501 of 4529 splits, 34 changes (max delta 5.384)
##      77.16 seconds: ML NNI round 11 of 24, 1601 of 4529 splits, 37 changes (max delta 5.384)
##      77.36 seconds: ML NNI round 11 of 24, 1701 of 4529 splits, 42 changes (max delta 5.384)
##      77.54 seconds: ML NNI round 11 of 24, 1801 of 4529 splits, 44 changes (max delta 5.384)
##      77.73 seconds: ML NNI round 11 of 24, 1901 of 4529 splits, 48 changes (max delta 5.384)
##      77.92 seconds: ML NNI round 11 of 24, 2001 of 4529 splits, 49 changes (max delta 5.384)
##      78.11 seconds: ML NNI round 11 of 24, 2101 of 4529 splits, 51 changes (max delta 5.384)
##      78.32 seconds: ML NNI round 11 of 24, 2201 of 4529 splits, 53 changes (max delta 5.384)
##      78.52 seconds: ML NNI round 11 of 24, 2301 of 4529 splits, 54 changes (max delta 5.384)
##      78.72 seconds: ML NNI round 11 of 24, 2401 of 4529 splits, 57 changes (max delta 5.384)
##      78.93 seconds: ML NNI round 11 of 24, 2501 of 4529 splits, 60 changes (max delta 5.384)
##      79.15 seconds: ML NNI round 11 of 24, 2601 of 4529 splits, 62 changes (max delta 5.384)
##      79.38 seconds: ML NNI round 11 of 24, 2701 of 4529 splits, 62 changes (max delta 5.384)
##      79.60 seconds: ML NNI round 11 of 24, 2801 of 4529 splits, 64 changes (max delta 5.384)
##      79.79 seconds: ML NNI round 11 of 24, 2901 of 4529 splits, 66 changes (max delta 5.384)
##      80.00 seconds: ML NNI round 11 of 24, 3001 of 4529 splits, 66 changes (max delta 5.384)
##      80.19 seconds: ML NNI round 11 of 24, 3101 of 4529 splits, 67 changes (max delta 5.384)
##      80.39 seconds: ML NNI round 11 of 24, 3201 of 4529 splits, 70 changes (max delta 5.384)
##      80.59 seconds: ML NNI round 11 of 24, 3301 of 4529 splits, 71 changes (max delta 5.384)
##      80.80 seconds: ML NNI round 11 of 24, 3401 of 4529 splits, 74 changes (max delta 5.384)
##      81.02 seconds: ML NNI round 11 of 24, 3501 of 4529 splits, 74 changes (max delta 5.384)
##      81.23 seconds: ML NNI round 11 of 24, 3601 of 4529 splits, 75 changes (max delta 5.384)
##      81.45 seconds: ML NNI round 11 of 24, 3701 of 4529 splits, 78 changes (max delta 5.384)
##      81.65 seconds: ML NNI round 11 of 24, 3801 of 4529 splits, 78 changes (max delta 5.384)
##      81.86 seconds: ML NNI round 11 of 24, 3901 of 4529 splits, 80 changes (max delta 5.384)
##      82.05 seconds: ML NNI round 11 of 24, 4001 of 4529 splits, 81 changes (max delta 5.384)
##      82.25 seconds: ML NNI round 11 of 24, 4101 of 4529 splits, 82 changes (max delta 5.384)
##      82.46 seconds: ML NNI round 11 of 24, 4201 of 4529 splits, 85 changes (max delta 5.384)
##      82.67 seconds: ML NNI round 11 of 24, 4301 of 4529 splits, 86 changes (max delta 5.384)
##      82.89 seconds: ML NNI round 11 of 24, 4401 of 4529 splits, 87 changes (max delta 5.384)
##      83.12 seconds: ML NNI round 11 of 24, 4501 of 4529 splits, 92 changes (max delta 5.384)
## ML-NNI round 11: LogLk = -258959.858 NNIs 93 max delta 5.38 Time 83.21 (final)
##      83.26 seconds: ML Lengths 101 of 4529 splits
##      83.37 seconds: ML Lengths 301 of 4529 splits
##      83.49 seconds: ML Lengths 501 of 4529 splits
##      83.61 seconds: ML Lengths 701 of 4529 splits
##      83.73 seconds: ML Lengths 901 of 4529 splits
##      83.85 seconds: ML Lengths 1101 of 4529 splits
##      83.96 seconds: ML Lengths 1301 of 4529 splits
##      84.07 seconds: ML Lengths 1501 of 4529 splits
##      84.18 seconds: ML Lengths 1701 of 4529 splits
##      84.29 seconds: ML Lengths 1901 of 4529 splits
##      84.39 seconds: ML Lengths 2101 of 4529 splits
##      84.51 seconds: ML Lengths 2301 of 4529 splits
##      84.62 seconds: ML Lengths 2501 of 4529 splits
##      84.73 seconds: ML Lengths 2701 of 4529 splits
##      84.84 seconds: ML Lengths 2901 of 4529 splits
##      84.95 seconds: ML Lengths 3101 of 4529 splits
##      85.06 seconds: ML Lengths 3301 of 4529 splits
##      85.17 seconds: ML Lengths 3501 of 4529 splits
##      85.29 seconds: ML Lengths 3701 of 4529 splits
##      85.39 seconds: ML Lengths 3901 of 4529 splits
##      85.50 seconds: ML Lengths 4101 of 4529 splits
##      85.61 seconds: ML Lengths 4301 of 4529 splits
##      85.73 seconds: ML Lengths 4501 of 4529 splits
## Optimize all lengths: LogLk = -258942.934 Time 85.78
##      86.04 seconds: ML split tests for    100 of   4528 internal splits
##      86.29 seconds: ML split tests for    200 of   4528 internal splits
##      86.54 seconds: ML split tests for    300 of   4528 internal splits
##      86.79 seconds: ML split tests for    400 of   4528 internal splits
##      87.04 seconds: ML split tests for    500 of   4528 internal splits
##      87.29 seconds: ML split tests for    600 of   4528 internal splits
##      87.54 seconds: ML split tests for    700 of   4528 internal splits
##      87.78 seconds: ML split tests for    800 of   4528 internal splits
##      88.02 seconds: ML split tests for    900 of   4528 internal splits
##      88.27 seconds: ML split tests for   1000 of   4528 internal splits
##      88.51 seconds: ML split tests for   1100 of   4528 internal splits
##      88.76 seconds: ML split tests for   1200 of   4528 internal splits
##      89.00 seconds: ML split tests for   1300 of   4528 internal splits
##      89.25 seconds: ML split tests for   1400 of   4528 internal splits
##      89.49 seconds: ML split tests for   1500 of   4528 internal splits
##      89.74 seconds: ML split tests for   1600 of   4528 internal splits
##      89.98 seconds: ML split tests for   1700 of   4528 internal splits
##      90.23 seconds: ML split tests for   1800 of   4528 internal splits
##      90.47 seconds: ML split tests for   1900 of   4528 internal splits
##      90.71 seconds: ML split tests for   2000 of   4528 internal splits
##      90.96 seconds: ML split tests for   2100 of   4528 internal splits
##      91.20 seconds: ML split tests for   2200 of   4528 internal splits
##      91.45 seconds: ML split tests for   2300 of   4528 internal splits
##      91.69 seconds: ML split tests for   2400 of   4528 internal splits
##      91.94 seconds: ML split tests for   2500 of   4528 internal splits
##      92.20 seconds: ML split tests for   2600 of   4528 internal splits
##      92.46 seconds: ML split tests for   2700 of   4528 internal splits
##      92.72 seconds: ML split tests for   2800 of   4528 internal splits
##      92.96 seconds: ML split tests for   2900 of   4528 internal splits
##      93.21 seconds: ML split tests for   3000 of   4528 internal splits
##      93.46 seconds: ML split tests for   3100 of   4528 internal splits
##      93.70 seconds: ML split tests for   3200 of   4528 internal splits
##      93.95 seconds: ML split tests for   3300 of   4528 internal splits
##      94.20 seconds: ML split tests for   3400 of   4528 internal splits
##      94.44 seconds: ML split tests for   3500 of   4528 internal splits
##      94.69 seconds: ML split tests for   3600 of   4528 internal splits
##      94.93 seconds: ML split tests for   3700 of   4528 internal splits
##      95.18 seconds: ML split tests for   3800 of   4528 internal splits
##      95.43 seconds: ML split tests for   3900 of   4528 internal splits
##      95.67 seconds: ML split tests for   4000 of   4528 internal splits
##      95.93 seconds: ML split tests for   4100 of   4528 internal splits
##      96.19 seconds: ML split tests for   4200 of   4528 internal splits
##      96.45 seconds: ML split tests for   4300 of   4528 internal splits
##      96.70 seconds: ML split tests for   4400 of   4528 internal splits
##      96.96 seconds: ML split tests for   4500 of   4528 internal splits
## Total time: 97.04 seconds Unique: 4531/4531 Bad splits: 20/4528 Worst delta-LogLk 1.754
## /local/workdir/sna49/moon_milk/moonmilk
## The working directory is /local/workdir/sna49/moon_milk/moonmilk
```



# Session Information

```r
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
##  package           * version    date (UTC) lib source
##  ade4                1.7-22     2023-02-06 [1] CRAN (R 4.3.2)
##  ape               * 5.7-1      2023-03-13 [2] CRAN (R 4.3.2)
##  aplot               0.2.2      2023-10-06 [1] CRAN (R 4.3.2)
##  Biobase             2.62.0     2023-10-24 [2] Bioconductor
##  BiocGenerics        0.48.1     2023-11-01 [2] Bioconductor
##  biomformat          1.30.0     2023-10-24 [1] Bioconductor
##  Biostrings          2.70.1     2023-10-25 [2] Bioconductor
##  bitops              1.0-7      2021-04-24 [2] CRAN (R 4.3.2)
##  bslib               0.5.1      2023-08-11 [2] CRAN (R 4.3.2)
##  cachem              1.0.8      2023-05-01 [2] CRAN (R 4.3.2)
##  callr               3.7.3      2022-11-02 [2] CRAN (R 4.3.2)
##  cli                 3.6.1      2023-03-23 [2] CRAN (R 4.3.2)
##  cluster             2.1.4      2022-08-22 [2] CRAN (R 4.3.2)
##  clusterGeneration   1.3.8      2023-08-16 [1] CRAN (R 4.3.2)
##  coda                0.19-4     2020-09-30 [2] CRAN (R 4.3.2)
##  codetools           0.2-19     2023-02-01 [2] CRAN (R 4.3.2)
##  colorspace          2.1-0      2023-01-23 [2] CRAN (R 4.3.2)
##  combinat            0.0-8      2012-10-29 [1] CRAN (R 4.3.2)
##  crayon              1.5.2      2022-09-29 [2] CRAN (R 4.3.2)
##  data.table          1.14.8     2023-02-17 [2] CRAN (R 4.3.2)
##  devtools            2.4.4      2022-07-20 [2] CRAN (R 4.2.1)
##  digest              0.6.33     2023-07-07 [2] CRAN (R 4.3.2)
##  doParallel          1.0.17     2022-02-07 [2] CRAN (R 4.3.2)
##  dplyr               1.1.3      2023-09-03 [2] CRAN (R 4.3.2)
##  ellipsis            0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
##  evaluate            0.23       2023-11-01 [2] CRAN (R 4.3.2)
##  expm                0.999-9    2024-01-11 [1] CRAN (R 4.3.2)
##  fansi               1.0.5      2023-10-08 [2] CRAN (R 4.3.2)
##  fastmap             1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
##  fastmatch           1.1-4      2023-08-18 [1] CRAN (R 4.3.2)
##  foreach             1.5.2      2022-02-02 [2] CRAN (R 4.3.2)
##  fs                  1.6.3      2023-07-20 [2] CRAN (R 4.3.2)
##  generics            0.1.3      2022-07-05 [2] CRAN (R 4.3.2)
##  GenomeInfoDb        1.38.0     2023-10-24 [2] Bioconductor
##  GenomeInfoDbData    1.2.11     2023-11-07 [2] Bioconductor
##  ggfun               0.1.4      2024-01-19 [1] CRAN (R 4.3.2)
##  ggplot2             3.5.0      2024-02-23 [2] CRAN (R 4.3.2)
##  ggplotify           0.1.2      2023-08-09 [1] CRAN (R 4.3.2)
##  ggtree            * 3.10.1     2024-02-25 [1] Bioconductor 3.18 (R 4.3.2)
##  glue                1.6.2      2022-02-24 [2] CRAN (R 4.3.2)
##  gridGraphics        0.5-1      2020-12-13 [1] CRAN (R 4.3.2)
##  gtable              0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
##  htmltools           0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
##  htmlwidgets         1.6.2      2023-03-17 [2] CRAN (R 4.3.2)
##  httpuv              1.6.12     2023-10-23 [2] CRAN (R 4.3.2)
##  igraph              1.5.1      2023-08-10 [2] CRAN (R 4.3.2)
##  IRanges             2.36.0     2023-10-24 [2] Bioconductor
##  iterators           1.0.14     2022-02-05 [2] CRAN (R 4.3.2)
##  jquerylib           0.1.4      2021-04-26 [2] CRAN (R 4.3.2)
##  jsonlite            1.8.7      2023-06-29 [2] CRAN (R 4.3.2)
##  knitr               1.45       2023-10-30 [2] CRAN (R 4.3.2)
##  later               1.3.1      2023-05-02 [2] CRAN (R 4.3.2)
##  lattice             0.21-9     2023-10-01 [2] CRAN (R 4.3.2)
##  lazyeval            0.2.2      2019-03-15 [2] CRAN (R 4.3.2)
##  lifecycle           1.0.3      2022-10-07 [2] CRAN (R 4.3.2)
##  magrittr            2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
##  maps              * 3.4.2      2023-12-15 [1] CRAN (R 4.3.2)
##  MASS                7.3-60     2023-05-04 [2] CRAN (R 4.3.2)
##  Matrix              1.6-1.1    2023-09-18 [2] CRAN (R 4.3.2)
##  memoise             2.0.1      2021-11-26 [2] CRAN (R 4.3.2)
##  mgcv                1.9-0      2023-07-11 [2] CRAN (R 4.3.2)
##  mime                0.12       2021-09-28 [2] CRAN (R 4.3.2)
##  miniUI              0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
##  mnormt              2.1.1      2022-09-26 [1] CRAN (R 4.3.2)
##  multtest            2.58.0     2023-10-24 [1] Bioconductor
##  munsell             0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
##  nlme                3.1-163    2023-08-09 [2] CRAN (R 4.3.2)
##  numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.2)
##  optimParallel       1.0-2      2021-02-11 [1] CRAN (R 4.3.2)
##  pacman              0.5.1      2019-03-11 [1] CRAN (R 4.3.2)
##  patchwork           1.2.0.9000 2024-03-11 [1] Github (thomasp85/patchwork@d943757)
##  permute             0.9-7      2022-01-27 [1] CRAN (R 4.3.2)
##  phangorn            2.11.1     2023-01-23 [1] CRAN (R 4.3.2)
##  phyloseq          * 1.46.0     2023-10-24 [1] Bioconductor
##  phytools          * 2.1-1      2024-01-09 [1] CRAN (R 4.3.2)
##  pillar              1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
##  pkgbuild            1.4.2      2023-06-26 [2] CRAN (R 4.3.2)
##  pkgconfig           2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
##  pkgload             1.3.3      2023-09-22 [2] CRAN (R 4.3.2)
##  plyr                1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
##  prettyunits         1.2.0      2023-09-24 [2] CRAN (R 4.3.2)
##  processx            3.8.2      2023-06-30 [2] CRAN (R 4.3.2)
##  profvis             0.3.8      2023-05-02 [2] CRAN (R 4.3.2)
##  promises            1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
##  ps                  1.7.5      2023-04-18 [2] CRAN (R 4.3.2)
##  purrr               1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
##  quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.3.2)
##  R6                  2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
##  RColorBrewer      * 1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
##  Rcpp                1.0.11     2023-07-06 [2] CRAN (R 4.3.2)
##  RCurl               1.98-1.14  2024-01-09 [1] CRAN (R 4.3.2)
##  remotes             2.4.2.1    2023-07-18 [2] CRAN (R 4.3.2)
##  reshape2            1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
##  rhdf5               2.46.1     2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
##  rhdf5filters        1.14.1     2023-11-06 [1] Bioconductor
##  Rhdf5lib            1.24.2     2024-02-07 [1] Bioconductor 3.18 (R 4.3.2)
##  rlang               1.1.2      2023-11-04 [2] CRAN (R 4.3.2)
##  rmarkdown           2.25       2023-09-18 [2] CRAN (R 4.3.2)
##  rstudioapi          0.15.0     2023-07-07 [2] CRAN (R 4.3.2)
##  S4Vectors           0.40.1     2023-10-26 [2] Bioconductor
##  sass                0.4.7      2023-07-15 [2] CRAN (R 4.3.2)
##  scales              1.3.0      2023-11-28 [2] CRAN (R 4.3.2)
##  scatterplot3d       0.3-44     2023-05-05 [1] CRAN (R 4.3.2)
##  sessioninfo         1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
##  shiny               1.7.5.1    2023-10-14 [2] CRAN (R 4.3.2)
##  stringi             1.7.12     2023-01-11 [2] CRAN (R 4.3.2)
##  stringr             1.5.0      2022-12-02 [2] CRAN (R 4.3.2)
##  survival            3.5-7      2023-08-14 [2] CRAN (R 4.3.2)
##  tibble              3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
##  tidyr               1.3.0      2023-01-24 [2] CRAN (R 4.3.2)
##  tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.3.2)
##  tidytree            0.4.6      2023-12-12 [1] CRAN (R 4.3.2)
##  treeio              1.26.0     2023-10-24 [1] Bioconductor
##  urlchecker          1.0.1      2021-11-30 [2] CRAN (R 4.3.2)
##  usethis             2.2.2      2023-07-06 [2] CRAN (R 4.3.2)
##  utf8                1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
##  vctrs               0.6.4      2023-10-12 [2] CRAN (R 4.3.2)
##  vegan               2.6-4      2022-10-11 [1] CRAN (R 4.3.2)
##  withr               2.5.2      2023-10-30 [2] CRAN (R 4.3.2)
##  xfun                0.41       2023-11-01 [2] CRAN (R 4.3.2)
##  xtable              1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
##  XVector             0.42.0     2023-10-24 [2] Bioconductor
##  yaml                2.3.7      2023-01-23 [2] CRAN (R 4.3.2)
##  yulab.utils         0.1.4      2024-01-28 [1] CRAN (R 4.3.2)
##  zlibbioc            1.48.0     2023-10-24 [2] Bioconductor
## 
##  [1] /home/sna49/R/x86_64-pc-linux-gnu-library/4.3
##  [2] /programs/R-4.3.2/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```
