---
title: "03a_Phylogenetics"
author: "sophiaaredas"
date: "2024-04-21"
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
load("data/02_PreProcessing/raw_preprocessed_physeq.RData")
raw_preprocessed_physeq
```

```
## Loading required package: phyloseq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2138 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 2138 taxa by 9 taxonomic ranks ]
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
## [2] "GACTACCGGGGTATCTAATCCCGTTTGCTCCCCTAGCTTTCGCGCCTCAGTGTCAGTGTCGGTCTAGGAAGCCGCCTTCGCCACCGGTGTTCCTCCGGATATCTACGCATTTCACCGCTACACCCGGAATTCCGCTTCCCTCTCCCGAACTCTAGCCCGCCAGTTTCCCGTGCCATTCCTCAGTTGAGCCGAGGGCTTTCACACGGGACTTAGCGGACCACCTACGCGCCCTTTACGCCCAGTAATTCCGAACAACGCTTGCCACCTCTGTATTACCGCGGCTGCTGGCACAGAGTTAGCCGTGGCTTCCTCCACCGGTACAGTCAATAGCCCGGCCTGTTCAGCCGTTCTACATTCGTCCCGGTCGACAGGGGTTTACGATCCGAAGACCTTCATCCCCCACGCGGCGTTGCTTCGTCAGGGTTTCCCCCATTGCGAAAAATTCCCCACTGCAGCCCCCCGTAGG"
## [3] ">ASV_2"                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
## [4] "CCTACGGGGGGCTGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCATTAACCTAATACGTTAGTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTAATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCAGTAGTC" 
## [5] ">ASV_3"                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
## [6] "CCTACGGGGGGCTGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCACTTTAAGTTGGGAGGAAGGGCATTAACCTAATACGTTAGTGTTTTGACGTTACCGACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTAATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCGGTAGTC"
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
    1 / 2138
  101 / 2138
  201 / 2138
  301 / 2138
  401 / 2138
  501 / 2138
  601 / 2138
  701 / 2138
  801 / 2138
  901 / 2138
 1001 / 2138
 1101 / 2138
 1201 / 2138
 1301 / 2138
 1401 / 2138
 1501 / 2138
 1601 / 2138
 1701 / 2138
 1801 / 2138
 1901 / 2138
 2001 / 2138
 2101 / 2138
## done.
## 
## Constructing a UPGMA tree (efffree=0) ... 
## 
    0 / 2138
   10 / 2138
   20 / 2138
   30 / 2138
   40 / 2138
   50 / 2138
   60 / 2138
   70 / 2138
   80 / 2138
   90 / 2138
  100 / 2138
  110 / 2138
  120 / 2138
  130 / 2138
  140 / 2138
  150 / 2138
  160 / 2138
  170 / 2138
  180 / 2138
  190 / 2138
  200 / 2138
  210 / 2138
  220 / 2138
  230 / 2138
  240 / 2138
  250 / 2138
  260 / 2138
  270 / 2138
  280 / 2138
  290 / 2138
  300 / 2138
  310 / 2138
  320 / 2138
  330 / 2138
  340 / 2138
  350 / 2138
  360 / 2138
  370 / 2138
  380 / 2138
  390 / 2138
  400 / 2138
  410 / 2138
  420 / 2138
  430 / 2138
  440 / 2138
  450 / 2138
  460 / 2138
  470 / 2138
  480 / 2138
  490 / 2138
  500 / 2138
  510 / 2138
  520 / 2138
  530 / 2138
  540 / 2138
  550 / 2138
  560 / 2138
  570 / 2138
  580 / 2138
  590 / 2138
  600 / 2138
  610 / 2138
  620 / 2138
  630 / 2138
  640 / 2138
  650 / 2138
  660 / 2138
  670 / 2138
  680 / 2138
  690 / 2138
  700 / 2138
  710 / 2138
  720 / 2138
  730 / 2138
  740 / 2138
  750 / 2138
  760 / 2138
  770 / 2138
  780 / 2138
  790 / 2138
  800 / 2138
  810 / 2138
  820 / 2138
  830 / 2138
  840 / 2138
  850 / 2138
  860 / 2138
  870 / 2138
  880 / 2138
  890 / 2138
  900 / 2138
  910 / 2138
  920 / 2138
  930 / 2138
  940 / 2138
  950 / 2138
  960 / 2138
  970 / 2138
  980 / 2138
  990 / 2138
 1000 / 2138
 1010 / 2138
 1020 / 2138
 1030 / 2138
 1040 / 2138
 1050 / 2138
 1060 / 2138
 1070 / 2138
 1080 / 2138
 1090 / 2138
 1100 / 2138
 1110 / 2138
 1120 / 2138
 1130 / 2138
 1140 / 2138
 1150 / 2138
 1160 / 2138
 1170 / 2138
 1180 / 2138
 1190 / 2138
 1200 / 2138
 1210 / 2138
 1220 / 2138
 1230 / 2138
 1240 / 2138
 1250 / 2138
 1260 / 2138
 1270 / 2138
 1280 / 2138
 1290 / 2138
 1300 / 2138
 1310 / 2138
 1320 / 2138
 1330 / 2138
 1340 / 2138
 1350 / 2138
 1360 / 2138
 1370 / 2138
 1380 / 2138
 1390 / 2138
 1400 / 2138
 1410 / 2138
 1420 / 2138
 1430 / 2138
 1440 / 2138
 1450 / 2138
 1460 / 2138
 1470 / 2138
 1480 / 2138
 1490 / 2138
 1500 / 2138
 1510 / 2138
 1520 / 2138
 1530 / 2138
 1540 / 2138
 1550 / 2138
 1560 / 2138
 1570 / 2138
 1580 / 2138
 1590 / 2138
 1600 / 2138
 1610 / 2138
 1620 / 2138
 1630 / 2138
 1640 / 2138
 1650 / 2138
 1660 / 2138
 1670 / 2138
 1680 / 2138
 1690 / 2138
 1700 / 2138
 1710 / 2138
 1720 / 2138
 1730 / 2138
 1740 / 2138
 1750 / 2138
 1760 / 2138
 1770 / 2138
 1780 / 2138
 1790 / 2138
 1800 / 2138
 1810 / 2138
 1820 / 2138
 1830 / 2138
 1840 / 2138
 1850 / 2138
 1860 / 2138
 1870 / 2138
 1880 / 2138
 1890 / 2138
 1900 / 2138
 1910 / 2138
 1920 / 2138
 1930 / 2138
 1940 / 2138
 1950 / 2138
 1960 / 2138
 1970 / 2138
 1980 / 2138
 1990 / 2138
 2000 / 2138
 2010 / 2138
 2020 / 2138
 2030 / 2138
 2040 / 2138
 2050 / 2138
 2060 / 2138
 2070 / 2138
 2080 / 2138
 2090 / 2138
 2100 / 2138
 2110 / 2138
 2120 / 2138
 2130 / 2138
## done.
## 
## Progressive alignment 1/2... 
## 
STEP     1 / 2137  f
STEP     2 / 2137  f
STEP     3 / 2137  f
STEP     4 / 2137  f
STEP     5 / 2137  f
STEP     6 / 2137  f
STEP     7 / 2137  f
STEP     8 / 2137  f
STEP     9 / 2137  f
STEP    10 / 2137  f
STEP    11 / 2137  f
STEP    12 / 2137  f
STEP    13 / 2137  f
STEP    14 / 2137  f
STEP    15 / 2137  f
STEP    16 / 2137  f
STEP    17 / 2137  f
STEP    18 / 2137  f
STEP    19 / 2137  f
STEP    20 / 2137  f
STEP    21 / 2137  f
STEP    22 / 2137  f
STEP    23 / 2137  f
STEP    24 / 2137  f
STEP    25 / 2137  f
STEP    26 / 2137  f
STEP    27 / 2137  f
STEP    28 / 2137  f
STEP    29 / 2137  f
STEP    30 / 2137  f
STEP    31 / 2137  f
STEP    32 / 2137  f
STEP    33 / 2137  f
STEP    34 / 2137  f
STEP    35 / 2137  f
STEP    36 / 2137  f
STEP    37 / 2137  f
STEP    38 / 2137  f
STEP    39 / 2137  f
STEP    40 / 2137  f
STEP    41 / 2137  f
STEP    42 / 2137  f
STEP    43 / 2137  f
STEP    44 / 2137  f
STEP    45 / 2137  f
STEP    46 / 2137  f
STEP    47 / 2137  f
STEP    48 / 2137  f
STEP    49 / 2137  f
STEP    50 / 2137  f
STEP    51 / 2137  f
STEP    52 / 2137  f
STEP    53 / 2137  f
STEP    54 / 2137  f
STEP    55 / 2137  f
STEP    56 / 2137  f
STEP    57 / 2137  f
STEP    58 / 2137  f
STEP    59 / 2137  f
STEP    60 / 2137  f
STEP    61 / 2137  f
STEP    62 / 2137  f
STEP    63 / 2137  f
STEP    64 / 2137  f
STEP    65 / 2137  f
STEP    66 / 2137  f
STEP    67 / 2137  f
STEP    68 / 2137  f
STEP    69 / 2137  f
STEP    70 / 2137  f
STEP    71 / 2137  f
STEP    72 / 2137  f
STEP    73 / 2137  f
STEP    74 / 2137  f
STEP    75 / 2137  f
STEP    76 / 2137  f
STEP    77 / 2137  f
STEP    78 / 2137  f
STEP    79 / 2137  f
STEP    80 / 2137  f
STEP    81 / 2137  f
STEP    82 / 2137  f
STEP    83 / 2137  f
STEP    84 / 2137  f
STEP    85 / 2137  f
STEP    86 / 2137  f
STEP    87 / 2137  f
STEP    88 / 2137  f
STEP    89 / 2137  f
STEP    90 / 2137  f
STEP    91 / 2137  f
STEP    92 / 2137  f
STEP    93 / 2137  f
STEP    94 / 2137  f
STEP    95 / 2137  f
STEP    96 / 2137  f
STEP    97 / 2137  f
STEP    98 / 2137  f
STEP    99 / 2137  f
STEP   100 / 2137  f
STEP   101 / 2137  f
STEP   102 / 2137  f
STEP   103 / 2137  f
STEP   104 / 2137  f
STEP   105 / 2137  f
STEP   106 / 2137  f
STEP   107 / 2137  f
STEP   108 / 2137  f
STEP   109 / 2137  f
STEP   110 / 2137  f
STEP   111 / 2137  f
STEP   112 / 2137  f
STEP   113 / 2137  f
STEP   114 / 2137  f
STEP   115 / 2137  f
STEP   116 / 2137  f
STEP   117 / 2137  f
STEP   118 / 2137  f
STEP   119 / 2137  f
STEP   120 / 2137  f
STEP   121 / 2137  f
STEP   122 / 2137  f
STEP   123 / 2137  f
STEP   124 / 2137  f
STEP   125 / 2137  f
STEP   126 / 2137  f
STEP   127 / 2137  f
STEP   128 / 2137  f
STEP   129 / 2137  f
STEP   130 / 2137  f
STEP   131 / 2137  f
STEP   132 / 2137  f
STEP   133 / 2137  f
STEP   134 / 2137  f
STEP   135 / 2137  f
STEP   136 / 2137  f
STEP   137 / 2137  f
STEP   138 / 2137  f
STEP   139 / 2137  f
STEP   140 / 2137  f
STEP   141 / 2137  f
STEP   142 / 2137  f
STEP   143 / 2137  f
STEP   144 / 2137  f
STEP   145 / 2137  f
STEP   146 / 2137  f
STEP   147 / 2137  f
STEP   148 / 2137  f
STEP   149 / 2137  f
STEP   150 / 2137  f
STEP   151 / 2137  f
STEP   152 / 2137  f
STEP   153 / 2137  f
STEP   154 / 2137  f
STEP   155 / 2137  f
STEP   156 / 2137  f
STEP   157 / 2137  f
STEP   158 / 2137  f
STEP   159 / 2137  f
STEP   160 / 2137  f
STEP   161 / 2137  f
STEP   162 / 2137  f
STEP   163 / 2137  f
STEP   164 / 2137  f
STEP   165 / 2137  f
STEP   166 / 2137  f
STEP   167 / 2137  f
STEP   168 / 2137  f
STEP   169 / 2137  f
STEP   170 / 2137  f
STEP   171 / 2137  f
STEP   172 / 2137  f
STEP   173 / 2137  f
STEP   174 / 2137  f
STEP   175 / 2137  f
STEP   176 / 2137  f
STEP   177 / 2137  f
STEP   178 / 2137  f
STEP   179 / 2137  f
STEP   180 / 2137  f
STEP   181 / 2137  f
STEP   182 / 2137  f
STEP   183 / 2137  f
STEP   184 / 2137  f
STEP   185 / 2137  f
STEP   186 / 2137  f
STEP   187 / 2137  f
STEP   188 / 2137  f
STEP   189 / 2137  f
STEP   190 / 2137  f
STEP   191 / 2137  f
STEP   192 / 2137  f
STEP   193 / 2137  f
STEP   194 / 2137  f
STEP   195 / 2137  f
STEP   196 / 2137  f
STEP   197 / 2137  f
STEP   198 / 2137  f
STEP   199 / 2137  f
STEP   200 / 2137  f
STEP   201 / 2137  f
STEP   202 / 2137  f
STEP   203 / 2137  f
STEP   204 / 2137  f
STEP   205 / 2137  f
STEP   206 / 2137  f
STEP   207 / 2137  f
STEP   208 / 2137  f
STEP   209 / 2137  f
STEP   210 / 2137  f
STEP   211 / 2137  f
STEP   212 / 2137  f
STEP   213 / 2137  f
STEP   214 / 2137  f
STEP   215 / 2137  f
STEP   216 / 2137  f
STEP   217 / 2137  f
STEP   218 / 2137  f
STEP   219 / 2137  f
STEP   220 / 2137  f
STEP   221 / 2137  f
STEP   222 / 2137  f
STEP   223 / 2137  f
STEP   224 / 2137  f
STEP   225 / 2137  f
STEP   226 / 2137  f
STEP   227 / 2137  f
STEP   228 / 2137  f
STEP   229 / 2137  f
STEP   230 / 2137  f
STEP   231 / 2137  f
STEP   232 / 2137  f
STEP   233 / 2137  f
STEP   234 / 2137  f
STEP   235 / 2137  f
STEP   236 / 2137  f
STEP   237 / 2137  f
STEP   238 / 2137  f
STEP   239 / 2137  f
STEP   240 / 2137  f
STEP   241 / 2137  f
STEP   242 / 2137  f
STEP   243 / 2137  f
STEP   244 / 2137  f
STEP   245 / 2137  f
STEP   246 / 2137  f
STEP   247 / 2137  f
STEP   248 / 2137  f
STEP   249 / 2137  f
STEP   250 / 2137  f
STEP   251 / 2137  f
STEP   252 / 2137  f
STEP   253 / 2137  f
STEP   254 / 2137  f
STEP   255 / 2137  f
STEP   256 / 2137  f
STEP   257 / 2137  f
STEP   258 / 2137  f
STEP   259 / 2137  f
STEP   260 / 2137  f
STEP   261 / 2137  f
STEP   262 / 2137  f
STEP   263 / 2137  f
STEP   264 / 2137  f
STEP   265 / 2137  f
STEP   266 / 2137  f
STEP   267 / 2137  f
STEP   268 / 2137  f
STEP   269 / 2137  f
STEP   270 / 2137  f
STEP   271 / 2137  f
STEP   272 / 2137  f
STEP   273 / 2137  f
STEP   274 / 2137  f
STEP   275 / 2137  f
STEP   276 / 2137  f
STEP   277 / 2137  f
STEP   278 / 2137  f
STEP   279 / 2137  f
STEP   280 / 2137  f
STEP   281 / 2137  f
STEP   282 / 2137  f
STEP   283 / 2137  f
STEP   284 / 2137  f
STEP   285 / 2137  f
STEP   286 / 2137  f
STEP   287 / 2137  f
STEP   288 / 2137  f
STEP   289 / 2137  f
STEP   290 / 2137  f
STEP   291 / 2137  f
STEP   292 / 2137  f
STEP   293 / 2137  f
STEP   294 / 2137  f
STEP   295 / 2137  f
STEP   296 / 2137  f
STEP   297 / 2137  f
STEP   298 / 2137  f
STEP   299 / 2137  f
STEP   300 / 2137  f
STEP   301 / 2137  f
STEP   302 / 2137  f
STEP   303 / 2137  f
STEP   304 / 2137  f
STEP   305 / 2137  f
STEP   306 / 2137  f
STEP   307 / 2137  f
STEP   308 / 2137  f
STEP   309 / 2137  f
STEP   310 / 2137  f
STEP   311 / 2137  f
STEP   312 / 2137  f
STEP   313 / 2137  f
STEP   314 / 2137  f
STEP   315 / 2137  f
STEP   316 / 2137  f
STEP   317 / 2137  f
STEP   318 / 2137  f
STEP   319 / 2137  f
STEP   320 / 2137  f
STEP   321 / 2137  f
STEP   322 / 2137  f
STEP   323 / 2137  f
STEP   324 / 2137  f
STEP   325 / 2137  f
STEP   326 / 2137  f
STEP   327 / 2137  f
STEP   328 / 2137  f
STEP   329 / 2137  f
STEP   330 / 2137  f
STEP   331 / 2137  f
STEP   332 / 2137  f
STEP   333 / 2137  f
STEP   334 / 2137  f
STEP   335 / 2137  f
STEP   336 / 2137  f
STEP   337 / 2137  f
STEP   338 / 2137  f
STEP   339 / 2137  f
STEP   340 / 2137  f
STEP   341 / 2137  f
STEP   342 / 2137  f
STEP   343 / 2137  f
STEP   344 / 2137  f
STEP   345 / 2137  f
STEP   346 / 2137  f
STEP   347 / 2137  f
STEP   348 / 2137  f
STEP   349 / 2137  f
STEP   350 / 2137  f
STEP   351 / 2137  f
STEP   352 / 2137  f
STEP   353 / 2137  f
STEP   354 / 2137  f
STEP   355 / 2137  f
STEP   356 / 2137  f
STEP   357 / 2137  f
STEP   358 / 2137  f
STEP   359 / 2137  f
STEP   360 / 2137  f
STEP   361 / 2137  f
STEP   362 / 2137  f
STEP   363 / 2137  f
STEP   364 / 2137  f
STEP   365 / 2137  f
STEP   366 / 2137  f
STEP   367 / 2137  f
STEP   368 / 2137  f
STEP   369 / 2137  f
STEP   370 / 2137  f
STEP   371 / 2137  f
STEP   372 / 2137  f
STEP   373 / 2137  f
STEP   374 / 2137  f
STEP   375 / 2137  f
STEP   376 / 2137  f
STEP   377 / 2137  f
STEP   378 / 2137  f
STEP   379 / 2137  f
STEP   380 / 2137  f
STEP   381 / 2137  f
STEP   382 / 2137  f
STEP   383 / 2137  f
STEP   384 / 2137  f
STEP   385 / 2137  f
STEP   386 / 2137  f
STEP   387 / 2137  f
STEP   388 / 2137  f
STEP   389 / 2137  f
STEP   390 / 2137  f
STEP   391 / 2137  f
STEP   392 / 2137  f
STEP   393 / 2137  f
STEP   394 / 2137  f
STEP   395 / 2137  f
STEP   396 / 2137  f
STEP   397 / 2137  f
STEP   398 / 2137  f
STEP   399 / 2137  f
STEP   400 / 2137  f
STEP   401 / 2137  f
STEP   402 / 2137  f
STEP   403 / 2137  f
STEP   404 / 2137  f
STEP   405 / 2137  f
STEP   406 / 2137  f
STEP   407 / 2137  f
STEP   408 / 2137  f
STEP   409 / 2137  f
STEP   410 / 2137  f
STEP   411 / 2137  f
STEP   412 / 2137  f
STEP   413 / 2137  f
STEP   414 / 2137  f
STEP   415 / 2137  f
STEP   416 / 2137  f
STEP   417 / 2137  f
STEP   418 / 2137  f
STEP   419 / 2137  f
STEP   420 / 2137  f
STEP   421 / 2137  f
STEP   422 / 2137  f
STEP   423 / 2137  f
STEP   424 / 2137  f
STEP   425 / 2137  f
STEP   426 / 2137  f
STEP   427 / 2137  f
STEP   428 / 2137  f
STEP   429 / 2137  f
STEP   430 / 2137  f
STEP   431 / 2137  f
STEP   432 / 2137  f
STEP   433 / 2137  f
STEP   434 / 2137  f
STEP   435 / 2137  f
STEP   436 / 2137  f
STEP   437 / 2137  f
STEP   438 / 2137  f
STEP   439 / 2137  f
STEP   440 / 2137  f
STEP   441 / 2137  f
STEP   442 / 2137  f
STEP   443 / 2137  f
STEP   444 / 2137  f
STEP   445 / 2137  f
STEP   446 / 2137  f
STEP   447 / 2137  f
STEP   448 / 2137  f
STEP   449 / 2137  f
STEP   450 / 2137  f
STEP   451 / 2137  f
STEP   452 / 2137  f
STEP   453 / 2137  f
STEP   454 / 2137  f
STEP   455 / 2137  f
STEP   456 / 2137  f
STEP   457 / 2137  f
STEP   458 / 2137  f
STEP   459 / 2137  f
STEP   460 / 2137  f
STEP   461 / 2137  f
STEP   462 / 2137  f
STEP   463 / 2137  f
STEP   464 / 2137  f
STEP   465 / 2137  f
STEP   466 / 2137  f
STEP   467 / 2137  f
STEP   468 / 2137  f
STEP   469 / 2137  f
STEP   470 / 2137  f
STEP   471 / 2137  f
STEP   472 / 2137  f
STEP   473 / 2137  f
STEP   474 / 2137  f
STEP   475 / 2137  f
STEP   476 / 2137  f
STEP   477 / 2137  f
STEP   478 / 2137  f
STEP   479 / 2137  f
STEP   480 / 2137  f
STEP   481 / 2137  f
STEP   482 / 2137  f
STEP   483 / 2137  f
STEP   484 / 2137  f
STEP   485 / 2137  f
STEP   486 / 2137  f
STEP   487 / 2137  f
STEP   488 / 2137  f
STEP   489 / 2137  f
STEP   490 / 2137  f
STEP   491 / 2137  f
STEP   492 / 2137  f
STEP   493 / 2137  f
STEP   494 / 2137  f
STEP   495 / 2137  f
STEP   496 / 2137  f
STEP   497 / 2137  f
STEP   498 / 2137  f
STEP   499 / 2137  f
STEP   500 / 2137  f
STEP   501 / 2137  f
STEP   601 / 2137  f
STEP   701 / 2137  f
STEP   801 / 2137  f
STEP   901 / 2137  f
STEP  1001 / 2137  f
STEP  1101 / 2137  f
STEP  1201 / 2137  f
STEP  1301 / 2137  f
STEP  1401 / 2137  f
STEP  1501 / 2137  f
STEP  1601 / 2137  f
## Reallocating..done. *alloclen = 1934
## 
STEP  1701 / 2137  f
STEP  1801 / 2137  f
STEP  1901 / 2137  f
STEP  2001 / 2137  f
STEP  2101 / 2137  f
## done.
## 
## Making a distance matrix from msa.. 
## 
    0 / 2138
  100 / 2138
  200 / 2138
  300 / 2138
  400 / 2138
  500 / 2138
  600 / 2138
  700 / 2138
  800 / 2138
  900 / 2138
 1000 / 2138
 1100 / 2138
 1200 / 2138
 1300 / 2138
 1400 / 2138
 1500 / 2138
 1600 / 2138
 1700 / 2138
 1800 / 2138
 1900 / 2138
 2000 / 2138
 2100 / 2138
## done.
## 
## Constructing a UPGMA tree (efffree=1) ... 
## 
    0 / 2138
   10 / 2138
   20 / 2138
   30 / 2138
   40 / 2138
   50 / 2138
   60 / 2138
   70 / 2138
   80 / 2138
   90 / 2138
  100 / 2138
  110 / 2138
  120 / 2138
  130 / 2138
  140 / 2138
  150 / 2138
  160 / 2138
  170 / 2138
  180 / 2138
  190 / 2138
  200 / 2138
  210 / 2138
  220 / 2138
  230 / 2138
  240 / 2138
  250 / 2138
  260 / 2138
  270 / 2138
  280 / 2138
  290 / 2138
  300 / 2138
  310 / 2138
  320 / 2138
  330 / 2138
  340 / 2138
  350 / 2138
  360 / 2138
  370 / 2138
  380 / 2138
  390 / 2138
  400 / 2138
  410 / 2138
  420 / 2138
  430 / 2138
  440 / 2138
  450 / 2138
  460 / 2138
  470 / 2138
  480 / 2138
  490 / 2138
  500 / 2138
  510 / 2138
  520 / 2138
  530 / 2138
  540 / 2138
  550 / 2138
  560 / 2138
  570 / 2138
  580 / 2138
  590 / 2138
  600 / 2138
  610 / 2138
  620 / 2138
  630 / 2138
  640 / 2138
  650 / 2138
  660 / 2138
  670 / 2138
  680 / 2138
  690 / 2138
  700 / 2138
  710 / 2138
  720 / 2138
  730 / 2138
  740 / 2138
  750 / 2138
  760 / 2138
  770 / 2138
  780 / 2138
  790 / 2138
  800 / 2138
  810 / 2138
  820 / 2138
  830 / 2138
  840 / 2138
  850 / 2138
  860 / 2138
  870 / 2138
  880 / 2138
  890 / 2138
  900 / 2138
  910 / 2138
  920 / 2138
  930 / 2138
  940 / 2138
  950 / 2138
  960 / 2138
  970 / 2138
  980 / 2138
  990 / 2138
 1000 / 2138
 1010 / 2138
 1020 / 2138
 1030 / 2138
 1040 / 2138
 1050 / 2138
 1060 / 2138
 1070 / 2138
 1080 / 2138
 1090 / 2138
 1100 / 2138
 1110 / 2138
 1120 / 2138
 1130 / 2138
 1140 / 2138
 1150 / 2138
 1160 / 2138
 1170 / 2138
 1180 / 2138
 1190 / 2138
 1200 / 2138
 1210 / 2138
 1220 / 2138
 1230 / 2138
 1240 / 2138
 1250 / 2138
 1260 / 2138
 1270 / 2138
 1280 / 2138
 1290 / 2138
 1300 / 2138
 1310 / 2138
 1320 / 2138
 1330 / 2138
 1340 / 2138
 1350 / 2138
 1360 / 2138
 1370 / 2138
 1380 / 2138
 1390 / 2138
 1400 / 2138
 1410 / 2138
 1420 / 2138
 1430 / 2138
 1440 / 2138
 1450 / 2138
 1460 / 2138
 1470 / 2138
 1480 / 2138
 1490 / 2138
 1500 / 2138
 1510 / 2138
 1520 / 2138
 1530 / 2138
 1540 / 2138
 1550 / 2138
 1560 / 2138
 1570 / 2138
 1580 / 2138
 1590 / 2138
 1600 / 2138
 1610 / 2138
 1620 / 2138
 1630 / 2138
 1640 / 2138
 1650 / 2138
 1660 / 2138
 1670 / 2138
 1680 / 2138
 1690 / 2138
 1700 / 2138
 1710 / 2138
 1720 / 2138
 1730 / 2138
 1740 / 2138
 1750 / 2138
 1760 / 2138
 1770 / 2138
 1780 / 2138
 1790 / 2138
 1800 / 2138
 1810 / 2138
 1820 / 2138
 1830 / 2138
 1840 / 2138
 1850 / 2138
 1860 / 2138
 1870 / 2138
 1880 / 2138
 1890 / 2138
 1900 / 2138
 1910 / 2138
 1920 / 2138
 1930 / 2138
 1940 / 2138
 1950 / 2138
 1960 / 2138
 1970 / 2138
 1980 / 2138
 1990 / 2138
 2000 / 2138
 2010 / 2138
 2020 / 2138
 2030 / 2138
 2040 / 2138
 2050 / 2138
 2060 / 2138
 2070 / 2138
 2080 / 2138
 2090 / 2138
 2100 / 2138
 2110 / 2138
 2120 / 2138
 2130 / 2138
## done.
## 
## Progressive alignment 2/2... 
## 
STEP     1 / 2137  f
STEP     2 / 2137  f
STEP     3 / 2137  f
STEP     4 / 2137  f
STEP     5 / 2137  f
STEP     6 / 2137  f
STEP     7 / 2137  f
STEP     8 / 2137  f
STEP     9 / 2137  f
STEP    10 / 2137  f
STEP    11 / 2137  f
STEP    12 / 2137  f
STEP    13 / 2137  f
STEP    14 / 2137  f
STEP    15 / 2137  f
STEP    16 / 2137  f
STEP    17 / 2137  f
STEP    18 / 2137  f
STEP    19 / 2137  f
STEP    20 / 2137  f
STEP    21 / 2137  f
STEP    22 / 2137  f
STEP    23 / 2137  f
STEP    24 / 2137  f
STEP    25 / 2137  f
STEP    26 / 2137  f
STEP    27 / 2137  f
STEP    28 / 2137  f
STEP    29 / 2137  f
STEP    30 / 2137  f
STEP    31 / 2137  f
STEP    32 / 2137  f
STEP    33 / 2137  f
STEP    34 / 2137  f
STEP    35 / 2137  f
STEP    36 / 2137  f
STEP    37 / 2137  f
STEP    38 / 2137  f
STEP    39 / 2137  f
STEP    40 / 2137  f
STEP    41 / 2137  f
STEP    42 / 2137  f
STEP    43 / 2137  f
STEP    44 / 2137  f
STEP    45 / 2137  f
STEP    46 / 2137  f
STEP    47 / 2137  f
STEP    48 / 2137  f
STEP    49 / 2137  f
STEP    50 / 2137  f
STEP    51 / 2137  f
STEP    52 / 2137  f
STEP    53 / 2137  f
STEP    54 / 2137  f
STEP    55 / 2137  f
STEP    56 / 2137  f
STEP    57 / 2137  f
STEP    58 / 2137  f
STEP    59 / 2137  f
STEP    60 / 2137  f
STEP    61 / 2137  f
STEP    62 / 2137  f
STEP    63 / 2137  f
STEP    64 / 2137  f
STEP    65 / 2137  f
STEP    66 / 2137  f
STEP    67 / 2137  f
STEP    68 / 2137  f
STEP    69 / 2137  f
STEP    70 / 2137  f
STEP    71 / 2137  f
STEP    72 / 2137  f
STEP    73 / 2137  f
STEP    74 / 2137  f
STEP    75 / 2137  f
STEP    76 / 2137  f
STEP    77 / 2137  f
STEP    78 / 2137  f
STEP    79 / 2137  f
STEP    80 / 2137  f
STEP    81 / 2137  f
STEP    82 / 2137  f
STEP    83 / 2137  f
STEP    84 / 2137  f
STEP    85 / 2137  f
STEP    86 / 2137  f
STEP    87 / 2137  f
STEP    88 / 2137  f
STEP    89 / 2137  f
STEP    90 / 2137  f
STEP    91 / 2137  f
STEP    92 / 2137  f
STEP    93 / 2137  f
STEP    94 / 2137  f
STEP    95 / 2137  f
STEP    96 / 2137  f
STEP    97 / 2137  f
STEP    98 / 2137  f
STEP    99 / 2137  f
STEP   100 / 2137  f
STEP   101 / 2137  f
STEP   102 / 2137  f
STEP   103 / 2137  f
STEP   104 / 2137  f
STEP   105 / 2137  f
STEP   106 / 2137  f
STEP   107 / 2137  f
STEP   108 / 2137  f
STEP   109 / 2137  f
STEP   110 / 2137  f
STEP   111 / 2137  f
STEP   112 / 2137  f
STEP   113 / 2137  f
STEP   114 / 2137  f
STEP   115 / 2137  f
STEP   116 / 2137  f
STEP   117 / 2137  f
STEP   118 / 2137  f
STEP   119 / 2137  f
STEP   120 / 2137  f
STEP   121 / 2137  f
STEP   122 / 2137  f
STEP   123 / 2137  f
STEP   124 / 2137  f
STEP   125 / 2137  f
STEP   126 / 2137  f
STEP   127 / 2137  f
STEP   128 / 2137  f
STEP   129 / 2137  f
STEP   130 / 2137  f
STEP   131 / 2137  f
STEP   132 / 2137  f
STEP   133 / 2137  f
STEP   134 / 2137  f
STEP   135 / 2137  f
STEP   136 / 2137  f
STEP   137 / 2137  f
STEP   138 / 2137  f
STEP   139 / 2137  f
STEP   140 / 2137  f
STEP   141 / 2137  f
STEP   142 / 2137  f
STEP   143 / 2137  f
STEP   144 / 2137  f
STEP   145 / 2137  f
STEP   146 / 2137  f
STEP   147 / 2137  f
STEP   148 / 2137  f
STEP   149 / 2137  f
STEP   150 / 2137  f
STEP   151 / 2137  f
STEP   152 / 2137  f
STEP   153 / 2137  f
STEP   154 / 2137  f
STEP   155 / 2137  f
STEP   156 / 2137  f
STEP   157 / 2137  f
STEP   158 / 2137  f
STEP   159 / 2137  f
STEP   160 / 2137  f
STEP   161 / 2137  f
STEP   162 / 2137  f
STEP   163 / 2137  f
STEP   164 / 2137  f
STEP   165 / 2137  f
STEP   166 / 2137  f
STEP   167 / 2137  f
STEP   168 / 2137  f
STEP   169 / 2137  f
STEP   170 / 2137  f
STEP   171 / 2137  f
STEP   172 / 2137  f
STEP   173 / 2137  f
STEP   174 / 2137  f
STEP   175 / 2137  f
STEP   176 / 2137  f
STEP   177 / 2137  f
STEP   178 / 2137  f
STEP   179 / 2137  f
STEP   180 / 2137  f
STEP   181 / 2137  f
STEP   182 / 2137  f
STEP   183 / 2137  f
STEP   184 / 2137  f
STEP   185 / 2137  f
STEP   186 / 2137  f
STEP   187 / 2137  f
STEP   188 / 2137  f
STEP   189 / 2137  f
STEP   190 / 2137  f
STEP   191 / 2137  f
STEP   192 / 2137  f
STEP   193 / 2137  f
STEP   194 / 2137  f
STEP   195 / 2137  f
STEP   196 / 2137  f
STEP   197 / 2137  f
STEP   198 / 2137  f
STEP   199 / 2137  f
STEP   200 / 2137  f
STEP   201 / 2137  f
STEP   202 / 2137  f
STEP   203 / 2137  f
STEP   204 / 2137  f
STEP   205 / 2137  f
STEP   206 / 2137  f
STEP   207 / 2137  f
STEP   208 / 2137  f
STEP   209 / 2137  f
STEP   210 / 2137  f
STEP   211 / 2137  f
STEP   212 / 2137  f
STEP   213 / 2137  f
STEP   214 / 2137  f
STEP   215 / 2137  f
STEP   216 / 2137  f
STEP   217 / 2137  f
STEP   218 / 2137  f
STEP   219 / 2137  f
STEP   220 / 2137  f
STEP   221 / 2137  f
STEP   222 / 2137  f
STEP   223 / 2137  f
STEP   224 / 2137  f
STEP   225 / 2137  f
STEP   226 / 2137  f
STEP   227 / 2137  f
STEP   228 / 2137  f
STEP   229 / 2137  f
STEP   230 / 2137  f
STEP   231 / 2137  f
STEP   232 / 2137  f
STEP   233 / 2137  f
STEP   234 / 2137  f
STEP   235 / 2137  f
STEP   236 / 2137  f
STEP   237 / 2137  f
STEP   238 / 2137  f
STEP   239 / 2137  f
STEP   240 / 2137  f
STEP   241 / 2137  f
STEP   242 / 2137  f
STEP   243 / 2137  f
STEP   244 / 2137  f
STEP   245 / 2137  f
STEP   246 / 2137  f
STEP   247 / 2137  f
STEP   248 / 2137  f
STEP   249 / 2137  f
STEP   250 / 2137  f
STEP   251 / 2137  f
STEP   252 / 2137  f
STEP   253 / 2137  f
STEP   254 / 2137  f
STEP   255 / 2137  f
STEP   256 / 2137  f
STEP   257 / 2137  f
STEP   258 / 2137  f
STEP   259 / 2137  f
STEP   260 / 2137  f
STEP   261 / 2137  f
STEP   262 / 2137  f
STEP   263 / 2137  f
STEP   264 / 2137  f
STEP   265 / 2137  f
STEP   266 / 2137  f
STEP   267 / 2137  f
STEP   268 / 2137  f
STEP   269 / 2137  f
STEP   270 / 2137  f
STEP   271 / 2137  f
STEP   272 / 2137  f
STEP   273 / 2137  f
STEP   274 / 2137  f
STEP   275 / 2137  f
STEP   276 / 2137  f
STEP   277 / 2137  f
STEP   278 / 2137  f
STEP   279 / 2137  f
STEP   280 / 2137  f
STEP   281 / 2137  f
STEP   282 / 2137  f
STEP   283 / 2137  f
STEP   284 / 2137  f
STEP   285 / 2137  f
STEP   286 / 2137  f
STEP   287 / 2137  f
STEP   288 / 2137  f
STEP   289 / 2137  f
STEP   290 / 2137  f
STEP   291 / 2137  f
STEP   292 / 2137  f
STEP   293 / 2137  f
STEP   294 / 2137  f
STEP   295 / 2137  f
STEP   296 / 2137  f
STEP   297 / 2137  f
STEP   298 / 2137  f
STEP   299 / 2137  f
STEP   300 / 2137  f
STEP   301 / 2137  f
STEP   302 / 2137  f
STEP   303 / 2137  f
STEP   304 / 2137  f
STEP   305 / 2137  f
STEP   306 / 2137  f
STEP   307 / 2137  f
STEP   308 / 2137  f
STEP   309 / 2137  f
STEP   310 / 2137  f
STEP   311 / 2137  f
STEP   312 / 2137  f
STEP   313 / 2137  f
STEP   314 / 2137  f
STEP   315 / 2137  f
STEP   316 / 2137  f
STEP   317 / 2137  f
STEP   318 / 2137  f
STEP   319 / 2137  f
STEP   320 / 2137  f
STEP   321 / 2137  f
STEP   322 / 2137  f
STEP   323 / 2137  f
STEP   324 / 2137  f
STEP   325 / 2137  f
STEP   326 / 2137  f
STEP   327 / 2137  f
STEP   328 / 2137  f
STEP   329 / 2137  f
STEP   330 / 2137  f
STEP   331 / 2137  f
STEP   332 / 2137  f
STEP   333 / 2137  f
STEP   334 / 2137  f
STEP   335 / 2137  f
STEP   336 / 2137  f
STEP   337 / 2137  f
STEP   338 / 2137  f
STEP   339 / 2137  f
STEP   340 / 2137  f
STEP   341 / 2137  f
STEP   342 / 2137  f
STEP   343 / 2137  f
STEP   344 / 2137  f
STEP   345 / 2137  f
STEP   346 / 2137  f
STEP   347 / 2137  f
STEP   348 / 2137  f
STEP   349 / 2137  f
STEP   350 / 2137  f
STEP   351 / 2137  f
STEP   352 / 2137  f
STEP   353 / 2137  f
STEP   354 / 2137  f
STEP   355 / 2137  f
STEP   356 / 2137  f
STEP   357 / 2137  f
STEP   358 / 2137  f
STEP   359 / 2137  f
STEP   360 / 2137  f
STEP   361 / 2137  f
STEP   362 / 2137  f
STEP   363 / 2137  f
STEP   364 / 2137  f
STEP   365 / 2137  f
STEP   366 / 2137  f
STEP   367 / 2137  f
STEP   368 / 2137  f
STEP   369 / 2137  f
STEP   370 / 2137  f
STEP   371 / 2137  f
STEP   372 / 2137  f
STEP   373 / 2137  f
STEP   374 / 2137  f
STEP   375 / 2137  f
STEP   376 / 2137  f
STEP   377 / 2137  f
STEP   378 / 2137  f
STEP   379 / 2137  f
STEP   380 / 2137  f
STEP   381 / 2137  f
STEP   382 / 2137  f
STEP   383 / 2137  f
STEP   384 / 2137  f
STEP   385 / 2137  f
STEP   386 / 2137  f
STEP   387 / 2137  f
STEP   388 / 2137  f
STEP   389 / 2137  f
STEP   390 / 2137  f
STEP   391 / 2137  f
STEP   392 / 2137  f
STEP   393 / 2137  f
STEP   394 / 2137  f
STEP   395 / 2137  f
STEP   396 / 2137  f
STEP   397 / 2137  f
STEP   398 / 2137  f
STEP   399 / 2137  f
STEP   400 / 2137  f
STEP   401 / 2137  f
STEP   402 / 2137  f
STEP   403 / 2137  f
STEP   404 / 2137  f
STEP   405 / 2137  f
STEP   406 / 2137  f
STEP   407 / 2137  f
STEP   408 / 2137  f
STEP   409 / 2137  f
STEP   410 / 2137  f
STEP   411 / 2137  f
STEP   412 / 2137  f
STEP   413 / 2137  f
STEP   414 / 2137  f
STEP   415 / 2137  f
STEP   416 / 2137  f
STEP   417 / 2137  f
STEP   418 / 2137  f
STEP   419 / 2137  f
STEP   420 / 2137  f
STEP   421 / 2137  f
STEP   422 / 2137  f
STEP   423 / 2137  f
STEP   424 / 2137  f
STEP   425 / 2137  f
STEP   426 / 2137  f
STEP   427 / 2137  f
STEP   428 / 2137  f
STEP   429 / 2137  f
STEP   430 / 2137  f
STEP   431 / 2137  f
STEP   432 / 2137  f
STEP   433 / 2137  f
STEP   434 / 2137  f
STEP   435 / 2137  f
STEP   436 / 2137  f
STEP   437 / 2137  f
STEP   438 / 2137  f
STEP   439 / 2137  f
STEP   440 / 2137  f
STEP   441 / 2137  f
STEP   442 / 2137  f
STEP   443 / 2137  f
STEP   444 / 2137  f
STEP   445 / 2137  f
STEP   446 / 2137  f
STEP   447 / 2137  f
STEP   448 / 2137  f
STEP   449 / 2137  f
STEP   450 / 2137  f
STEP   451 / 2137  f
STEP   452 / 2137  f
STEP   453 / 2137  f
STEP   454 / 2137  f
STEP   455 / 2137  f
STEP   456 / 2137  f
STEP   457 / 2137  f
STEP   458 / 2137  f
STEP   459 / 2137  f
STEP   460 / 2137  f
STEP   461 / 2137  f
STEP   462 / 2137  f
STEP   463 / 2137  f
STEP   464 / 2137  f
STEP   465 / 2137  f
STEP   466 / 2137  f
STEP   467 / 2137  f
STEP   468 / 2137  f
STEP   469 / 2137  f
STEP   470 / 2137  f
STEP   471 / 2137  f
STEP   472 / 2137  f
STEP   473 / 2137  f
STEP   474 / 2137  f
STEP   475 / 2137  f
STEP   476 / 2137  f
STEP   477 / 2137  f
STEP   478 / 2137  f
STEP   479 / 2137  f
STEP   480 / 2137  f
STEP   481 / 2137  f
STEP   482 / 2137  f
STEP   483 / 2137  f
STEP   484 / 2137  f
STEP   485 / 2137  f
STEP   486 / 2137  f
STEP   487 / 2137  f
STEP   488 / 2137  f
STEP   489 / 2137  f
STEP   490 / 2137  f
STEP   491 / 2137  f
STEP   492 / 2137  f
STEP   493 / 2137  f
STEP   494 / 2137  f
STEP   495 / 2137  f
STEP   496 / 2137  f
STEP   497 / 2137  f
STEP   498 / 2137  f
STEP   499 / 2137  f
STEP   500 / 2137  f
STEP   501 / 2137  f
STEP   601 / 2137  f
STEP   701 / 2137  f
STEP   801 / 2137  f
STEP   901 / 2137  f
STEP  1001 / 2137  f
STEP  1101 / 2137  f
STEP  1201 / 2137  f
STEP  1301 / 2137  f
STEP  1401 / 2137  f
STEP  1501 / 2137  f
STEP  1601 / 2137  f
STEP  1701 / 2137  f
STEP  1801 / 2137  f
STEP  1901 / 2137  f
## Reallocating..done. *alloclen = 1934
## 
STEP  2001 / 2137  f
STEP  2101 / 2137  f
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
##       0.10 seconds: Top hits for    735 of   2138 seqs (at seed    400)
##       0.20 seconds: Top hits for   1604 of   2138 seqs (at seed   1200)
##       0.33 seconds: Joined    100 of   2135
##       0.47 seconds: Joined    300 of   2135
##       0.66 seconds: Joined    500 of   2135
##       0.77 seconds: Joined    600 of   2135
##       0.88 seconds: Joined    800 of   2135
##       1.02 seconds: Joined   1000 of   2135
##       1.17 seconds: Joined   1300 of   2135
##       1.34 seconds: Joined   1500 of   2135
##       1.50 seconds: Joined   1800 of   2135
##       1.61 seconds: Joined   2000 of   2135
## Initial topology in 1.69 seconds
## Refining topology: 44 rounds ME-NNIs, 2 rounds ME-SPRs, 22 rounds ML-NNIs
##       1.71 seconds: ME NNI round 1 of 44, 601 of 2136 splits, 112 changes (max delta 0.021)
##       1.81 seconds: ME NNI round 2 of 44, 1001 of 2136 splits, 118 changes (max delta 0.020)
##       1.92 seconds: ME NNI round 3 of 44, 1401 of 2136 splits, 114 changes (max delta 0.023)
##       2.02 seconds: ME NNI round 6 of 44, 101 of 2136 splits, 12 changes (max delta 0.008)
##       2.14 seconds: SPR round   1 of   2, 101 of 4274 nodes
##       2.28 seconds: SPR round   1 of   2, 301 of 4274 nodes
##       2.44 seconds: SPR round   1 of   2, 501 of 4274 nodes
##       2.57 seconds: SPR round   1 of   2, 701 of 4274 nodes
##       2.71 seconds: SPR round   1 of   2, 901 of 4274 nodes
##       2.84 seconds: SPR round   1 of   2, 1101 of 4274 nodes
##       2.98 seconds: SPR round   1 of   2, 1301 of 4274 nodes
##       3.12 seconds: SPR round   1 of   2, 1501 of 4274 nodes
##       3.28 seconds: SPR round   1 of   2, 1701 of 4274 nodes
##       3.42 seconds: SPR round   1 of   2, 1901 of 4274 nodes
##       3.56 seconds: SPR round   1 of   2, 2101 of 4274 nodes
##       3.69 seconds: SPR round   1 of   2, 2301 of 4274 nodes
##       3.83 seconds: SPR round   1 of   2, 2501 of 4274 nodes
##       3.96 seconds: SPR round   1 of   2, 2701 of 4274 nodes
##       4.11 seconds: SPR round   1 of   2, 2901 of 4274 nodes
##       4.24 seconds: SPR round   1 of   2, 3101 of 4274 nodes
##       4.38 seconds: SPR round   1 of   2, 3301 of 4274 nodes
##       4.51 seconds: SPR round   1 of   2, 3501 of 4274 nodes
##       4.64 seconds: SPR round   1 of   2, 3701 of 4274 nodes
##       4.79 seconds: SPR round   1 of   2, 3901 of 4274 nodes
##       4.93 seconds: SPR round   1 of   2, 4101 of 4274 nodes
##       5.06 seconds: ME NNI round 15 of 44, 1 of 2136 splits
##       5.16 seconds: ME NNI round 16 of 44, 501 of 2136 splits, 2 changes (max delta 0.001)
##       5.32 seconds: SPR round   2 of   2, 101 of 4274 nodes
##       5.45 seconds: SPR round   2 of   2, 301 of 4274 nodes
##       5.58 seconds: SPR round   2 of   2, 501 of 4274 nodes
##       5.74 seconds: SPR round   2 of   2, 701 of 4274 nodes
##       5.88 seconds: SPR round   2 of   2, 901 of 4274 nodes
##       6.02 seconds: SPR round   2 of   2, 1101 of 4274 nodes
##       6.17 seconds: SPR round   2 of   2, 1301 of 4274 nodes
##       6.29 seconds: SPR round   2 of   2, 1501 of 4274 nodes
##       6.44 seconds: SPR round   2 of   2, 1701 of 4274 nodes
##       6.56 seconds: SPR round   2 of   2, 1901 of 4274 nodes
##       6.69 seconds: SPR round   2 of   2, 2101 of 4274 nodes
##       6.83 seconds: SPR round   2 of   2, 2301 of 4274 nodes
##       6.96 seconds: SPR round   2 of   2, 2501 of 4274 nodes
##       7.11 seconds: SPR round   2 of   2, 2701 of 4274 nodes
##       7.22 seconds: SPR round   2 of   2, 2901 of 4274 nodes
##       7.35 seconds: SPR round   2 of   2, 3101 of 4274 nodes
##       7.50 seconds: SPR round   2 of   2, 3301 of 4274 nodes
##       7.64 seconds: SPR round   2 of   2, 3501 of 4274 nodes
##       7.79 seconds: SPR round   2 of   2, 3701 of 4274 nodes
##       7.94 seconds: SPR round   2 of   2, 3901 of 4274 nodes
##       8.07 seconds: SPR round   2 of   2, 4101 of 4274 nodes
##       8.19 seconds: ME NNI round 29 of 44, 1 of 2136 splits
##       8.30 seconds: ME NNI round 30 of 44, 601 of 2136 splits, 0 changes
## Total branch-length 49.127 after 8.46 sec
##       8.47 seconds: ML Lengths 1 of 2136 splits
##       8.58 seconds: ML Lengths 301 of 2136 splits
##       8.70 seconds: ML Lengths 601 of 2136 splits
##       8.81 seconds: ML Lengths 901 of 2136 splits
##       8.92 seconds: ML Lengths 1201 of 2136 splits
##       9.03 seconds: ML Lengths 1501 of 2136 splits
##       9.15 seconds: ML Lengths 1801 of 2136 splits
##       9.26 seconds: ML Lengths 2101 of 2136 splits
##       9.39 seconds: ML NNI round 1 of 22, 101 of 2136 splits, 15 changes (max delta 11.697)
##       9.51 seconds: ML NNI round 1 of 22, 201 of 2136 splits, 36 changes (max delta 11.697)
##       9.64 seconds: ML NNI round 1 of 22, 301 of 2136 splits, 46 changes (max delta 11.697)
##       9.77 seconds: ML NNI round 1 of 22, 401 of 2136 splits, 63 changes (max delta 11.697)
##       9.91 seconds: ML NNI round 1 of 22, 501 of 2136 splits, 84 changes (max delta 11.697)
##      10.04 seconds: ML NNI round 1 of 22, 601 of 2136 splits, 102 changes (max delta 11.697)
##      10.16 seconds: ML NNI round 1 of 22, 701 of 2136 splits, 113 changes (max delta 11.697)
##      10.28 seconds: ML NNI round 1 of 22, 801 of 2136 splits, 126 changes (max delta 11.697)
##      10.40 seconds: ML NNI round 1 of 22, 901 of 2136 splits, 143 changes (max delta 12.535)
##      10.52 seconds: ML NNI round 1 of 22, 1001 of 2136 splits, 163 changes (max delta 12.535)
##      10.65 seconds: ML NNI round 1 of 22, 1101 of 2136 splits, 177 changes (max delta 12.535)
##      10.77 seconds: ML NNI round 1 of 22, 1201 of 2136 splits, 188 changes (max delta 12.535)
##      10.91 seconds: ML NNI round 1 of 22, 1301 of 2136 splits, 204 changes (max delta 12.535)
##      11.02 seconds: ML NNI round 1 of 22, 1401 of 2136 splits, 222 changes (max delta 12.535)
##      11.15 seconds: ML NNI round 1 of 22, 1501 of 2136 splits, 236 changes (max delta 12.535)
##      11.28 seconds: ML NNI round 1 of 22, 1601 of 2136 splits, 249 changes (max delta 15.017)
##      11.40 seconds: ML NNI round 1 of 22, 1701 of 2136 splits, 258 changes (max delta 15.017)
##      11.52 seconds: ML NNI round 1 of 22, 1801 of 2136 splits, 271 changes (max delta 15.017)
##      11.65 seconds: ML NNI round 1 of 22, 1901 of 2136 splits, 282 changes (max delta 15.017)
##      11.78 seconds: ML NNI round 1 of 22, 2001 of 2136 splits, 297 changes (max delta 15.017)
##      11.91 seconds: ML NNI round 1 of 22, 2101 of 2136 splits, 313 changes (max delta 15.017)
## ML-NNI round 1: LogLk = -142287.487 NNIs 317 max delta 15.02 Time 11.98
##      12.74 seconds: Optimizing GTR model, step 2 of 12
##      13.68 seconds: Optimizing GTR model, step 3 of 12
##      14.21 seconds: Optimizing GTR model, step 4 of 12
##      15.16 seconds: Optimizing GTR model, step 5 of 12
##      16.11 seconds: Optimizing GTR model, step 6 of 12
##      17.08 seconds: Optimizing GTR model, step 7 of 12
##      17.92 seconds: Optimizing GTR model, step 8 of 12
##      18.47 seconds: Optimizing GTR model, step 9 of 12
##      19.18 seconds: Optimizing GTR model, step 10 of 12
##      19.76 seconds: Optimizing GTR model, step 11 of 12
##      20.50 seconds: Optimizing GTR model, step 12 of 12
## GTR Frequencies: 0.2119 0.3009 0.2493 0.2379
## GTR rates(ac ag at cg ct gt) 1.3353 3.2040 1.6103 0.8988 2.3946 1.0000
##      21.16 seconds: ML Lengths 1 of 2136 splits
##      21.26 seconds: ML Lengths 201 of 2136 splits
##      21.41 seconds: ML Lengths 501 of 2136 splits
##      21.55 seconds: ML Lengths 801 of 2136 splits
##      21.70 seconds: ML Lengths 1101 of 2136 splits
##      21.85 seconds: ML Lengths 1401 of 2136 splits
##      22.00 seconds: ML Lengths 1701 of 2136 splits
##      22.15 seconds: ML Lengths 2001 of 2136 splits
##      22.29 seconds: Site likelihoods with rate category 1 of 20
##      22.42 seconds: Site likelihoods with rate category 3 of 20
##      22.56 seconds: Site likelihoods with rate category 5 of 20
##      22.69 seconds: Site likelihoods with rate category 7 of 20
##      22.84 seconds: Site likelihoods with rate category 9 of 20
##      22.97 seconds: Site likelihoods with rate category 11 of 20
##      23.11 seconds: Site likelihoods with rate category 13 of 20
##      23.24 seconds: Site likelihoods with rate category 15 of 20
##      23.38 seconds: Site likelihoods with rate category 17 of 20
##      23.51 seconds: Site likelihoods with rate category 19 of 20
## Switched to using 20 rate categories (CAT approximation)
## Rate categories were divided by 1.175 so that average rate = 1.0
## CAT-based log-likelihoods may not be comparable across runs
## Use -gamma for approximate but comparable Gamma(20) log-likelihoods
##      23.70 seconds: ML NNI round 2 of 22, 1 of 2136 splits
##      23.84 seconds: ML NNI round 2 of 22, 101 of 2136 splits, 7 changes (max delta 4.222)
##      23.99 seconds: ML NNI round 2 of 22, 201 of 2136 splits, 21 changes (max delta 5.966)
##      24.16 seconds: ML NNI round 2 of 22, 301 of 2136 splits, 33 changes (max delta 5.966)
##      24.31 seconds: ML NNI round 2 of 22, 401 of 2136 splits, 40 changes (max delta 5.966)
##      24.49 seconds: ML NNI round 2 of 22, 501 of 2136 splits, 56 changes (max delta 5.966)
##      24.64 seconds: ML NNI round 2 of 22, 601 of 2136 splits, 70 changes (max delta 5.966)
##      24.79 seconds: ML NNI round 2 of 22, 701 of 2136 splits, 83 changes (max delta 6.443)
##      24.95 seconds: ML NNI round 2 of 22, 801 of 2136 splits, 98 changes (max delta 6.443)
##      25.10 seconds: ML NNI round 2 of 22, 901 of 2136 splits, 113 changes (max delta 6.443)
##      25.25 seconds: ML NNI round 2 of 22, 1001 of 2136 splits, 120 changes (max delta 6.443)
##      25.41 seconds: ML NNI round 2 of 22, 1101 of 2136 splits, 134 changes (max delta 6.443)
##      25.56 seconds: ML NNI round 2 of 22, 1201 of 2136 splits, 146 changes (max delta 6.443)
##      25.70 seconds: ML NNI round 2 of 22, 1301 of 2136 splits, 153 changes (max delta 8.316)
##      25.84 seconds: ML NNI round 2 of 22, 1401 of 2136 splits, 164 changes (max delta 8.316)
##      26.01 seconds: ML NNI round 2 of 22, 1501 of 2136 splits, 174 changes (max delta 8.316)
##      26.17 seconds: ML NNI round 2 of 22, 1601 of 2136 splits, 189 changes (max delta 8.316)
##      26.30 seconds: ML NNI round 2 of 22, 1701 of 2136 splits, 197 changes (max delta 8.316)
##      26.44 seconds: ML NNI round 2 of 22, 1801 of 2136 splits, 206 changes (max delta 8.316)
##      26.59 seconds: ML NNI round 2 of 22, 1901 of 2136 splits, 222 changes (max delta 8.316)
##      26.75 seconds: ML NNI round 2 of 22, 2001 of 2136 splits, 236 changes (max delta 8.316)
##      26.89 seconds: ML NNI round 2 of 22, 2101 of 2136 splits, 246 changes (max delta 8.316)
## ML-NNI round 2: LogLk = -117490.895 NNIs 250 max delta 8.32 Time 26.95
##      27.12 seconds: ML NNI round 3 of 22, 101 of 2136 splits, 10 changes (max delta 1.001)
##      27.26 seconds: ML NNI round 3 of 22, 201 of 2136 splits, 17 changes (max delta 2.419)
##      27.41 seconds: ML NNI round 3 of 22, 301 of 2136 splits, 25 changes (max delta 3.247)
##      27.53 seconds: ML NNI round 3 of 22, 401 of 2136 splits, 28 changes (max delta 3.247)
##      27.68 seconds: ML NNI round 3 of 22, 501 of 2136 splits, 29 changes (max delta 3.247)
##      27.82 seconds: ML NNI round 3 of 22, 601 of 2136 splits, 33 changes (max delta 3.247)
##      27.98 seconds: ML NNI round 3 of 22, 701 of 2136 splits, 42 changes (max delta 3.247)
##      28.15 seconds: ML NNI round 3 of 22, 801 of 2136 splits, 46 changes (max delta 3.247)
##      28.31 seconds: ML NNI round 3 of 22, 901 of 2136 splits, 58 changes (max delta 3.247)
##      28.44 seconds: ML NNI round 3 of 22, 1001 of 2136 splits, 62 changes (max delta 3.247)
##      28.58 seconds: ML NNI round 3 of 22, 1101 of 2136 splits, 67 changes (max delta 3.247)
##      28.72 seconds: ML NNI round 3 of 22, 1201 of 2136 splits, 78 changes (max delta 3.247)
##      28.86 seconds: ML NNI round 3 of 22, 1301 of 2136 splits, 90 changes (max delta 5.076)
##      29.00 seconds: ML NNI round 3 of 22, 1401 of 2136 splits, 94 changes (max delta 5.076)
## ML-NNI round 3: LogLk = -117415.173 NNIs 96 max delta 5.08 Time 29.08
##      29.23 seconds: ML NNI round 4 of 22, 101 of 2136 splits, 6 changes (max delta 0.778)
##      29.38 seconds: ML NNI round 4 of 22, 201 of 2136 splits, 11 changes (max delta 7.644)
##      29.54 seconds: ML NNI round 4 of 22, 301 of 2136 splits, 27 changes (max delta 10.662)
##      29.68 seconds: ML NNI round 4 of 22, 401 of 2136 splits, 30 changes (max delta 10.662)
##      29.81 seconds: ML NNI round 4 of 22, 501 of 2136 splits, 33 changes (max delta 10.662)
##      29.95 seconds: ML NNI round 4 of 22, 601 of 2136 splits, 39 changes (max delta 10.662)
##      30.10 seconds: ML NNI round 4 of 22, 701 of 2136 splits, 48 changes (max delta 10.662)
##      30.22 seconds: ML NNI round 4 of 22, 801 of 2136 splits, 50 changes (max delta 10.662)
##      30.35 seconds: ML NNI round 4 of 22, 901 of 2136 splits, 54 changes (max delta 10.662)
## ML-NNI round 4: LogLk = -117358.101 NNIs 56 max delta 10.66 Time 30.49
##      30.48 seconds: ML NNI round 5 of 22, 1 of 2136 splits
##      30.61 seconds: ML NNI round 5 of 22, 101 of 2136 splits, 3 changes (max delta 0.185)
##      30.76 seconds: ML NNI round 5 of 22, 201 of 2136 splits, 12 changes (max delta 7.978)
##      30.93 seconds: ML NNI round 5 of 22, 301 of 2136 splits, 16 changes (max delta 7.978)
##      31.07 seconds: ML NNI round 5 of 22, 401 of 2136 splits, 23 changes (max delta 7.978)
## ML-NNI round 5: LogLk = -117336.432 NNIs 25 max delta 7.98 Time 31.19
##      31.18 seconds: ML NNI round 6 of 22, 1 of 2136 splits
##      31.30 seconds: ML NNI round 6 of 22, 101 of 2136 splits, 2 changes (max delta 0.020)
##      31.45 seconds: ML NNI round 6 of 22, 201 of 2136 splits, 7 changes (max delta 2.015)
##      31.59 seconds: ML NNI round 6 of 22, 301 of 2136 splits, 14 changes (max delta 2.015)
## ML-NNI round 6: LogLk = -117331.835 NNIs 14 max delta 2.01 Time 31.64
##      31.77 seconds: ML NNI round 7 of 22, 101 of 2136 splits, 5 changes (max delta 0.728)
##      31.90 seconds: ML NNI round 7 of 22, 201 of 2136 splits, 11 changes (max delta 0.728)
## ML-NNI round 7: LogLk = -117330.407 NNIs 11 max delta 0.73 Time 31.93
##      32.08 seconds: ML NNI round 8 of 22, 101 of 2136 splits, 7 changes (max delta 0.228)
## ML-NNI round 8: LogLk = -117329.612 NNIs 8 max delta 0.65 Time 32.15
## ML-NNI round 9: LogLk = -117330.205 NNIs 3 max delta 0.21 Time 32.28
## Turning off heuristics for final round of ML NNIs (converged)
##      32.28 seconds: ML NNI round 10 of 22, 1 of 2136 splits
##      32.48 seconds: ML NNI round 10 of 22, 101 of 2136 splits, 2 changes (max delta 0.124)
##      32.69 seconds: ML NNI round 10 of 22, 201 of 2136 splits, 5 changes (max delta 0.124)
##      32.88 seconds: ML NNI round 10 of 22, 301 of 2136 splits, 11 changes (max delta 1.640)
##      33.07 seconds: ML NNI round 10 of 22, 401 of 2136 splits, 11 changes (max delta 1.640)
##      33.26 seconds: ML NNI round 10 of 22, 501 of 2136 splits, 13 changes (max delta 6.651)
##      33.44 seconds: ML NNI round 10 of 22, 601 of 2136 splits, 14 changes (max delta 6.651)
##      33.62 seconds: ML NNI round 10 of 22, 701 of 2136 splits, 18 changes (max delta 6.651)
##      33.81 seconds: ML NNI round 10 of 22, 801 of 2136 splits, 19 changes (max delta 6.651)
##      34.01 seconds: ML NNI round 10 of 22, 901 of 2136 splits, 21 changes (max delta 6.651)
##      34.21 seconds: ML NNI round 10 of 22, 1001 of 2136 splits, 26 changes (max delta 6.651)
##      34.43 seconds: ML NNI round 10 of 22, 1101 of 2136 splits, 29 changes (max delta 6.651)
##      34.61 seconds: ML NNI round 10 of 22, 1201 of 2136 splits, 29 changes (max delta 6.651)
##      34.80 seconds: ML NNI round 10 of 22, 1301 of 2136 splits, 30 changes (max delta 6.651)
##      34.99 seconds: ML NNI round 10 of 22, 1401 of 2136 splits, 31 changes (max delta 6.651)
##      35.17 seconds: ML NNI round 10 of 22, 1501 of 2136 splits, 36 changes (max delta 6.651)
##      35.36 seconds: ML NNI round 10 of 22, 1601 of 2136 splits, 41 changes (max delta 6.651)
##      35.57 seconds: ML NNI round 10 of 22, 1701 of 2136 splits, 41 changes (max delta 6.651)
##      35.76 seconds: ML NNI round 10 of 22, 1801 of 2136 splits, 43 changes (max delta 6.651)
##      35.94 seconds: ML NNI round 10 of 22, 1901 of 2136 splits, 45 changes (max delta 6.651)
##      36.14 seconds: ML NNI round 10 of 22, 2001 of 2136 splits, 48 changes (max delta 6.651)
##      36.33 seconds: ML NNI round 10 of 22, 2101 of 2136 splits, 49 changes (max delta 6.651)
## ML-NNI round 10: LogLk = -117189.113 NNIs 50 max delta 6.65 Time 36.42 (final)
##      36.47 seconds: ML Lengths 101 of 2136 splits
##      36.58 seconds: ML Lengths 301 of 2136 splits
##      36.69 seconds: ML Lengths 501 of 2136 splits
##      36.79 seconds: ML Lengths 701 of 2136 splits
##      36.90 seconds: ML Lengths 901 of 2136 splits
##      37.00 seconds: ML Lengths 1101 of 2136 splits
##      37.11 seconds: ML Lengths 1301 of 2136 splits
##      37.21 seconds: ML Lengths 1501 of 2136 splits
##      37.32 seconds: ML Lengths 1701 of 2136 splits
##      37.42 seconds: ML Lengths 1901 of 2136 splits
##      37.53 seconds: ML Lengths 2101 of 2136 splits
## Optimize all lengths: LogLk = -117181.766 Time 37.56
##      37.79 seconds: ML split tests for    100 of   2135 internal splits
##      38.03 seconds: ML split tests for    200 of   2135 internal splits
##      38.27 seconds: ML split tests for    300 of   2135 internal splits
##      38.50 seconds: ML split tests for    400 of   2135 internal splits
##      38.72 seconds: ML split tests for    500 of   2135 internal splits
##      38.95 seconds: ML split tests for    600 of   2135 internal splits
##      39.18 seconds: ML split tests for    700 of   2135 internal splits
##      39.41 seconds: ML split tests for    800 of   2135 internal splits
##      39.64 seconds: ML split tests for    900 of   2135 internal splits
##      39.87 seconds: ML split tests for   1000 of   2135 internal splits
##      40.10 seconds: ML split tests for   1100 of   2135 internal splits
##      40.33 seconds: ML split tests for   1200 of   2135 internal splits
##      40.56 seconds: ML split tests for   1300 of   2135 internal splits
##      40.79 seconds: ML split tests for   1400 of   2135 internal splits
##      41.02 seconds: ML split tests for   1500 of   2135 internal splits
##      41.25 seconds: ML split tests for   1600 of   2135 internal splits
##      41.49 seconds: ML split tests for   1700 of   2135 internal splits
##      41.73 seconds: ML split tests for   1800 of   2135 internal splits
##      41.96 seconds: ML split tests for   1900 of   2135 internal splits
##      42.20 seconds: ML split tests for   2000 of   2135 internal splits
##      42.44 seconds: ML split tests for   2100 of   2135 internal splits
## Total time: 42.53 seconds Unique: 2138/2138 Bad splits: 11/2135 Worst delta-LogLk 2.181
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
##  date     2024-04-21
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
