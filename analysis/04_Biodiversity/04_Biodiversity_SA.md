---
title: "04_PreProcessing"
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



# Setting the Environment 

## Set the seed 

```r
set.seed(09091999)
```

## Load Libraries 

```r
pacman::p_load(tidyverse, devtools, patchwork, iNEXT, phyloseq,
               install = FALSE)
```

## Load in Data 

```r
load("/local/workdir/sna49/moon_milk/moonmilk/data/02_PreProcessing/raw_preprocessed_physeq.RData")
raw_preprocessed_physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4531 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 4531 taxa by 9 taxonomic ranks ]
```

```r
# Intiution Check 
min(sample_sums(raw_preprocessed_physeq))
```

```
## [1] 13720
```

```r
# Make metadata dataframe
metadata_df <-
  raw_preprocessed_physeq %>%
  sample_data() %>%
  data.frame()

#View(metadata_df)

# Setting colors for caves
cave_colors <- c(
  "Fata Apei Cave " = "dodgerblue4",
  "Ferice Cave " = "#D8BFD8",
  "Nestor Cave " = "#D9CC3C",
  "Tausoare Cave " = "#A0E0BA")
```


# Goals

1. Calculate the Hill Diversity of the samples. 
2. Evaluate the rarefaction curves. 
3. Evaluate the Diversity values. 
4. Makes notes of specific samples and their seq depth. 

# Diversity Calculations with iNEXT 


```r
# prepare input data 
iNEXT_input_df <- 
  raw_preprocessed_physeq %>%
  otu_table() %>%
  data.frame()

# Quick check
dim(iNEXT_input_df)
```

```
## [1] 4531   10
```

```r
# Run iNEXT: Calculate the Hill Numbers 
# Note that: Species in ROWS, Samples in COLUMNS 
# Remember to set the seed! 

iNEXT_data <- iNEXT(iNEXT_input_df, 
                    q = c(0,1,2), datatype = "abundance")

# Save the file
save(iNEXT_data, file = "/local/workdir/sna49/moon_milk/moonmilk/data/04_Biodiversity/iNEXT_data.RData")
```

# Evaluate the Diversity! 

```r
load("/local/workdir/sna49/moon_milk/moonmilk/data/04_Biodiversity/iNEXT_data.RData")
str(iNEXT_data)
```

```
## List of 3
##  $ DataInfo:'data.frame':	10 obs. of  14 variables:
##   ..$ Assemblage: chr [1:10] "ERR11588428" "ERR11588429" "ERR11588430" "ERR11588431" ...
##   ..$ n         : num [1:10] 21236 20768 20697 29177 28614 ...
##   ..$ S.obs     : num [1:10] 871 914 896 244 321 192 451 431 597 408
##   ..$ SC        : num [1:10] 1 0.999 1 1 1 ...
##   ..$ f1        : num [1:10] 8 17 9 0 0 0 2 0 5 4
##   ..$ f2        : num [1:10] 56 65 63 0 0 0 2 4 26 2
##   ..$ f3        : num [1:10] 22 21 19 0 1 0 2 2 14 5
##   ..$ f4        : num [1:10] 24 24 18 0 4 0 5 3 10 8
##   ..$ f5        : num [1:10] 31 42 32 3 6 0 4 9 10 14
##   ..$ f6        : num [1:10] 25 41 28 3 7 0 7 3 24 15
##   ..$ f7        : num [1:10] 26 27 25 7 3 1 5 8 26 14
##   ..$ f8        : num [1:10] 20 39 30 5 11 0 10 12 38 23
##   ..$ f9        : num [1:10] 28 24 29 5 5 2 7 10 19 17
##   ..$ f10       : num [1:10] 34 21 32 3 5 1 11 8 30 22
##  $ iNextEst:List of 2
##   ..$ size_based    :'data.frame':	1200 obs. of  10 variables:
##   .. ..$ Assemblage: chr [1:1200] "ERR11588428" "ERR11588428" "ERR11588428" "ERR11588428" ...
##   .. ..$ m         : num [1:1200] 1 1180 2360 3540 4719 ...
##   .. ..$ Method    : chr [1:1200] "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##   .. ..$ Order.q   : num [1:1200] 0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ qD        : num [1:1200] 1 533 677 742 779 ...
##   .. ..$ qD.LCL    : num [1:1200] 1 530 673 737 774 ...
##   .. ..$ qD.UCL    : num [1:1200] 1 536 681 747 784 ...
##   .. ..$ SC        : num [1:1200] 0.00192 0.80536 0.92403 0.96006 0.97554 ...
##   .. ..$ SC.LCL    : num [1:1200] 0.0019 0.8034 0.9227 0.9591 0.9748 ...
##   .. ..$ SC.UCL    : num [1:1200] 0.00195 0.80735 0.92536 0.96103 0.97626 ...
##   ..$ coverage_based:'data.frame':	948 obs. of  8 variables:
##   .. ..$ Assemblage: chr [1:948] "ERR11588428" "ERR11588428" "ERR11588428" "ERR11588428" ...
##   .. ..$ SC        : num [1:948] 0.00192 0.80536 0.92403 0.96006 0.97554 ...
##   .. ..$ m         : num [1:948] 1 1180 2360 3540 4719 ...
##   .. ..$ Method    : chr [1:948] "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##   .. ..$ Order.q   : num [1:948] 0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ qD        : num [1:948] 1 533 677 742 779 ...
##   .. ..$ qD.LCL    : num [1:948] 1 529 672 736 773 ...
##   .. ..$ qD.UCL    : num [1:948] 1 538 682 748 785 ...
##  $ AsyEst  :'data.frame':	30 obs. of  7 variables:
##   ..$ Assemblage: chr [1:30] "ERR11588428" "ERR11588428" "ERR11588428" "ERR11588429" ...
##   ..$ Diversity : chr [1:30] "Species richness" "Shannon diversity" "Simpson diversity" "Species richness" ...
##   ..$ Observed  : num [1:30] 871 626 508 914 634 ...
##   ..$ Estimator : num [1:30] 872 639 520 916 648 ...
##   ..$ s.e.      : num [1:30] 6.3 3.22 4.35 6.91 3.12 ...
##   ..$ LCL       : num [1:30] 871 633 511 914 642 ...
##   ..$ UCL       : num [1:30] 884 646 528 930 654 ...
##  - attr(*, "class")= chr "iNEXT"
```

```r
typeof(iNEXT_data)
```

```
## [1] "list"
```

# Plot Diversity 

```r
# Prepare Colors 
color_df <- 
  iNEXT_input_df %>%
  colnames() %>%
  data.frame()
# Check
head(color_df)
```

```
##             .
## 1 ERR11588428
## 2 ERR11588429
## 3 ERR11588430
## 4 ERR11588431
## 5 ERR11588432
## 6 ERR11588433
```

```r
# Rename the column 
colnames(color_df)[1] <- "names"
# Check
head(color_df)
```

```
##         names
## 1 ERR11588428
## 2 ERR11588429
## 3 ERR11588430
## 4 ERR11588431
## 5 ERR11588432
## 6 ERR11588433
```

```r
# Make a helper dataframe for plotting with colors 
iNEXT_color_df <- 
  color_df %>%
  # Merge with metadata
  left_join(metadata_df, by = "names") %>%
  # Merge with colors for plotting with ggiNEXT
  left_join(data.frame(cave_colors = cave_colors,
            cave = names(cave_colors)),
            by = "cave")
```

# Plot Rarefaction with `ggiNEXT`


```r
# Plot rarefaction! 
# rarefaction/extrapolation curve, type = 1 

# Order q: 
  # 0 = Richness/ Number of Total taxa
  # 1 = Exponential Shannon / Number of "Common" taxa
  # 2 = Inverse Simpson / Number of "Dominant" taxa 

ggiNEXT(iNEXT_data, type = 1, facet.var = "Order.q") + 
  facet_wrap(~Order.q, scales = "fixed") + 
  scale_color_manual(values = iNEXT_color_df$cave_colors, guide = FALSE) + 
  scale_fill_manual(values = iNEXT_color_df$cave_colors, guide = FALSE) + 
  scale_shape_manual(values = base::rep(17, nsamples(raw_preprocessed_physeq)),
                     guide = FALSE) +
  labs(x = "Number of Sequences (Library Size)", 
       y = "Effective Number of ASVs") + 
  theme_bw() 
```

```
## Scale for colour is already present.
## Adding another scale for colour, which will replace the existing scale.
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

<img src="../figures/04_Biodiversity/ggiNEXT-1.png" style="display: block; margin: auto;" />

Remembering that an Order q of:  

  - 0 = Richness/ Number of Total taxa
  - 1 = Exponential Shannon / Number of "Common" taxa
  - 2 = Inverse Simpson / Number of "Dominant" taxa 
  
*Also note that iNEXT will extrapolate the number of sequences to double the sequencing depth of the sample!*

Finally, let's draw some conclusions from the plot above (we discussed in class & I'll summarize here):  

1. The best gift of a diversity measure: *All three diversity measures are on the same x-axis: The Effective Number of Species.*  
2. As q increases from 0 to 1 to 2, the abundance of the species is weighted more, so the rarefaction curve platueaus earlier. **This is important because it indicates that richess is the most sensitive to sequencing depth!** However, Shannon and even less sensitive is Simpson!  


# Manually plot Diversity 

## Rarefaction

```r
iNEXT_manual_df <- 
  iNEXT_data$iNextEst$size_based %>%
  dplyr::rename(names = Assemblage) %>%
  # Fix the samples names 
  # join with metadata 
  left_join(., metadata_df, by = "names") %>%
  # Add colors to data frame
  left_join(., data.frame(cave_colors = cave_colors,
                          cave = names(cave_colors)),
            by = "cave") 
# Inspect 
dim(iNEXT_manual_df)
```

```
## [1] 1200   20
```

```r
str(iNEXT_manual_df)
```

```
## 'data.frame':	1200 obs. of  20 variables:
##  $ names              : chr  "ERR11588428" "ERR11588428" "ERR11588428" "ERR11588428" ...
##  $ m                  : num  1 1180 2360 3540 4719 ...
##  $ Method             : chr  "Rarefaction" "Rarefaction" "Rarefaction" "Rarefaction" ...
##  $ Order.q            : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ qD                 : num  1 533 677 742 779 ...
##  $ qD.LCL             : num  1 530 673 737 774 ...
##  $ qD.UCL             : num  1 536 681 747 784 ...
##  $ SC                 : num  0.00192 0.80536 0.92403 0.96006 0.97554 ...
##  $ SC.LCL             : num  0.0019 0.8034 0.9227 0.9591 0.9748 ...
##  $ SC.UCL             : num  0.00195 0.80735 0.92536 0.96103 0.97626 ...
##  $ input              : num  96845 96845 96845 96845 96845 ...
##  $ filtered           : num  57420 57420 57420 57420 57420 ...
##  $ denoisedF          : num  53646 53646 53646 53646 53646 ...
##  $ denoisedR          : num  52928 52928 52928 52928 52928 ...
##  $ merged             : num  39698 39698 39698 39698 39698 ...
##  $ nochim             : num  21236 21236 21236 21236 21236 ...
##  $ perc_reads_retained: num  21.9 21.9 21.9 21.9 21.9 ...
##  $ cave               : chr  "Fata Apei Cave " "Fata Apei Cave " "Fata Apei Cave " "Fata Apei Cave " ...
##  $ sample             : chr  " FA_1" " FA_1" " FA_1" " FA_1" ...
##  $ cave_colors        : chr  "dodgerblue4" "dodgerblue4" "dodgerblue4" "dodgerblue4" ...
```

```r
# Plot it - Rarefaction Curve 
iNEXT_manual_df %>%
  # Filter out rows that are calcaulted by rarefaction from iNEXT
  dplyr::filter(Method == "Rarefaction") %>%
  # Now, let's draw the plot, be sure to group by the samples if using geom_line!
  ggplot(aes(x = m, y= qD, color = cave, group = names)) + 
  geom_line() + 
  # Facet with the cave to see the samples better 
  facet_grid(Order.q~cave, scales = "free") + 
  scale_color_manual(values = cave_colors) + 
  labs(x = "Number of Sequences (Library Size)", 
       y = "Effective Number of ASVs") + 
  theme_bw() + 
  #guides(color = guide_legend(nrow = 2, byrow = TRUE)) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.title = element_blank()) 
```

<img src="../figures/04_Biodiversity/iNEXT-manual-1.png" style="display: block; margin: auto;" />


While the rarefaction curves from `ggiNEXT` agglomerated the data together, it was hard to see individual samples. Above, we can start to tease apart specific samples within each cave that appear higher. 

*Please interpret your results!!*

Some take-aways:  

- Richness:  
    - The richness is higher in Feta Apei 
    - all the other caves have similar but lower richness
    - in the Ferice cave the samples appear to have an "upper" and "lower" limit and split across ASVs. 
    - There is one line for the Nestor and Taosaure cave
    
- Shannon:
  - Feta Apei had the highest number of common taxa.
  - Then Nestor cave had the second highest amount of common taxa
  - The lowest amount of common taxa was seen in Taosaure cave. 
  - The Ferice cave is definitely on the lower end but there is a portion where samples have high amounts of common taxa but there are also samples that plateu early on indicating low number of common taxa in less ASVs (around 100 ASVs according to the bar plot)
    
- Simpson's Diversity:  
    - Again Feta Apei had the highest simpsons diversity meaning that there are more rare taxa 
    - This was followed by Nestor cave with the second highest Simpson's index
    - Ferice cave has samples that are higher which indicate ASVs with more rare taxa around the 200ish ASV amount. But the samples with less ASVs plateu so quickly and have a really flat line indicating lessss rare taxa. 
    - Taosaure cave also had a lowly abundant amount of rare taxa
    
# Boxplots of Diversity 

Now, let's zoom in on a bit of a better view of the data by visualizing it in a boxplot...  


```r
# Boxplots by station 
iNEXT_manual_df %>%
  ggplot(aes(x = cave, y = qD, fill = cave, color = cave)) + 
  facet_wrap(~Order.q, scales = "free") + 
  geom_jitter(size = 0.5, alpha = 0.5) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = cave_colors) + 
  scale_fill_manual(values = cave_colors) + 
  labs(y = "Effective Number of ASVs") + 
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
```

<img src="../figures/04_Biodiversity/div-boxplot-1.png" style="display: block; margin: auto;" />

We can draw that:  

- Across the board Fata Apei has the highest Richness, Shannon, and Simpson indeces 


- Richness is different between caves it looks significant but we need to tst. Fata Apei has the highest richness. The other caves have similar but lower richness. 


- Ferice cave apears to have the most variation between samples but are lower in Richness, Shannon, and Simpson

 
# Environmental Variable Check 
We do not have environmental data we only have Cave data 


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
##  package          * version    date (UTC) lib source
##  ade4               1.7-22     2023-02-06 [1] CRAN (R 4.3.2)
##  ape                5.7-1      2023-03-13 [2] CRAN (R 4.3.2)
##  Biobase            2.62.0     2023-10-24 [2] Bioconductor
##  BiocGenerics       0.48.1     2023-11-01 [2] Bioconductor
##  biomformat         1.30.0     2023-10-24 [1] Bioconductor
##  Biostrings         2.70.1     2023-10-25 [2] Bioconductor
##  bitops             1.0-7      2021-04-24 [2] CRAN (R 4.3.2)
##  bslib              0.5.1      2023-08-11 [2] CRAN (R 4.3.2)
##  cachem             1.0.8      2023-05-01 [2] CRAN (R 4.3.2)
##  callr              3.7.3      2022-11-02 [2] CRAN (R 4.3.2)
##  cli                3.6.1      2023-03-23 [2] CRAN (R 4.3.2)
##  cluster            2.1.4      2022-08-22 [2] CRAN (R 4.3.2)
##  codetools          0.2-19     2023-02-01 [2] CRAN (R 4.3.2)
##  colorspace         2.1-0      2023-01-23 [2] CRAN (R 4.3.2)
##  crayon             1.5.2      2022-09-29 [2] CRAN (R 4.3.2)
##  data.table         1.14.8     2023-02-17 [2] CRAN (R 4.3.2)
##  devtools         * 2.4.4      2022-07-20 [2] CRAN (R 4.2.1)
##  digest             0.6.33     2023-07-07 [2] CRAN (R 4.3.2)
##  dplyr            * 1.1.3      2023-09-03 [2] CRAN (R 4.3.2)
##  ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
##  evaluate           0.23       2023-11-01 [2] CRAN (R 4.3.2)
##  fansi              1.0.5      2023-10-08 [2] CRAN (R 4.3.2)
##  farver             2.1.1      2022-07-06 [2] CRAN (R 4.3.2)
##  fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
##  forcats          * 1.0.0      2023-01-29 [1] CRAN (R 4.3.2)
##  foreach            1.5.2      2022-02-02 [2] CRAN (R 4.3.2)
##  fs                 1.6.3      2023-07-20 [2] CRAN (R 4.3.2)
##  generics           0.1.3      2022-07-05 [2] CRAN (R 4.3.2)
##  GenomeInfoDb       1.38.0     2023-10-24 [2] Bioconductor
##  GenomeInfoDbData   1.2.11     2023-11-07 [2] Bioconductor
##  ggplot2          * 3.5.0      2024-02-23 [2] CRAN (R 4.3.2)
##  glue               1.6.2      2022-02-24 [2] CRAN (R 4.3.2)
##  gtable             0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
##  highr              0.10       2022-12-22 [2] CRAN (R 4.3.2)
##  hms                1.1.3      2023-03-21 [1] CRAN (R 4.3.2)
##  htmltools          0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
##  htmlwidgets        1.6.2      2023-03-17 [2] CRAN (R 4.3.2)
##  httpuv             1.6.12     2023-10-23 [2] CRAN (R 4.3.2)
##  igraph             1.5.1      2023-08-10 [2] CRAN (R 4.3.2)
##  iNEXT            * 3.0.1      2024-03-24 [1] CRAN (R 4.3.2)
##  IRanges            2.36.0     2023-10-24 [2] Bioconductor
##  iterators          1.0.14     2022-02-05 [2] CRAN (R 4.3.2)
##  jquerylib          0.1.4      2021-04-26 [2] CRAN (R 4.3.2)
##  jsonlite           1.8.7      2023-06-29 [2] CRAN (R 4.3.2)
##  knitr              1.45       2023-10-30 [2] CRAN (R 4.3.2)
##  labeling           0.4.3      2023-08-29 [2] CRAN (R 4.3.2)
##  later              1.3.1      2023-05-02 [2] CRAN (R 4.3.2)
##  lattice            0.21-9     2023-10-01 [2] CRAN (R 4.3.2)
##  lifecycle          1.0.3      2022-10-07 [2] CRAN (R 4.3.2)
##  lubridate        * 1.9.3      2023-09-27 [1] CRAN (R 4.3.2)
##  magrittr           2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
##  MASS               7.3-60     2023-05-04 [2] CRAN (R 4.3.2)
##  Matrix             1.6-1.1    2023-09-18 [2] CRAN (R 4.3.2)
##  memoise            2.0.1      2021-11-26 [2] CRAN (R 4.3.2)
##  mgcv               1.9-0      2023-07-11 [2] CRAN (R 4.3.2)
##  mime               0.12       2021-09-28 [2] CRAN (R 4.3.2)
##  miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
##  multtest           2.58.0     2023-10-24 [1] Bioconductor
##  munsell            0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
##  nlme               3.1-163    2023-08-09 [2] CRAN (R 4.3.2)
##  pacman             0.5.1      2019-03-11 [1] CRAN (R 4.3.2)
##  patchwork        * 1.2.0.9000 2024-03-11 [1] Github (thomasp85/patchwork@d943757)
##  permute            0.9-7      2022-01-27 [1] CRAN (R 4.3.2)
##  phyloseq         * 1.46.0     2023-10-24 [1] Bioconductor
##  pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
##  pkgbuild           1.4.2      2023-06-26 [2] CRAN (R 4.3.2)
##  pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
##  pkgload            1.3.3      2023-09-22 [2] CRAN (R 4.3.2)
##  plyr               1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
##  prettyunits        1.2.0      2023-09-24 [2] CRAN (R 4.3.2)
##  processx           3.8.2      2023-06-30 [2] CRAN (R 4.3.2)
##  profvis            0.3.8      2023-05-02 [2] CRAN (R 4.3.2)
##  promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
##  ps                 1.7.5      2023-04-18 [2] CRAN (R 4.3.2)
##  purrr            * 1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
##  R6                 2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
##  Rcpp               1.0.11     2023-07-06 [2] CRAN (R 4.3.2)
##  RCurl              1.98-1.14  2024-01-09 [1] CRAN (R 4.3.2)
##  readr            * 2.1.5      2024-01-10 [1] CRAN (R 4.3.2)
##  remotes            2.4.2.1    2023-07-18 [2] CRAN (R 4.3.2)
##  reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
##  rhdf5              2.46.1     2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
##  rhdf5filters       1.14.1     2023-11-06 [1] Bioconductor
##  Rhdf5lib           1.24.2     2024-02-07 [1] Bioconductor 3.18 (R 4.3.2)
##  rlang              1.1.2      2023-11-04 [2] CRAN (R 4.3.2)
##  rmarkdown          2.25       2023-09-18 [2] CRAN (R 4.3.2)
##  rstudioapi         0.15.0     2023-07-07 [2] CRAN (R 4.3.2)
##  S4Vectors          0.40.1     2023-10-26 [2] Bioconductor
##  sass               0.4.7      2023-07-15 [2] CRAN (R 4.3.2)
##  scales             1.3.0      2023-11-28 [2] CRAN (R 4.3.2)
##  sessioninfo        1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
##  shiny              1.7.5.1    2023-10-14 [2] CRAN (R 4.3.2)
##  stringi            1.7.12     2023-01-11 [2] CRAN (R 4.3.2)
##  stringr          * 1.5.0      2022-12-02 [2] CRAN (R 4.3.2)
##  survival           3.5-7      2023-08-14 [2] CRAN (R 4.3.2)
##  tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
##  tidyr            * 1.3.0      2023-01-24 [2] CRAN (R 4.3.2)
##  tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.3.2)
##  tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.3.2)
##  timechange         0.3.0      2024-01-18 [1] CRAN (R 4.3.2)
##  tzdb               0.4.0      2023-05-12 [1] CRAN (R 4.3.2)
##  urlchecker         1.0.1      2021-11-30 [2] CRAN (R 4.3.2)
##  usethis          * 2.2.2      2023-07-06 [2] CRAN (R 4.3.2)
##  utf8               1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
##  vctrs              0.6.4      2023-10-12 [2] CRAN (R 4.3.2)
##  vegan              2.6-4      2022-10-11 [1] CRAN (R 4.3.2)
##  withr              2.5.2      2023-10-30 [2] CRAN (R 4.3.2)
##  xfun               0.41       2023-11-01 [2] CRAN (R 4.3.2)
##  xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
##  XVector            0.42.0     2023-10-24 [2] Bioconductor
##  yaml               2.3.7      2023-01-23 [2] CRAN (R 4.3.2)
##  zlibbioc           1.48.0     2023-10-24 [2] Bioconductor
## 
##  [1] /home/sna49/R/x86_64-pc-linux-gnu-library/4.3
##  [2] /programs/R-4.3.2/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```

