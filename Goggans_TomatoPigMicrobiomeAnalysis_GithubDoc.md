Goggans et al., INSERT JOURNAL 2022, Microbiome Analysis Code
================
Mallory Goggans, Cristian Quiroz-Moreno, Emma Bilbrey and Jessica
Cooperstone

# Introduction

INSERT ABSTRACT

Our final taxa used in this anaylsis: - included from Bacteria, Archaea,
Eukaryota and viruses - removed Chordata, Arthropoda, Cnidaria,
Porifera, Echinodermata, Streptophyta, Platyhelminthes because they are
implausible in our biological system (pig gut microbiome) - remove
genera/phyla that have more than 20 zeroes/33.33% missing values across
our dataset of n=60

In the end: Our final dataset has 45 phyla and 755 genera.

### Load libraries

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages('devtools')

BiocManager::install("ALDEx2")
```

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'ALDEx2'

``` r
BiocManager::install("phyloseq")
```

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'phyloseq'

``` r
remotes::install_github("cpauvert/psadd")
BiocManager::install("multtest")
```

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'multtest'

``` r
library(devtools)

devtools::install_github("gauravsk/ranacapa")
```

``` r
# analysis packages
library(ALDEx2) # for univariate analysis
library(rstatix) # for ANOVA
library(vegan) # for beta and alpha diversity
```

    ## Warning: package 'permute' was built under R version 4.0.5

``` r
library(phyloseq) # for krona plots and rarefaction curves
library(psadd) # additions to phyloseq package for microbiome analysis
library(ranacapa) # Utility Functions  for Simple Environmental Visualizations

# functionality packages
library(data.table) # for nicer transposing
library(here) # for directory management
library(knitr) # for knitting and for kable()
library(tidyverse) # for wrangling and plotting
```

    ## Warning: package 'tidyr' was built under R version 4.0.5

    ## Warning: package 'dplyr' was built under R version 4.0.5

``` r
library(readxl) # for reading Excel files
```

### Set seed

Some of our analyses include permutations, so let’s set a seed so we get
consistent results each time we run.

``` r
set.seed(2021) # hoping this seed is better than 2020 :)
```

### Read in metadata

Input files can be found as supplementary information in:

-   UPDATE WITH MALLORY’S PAPER INFO

The data read in chunk below enables loading our data without any
outside-of-R handling. In “Metadata” tab of Supplementary Information.

``` r
# upload metadata
AllSamples.Metadata <- read_excel("Goggans_etal_2021_tomato_pig_microbiome_WGS.xlsx",
                                       sheet = "TableS2.SampleMetadata")

str(AllSamples.Metadata)
```

    ## tibble [60 × 5] (S3: tbl_df/tbl/data.frame)
    ##  $ Sample_Name       : chr [1:60] "ShotgunWGS-ControlPig6GutMicrobiome-Day14" "ShotgunWGS-ControlPig8GutMicrobiome-Day0" "ShotgunWGS-ControlPig3GutMicrobiome-Day14" "ShotgunWGS-TomatoPig14GutMicrobiome-Day7" ...
    ##  $ Pig               : num [1:60] 6 8 3 14 5 18 16 10 2 18 ...
    ##  $ Diet              : chr [1:60] "Control" "Control" "Control" "Tomato" ...
    ##  $ Time_Point        : chr [1:60] "Day 14" "Day 0" "Day 14" "Day 7" ...
    ##  $ Diet_By_Time_Point: chr [1:60] "Control Day 14" "Control Day 0" "Control Day 14" "Tomato Day 7" ...

``` r
# convert Pig, Diet, Time_Point, Diet_By_Time_Point to factors
# and set levels/order
AllSamples.Metadata$Pig <- as.factor(AllSamples.Metadata$Pig)
AllSamples.Metadata$Diet <- as.factor(AllSamples.Metadata$Diet)
AllSamples.Metadata$Time_Point <- factor(AllSamples.Metadata$Time_Point,
                                         levels = c("Day 0", "Day 7", "Day 14"))
AllSamples.Metadata$Diet_By_Time_Point <- 
  factor(AllSamples.Metadata$Diet_By_Time_Point,
         levels = c("Control Day 0", 
                  "Control Day 7", 
                  "Control Day 14", 
                  "Tomato Day 0", 
                  "Tomato Day 7", 
                  "Tomato Day 14"))

# check
str(AllSamples.Metadata)
```

    ## tibble [60 × 5] (S3: tbl_df/tbl/data.frame)
    ##  $ Sample_Name       : chr [1:60] "ShotgunWGS-ControlPig6GutMicrobiome-Day14" "ShotgunWGS-ControlPig8GutMicrobiome-Day0" "ShotgunWGS-ControlPig3GutMicrobiome-Day14" "ShotgunWGS-TomatoPig14GutMicrobiome-Day7" ...
    ##  $ Pig               : Factor w/ 20 levels "1","2","3","4",..: 6 8 3 14 5 18 16 10 2 18 ...
    ##  $ Diet              : Factor w/ 2 levels "Control","Tomato": 1 1 1 2 1 2 2 1 1 2 ...
    ##  $ Time_Point        : Factor w/ 3 levels "Day 0","Day 7",..: 3 1 3 2 2 2 2 2 1 1 ...
    ##  $ Diet_By_Time_Point: Factor w/ 6 levels "Control Day 0",..: 3 1 3 5 2 5 5 2 1 4 ...

# Genera-level annotation

Read in genera level data, annotated from MG-RAST. In “Genera” tab of
Supplementary Information.

``` r
Genus.AllSamples.Counts <- read_excel("Goggans_etal_2021_tomato_pig_microbiome_WGS.xlsx",
                                       sheet = "TableS4.Genera")

str(Genus.AllSamples.Counts)
```

    ## tibble [1,085 × 66] (S3: tbl_df/tbl/data.frame)
    ##  $ domain                                    : chr [1:1085] "Viruses" "Bacteria" "Eukaryota" "Bacteria" ...
    ##  $ phylum                                    : chr [1:1085] "unclassified (derived from Viruses)" "Firmicutes" "unclassified (derived from Eukaryota)" "Cyanobacteria" ...
    ##  $ class                                     : chr [1:1085] "unclassified (derived from Viruses)" "Bacilli" "unclassified (derived from Eukaryota)" "unclassified (derived from Cyanobacteria)" ...
    ##  $ order                                     : chr [1:1085] "Caudovirales" "Lactobacillales" "unclassified (derived from Eukaryota)" "unclassified (derived from Cyanobacteria)" ...
    ##  $ family                                    : chr [1:1085] "Podoviridae" "Aerococcaceae" "unclassified (derived from Eukaryota)" "unclassified (derived from Cyanobacteria)" ...
    ##  $ genus                                     : chr [1:1085] "AHJD-like viruses" "Abiotrophia" "Acanthamoeba" "Acaryochloris" ...
    ##  $ ShotgunWGS-ControlPig6GutMicrobiome-Day14 : num [1:1085] 29 5067 0 271 1988 ...
    ##  $ ShotgunWGS-ControlPig8GutMicrobiome-Day0  : num [1:1085] 0 5661 0 416 2981 ...
    ##  $ ShotgunWGS-ControlPig3GutMicrobiome-Day14 : num [1:1085] 153 4117 0 267 2071 ...
    ##  $ ShotgunWGS-TomatoPig14GutMicrobiome-Day7  : num [1:1085] 0 1576 1 131 1012 ...
    ##  $ ShotgunWGS-ControlPig5GutMicrobiome-Day7  : num [1:1085] 14 3708 0 230 1991 ...
    ##  $ ShotgunWGS-TomatoPig18GutMicrobiome-Day7  : num [1:1085] 1 1159 0 146 585 ...
    ##  $ ShotgunWGS-TomatoPig16GutMicrobiome-Day7  : num [1:1085] 2 2495 0 133 1538 ...
    ##  $ ShotgunWGS-ControlPig10GutMicrobiome-Day7 : num [1:1085] 0 1636 0 141 812 ...
    ##  $ ShotgunWGS-ControlPig2GutMicrobiome-Day0  : num [1:1085] 0 4534 0 338 2670 ...
    ##  $ ShotgunWGS-TomatoPig18GutMicrobiome-Day0  : num [1:1085] 0 2964 1 272 1665 ...
    ##  $ ShotgunWGS-ControlPig10GutMicrobiome-Day0 : num [1:1085] 0 3197 0 264 1411 ...
    ##  $ ShotgunWGS-ControlPig7GutMicrobiome-Day0  : num [1:1085] 0 2513 0 263 1652 ...
    ##  $ ShotgunWGS-ControlPig8GutMicrobiome-Day14 : num [1:1085] 342 4231 0 274 1795 ...
    ##  $ ShotgunWGS-TomatoPig11GutMicrobiome-Day0  : num [1:1085] 0 3101 0 237 2160 ...
    ##  $ ShotgunWGS-TomatoPig19GutMicrobiome-Day0  : num [1:1085] 0 3274 0 228 1729 ...
    ##  $ ShotgunWGS-TomatoPig17GutMicrobiome-Day14 : num [1:1085] 6 1424 0 83 683 ...
    ##  $ ShotgunWGS-ControlPig9GutMicrobiome-Day14 : num [1:1085] 131 3337 0 328 1722 ...
    ##  $ ShotgunWGS-ControlPig10GutMicrobiome-Day14: num [1:1085] 86 3383 0 238 1976 ...
    ##  $ ShotgunWGS-TomatoPig19GutMicrobiome-Day7  : num [1:1085] 0 1849 0 120 940 ...
    ##  $ ShotgunWGS-ControlPig5GutMicrobiome-Day14 : num [1:1085] 76 3864 0 363 2395 ...
    ##  $ ShotgunWGS-ControlPig2GutMicrobiome-Day7  : num [1:1085] 1 5590 0 306 4493 ...
    ##  $ ShotgunWGS-ControlPig6GutMicrobiome-Day7  : num [1:1085] 0 3120 0 201 1273 ...
    ##  $ ShotgunWGS-TomatoPig12GutMicrobiome-Day0  : num [1:1085] 0 2599 0 190 1451 ...
    ##  $ ShotgunWGS-TomatoPig14GutMicrobiome-Day0  : num [1:1085] 1 1453 0 70 846 ...
    ##  $ ShotgunWGS-ControlPig7GutMicrobiome-Day14 : num [1:1085] 67 2906 0 248 1870 ...
    ##  $ ShotgunWGS-TomatoPig11GutMicrobiome-Day14 : num [1:1085] 12 973 0 79 542 16 186 0 185 82 ...
    ##  $ ShotgunWGS-TomatoPig20GutMicrobiome-Day0  : num [1:1085] 0 3682 0 211 2232 ...
    ##  $ ShotgunWGS-ControlPig9GutMicrobiome-Day0  : num [1:1085] 2 2717 1 160 1547 ...
    ##  $ ShotgunWGS-TomatoPig11GutMicrobiome-Day7  : num [1:1085] 0 375 0 31 227 7 69 0 89 29 ...
    ##  $ ShotgunWGS-TomatoPig13GutMicrobiome-Day7  : num [1:1085] 0 2158 0 159 1774 ...
    ##  $ ShotgunWGS-TomatoPig17GutMicrobiome-Day0  : num [1:1085] 0 1409 0 197 762 ...
    ##  $ ShotgunWGS-TomatoPig19GutMicrobiome-Day14 : num [1:1085] 89 1059 0 81 580 ...
    ##  $ ShotgunWGS-TomatoPig13GutMicrobiome-Day0  : num [1:1085] 1 3634 0 207 2188 ...
    ##  $ ShotgunWGS-ControlPig2GutMicrobiome-Day14 : num [1:1085] 106 6111 0 386 3446 ...
    ##  $ ShotgunWGS-ControlPig1GutMicrobiome-Day7  : num [1:1085] 0 3815 1 190 1775 ...
    ##  $ ShotgunWGS-TomatoPig15GutMicrobiome-Day7  : num [1:1085] 1 1126 0 67 974 ...
    ##  $ ShotgunWGS-TomatoPig15GutMicrobiome-Day0  : num [1:1085] 0 3134 0 224 2207 ...
    ##  $ ShotgunWGS-TomatoPig12GutMicrobiome-Day7  : num [1:1085] 0 2376 0 144 1437 ...
    ##  $ ShotgunWGS-TomatoPig14GutMicrobiome-Day14 : num [1:1085] 0 1079 0 61 469 ...
    ##  $ ShotgunWGS-TomatoPig20GutMicrobiome-Day14 : num [1:1085] 47 926 0 61 486 18 205 0 193 41 ...
    ##  $ ShotgunWGS-ControlPig1GutMicrobiome-Day0  : num [1:1085] 0 5545 0 310 2638 ...
    ##  $ ShotgunWGS-ControlPig4GutMicrobiome-Day14 : num [1:1085] 102 3677 0 270 1919 ...
    ##  $ ShotgunWGS-ControlPig6GutMicrobiome-Day0  : num [1:1085] 4 2687 0 200 1640 ...
    ##  $ ShotgunWGS-TomatoPig16GutMicrobiome-Day0  : num [1:1085] 0 2959 0 176 1599 ...
    ##  $ ShotgunWGS-TomatoPig16GutMicrobiome-Day14 : num [1:1085] 3 973 0 98 609 34 222 0 255 149 ...
    ##  $ ShotgunWGS-TomatoPig18GutMicrobiome-Day14 : num [1:1085] 11 1075 0 86 446 ...
    ##  $ ShotgunWGS-ControlPig7GutMicrobiome-Day7  : num [1:1085] 0 1587 0 103 654 ...
    ##  $ ShotgunWGS-ControlPig4GutMicrobiome-Day7  : num [1:1085] 0 1709 0 165 1059 ...
    ##  $ ShotgunWGS-TomatoPig13GutMicrobiome-Day14 : num [1:1085] 0 1021 0 74 685 ...
    ##  $ ShotgunWGS-ControlPig8GutMicrobiome-Day7  : num [1:1085] 12 3035 0 259 1579 ...
    ##  $ ShotgunWGS-TomatoPig15GutMicrobiome-Day14 : num [1:1085] 17 1660 0 140 892 35 366 0 410 224 ...
    ##  $ ShotgunWGS-TomatoPig12GutMicrobiome-Day14 : num [1:1085] 19 2138 0 129 1334 ...
    ##  $ ShotgunWGS-TomatoPig20GutMicrobiome-Day7  : num [1:1085] 0 1699 0 121 645 ...
    ##  $ ShotgunWGS-ControlPig1GutMicrobiome-Day14 : num [1:1085] 14 3895 0 338 2158 ...
    ##  $ ShotgunWGS-ControlPig3GutMicrobiome-Day0  : num [1:1085] 0 4578 0 367 2356 ...
    ##  $ ShotgunWGS-ControlPig5GutMicrobiome-Day0  : num [1:1085] 0 4842 0 274 2894 ...
    ##  $ ShotgunWGS-ControlPig4GutMicrobiome-Day0  : num [1:1085] 0 4439 0 261 2733 ...
    ##  $ ShotgunWGS-ControlPig9GutMicrobiome-Day7  : num [1:1085] 0 703 0 54 384 6 159 0 200 45 ...
    ##  $ ShotgunWGS-ControlPig3GutMicrobiome-Day7  : num [1:1085] 6 4833 0 347 3158 ...
    ##  $ ShotgunWGS-TomatoPig17GutMicrobime-Day7   : num [1:1085] 0 1713 0 136 782 ...

## Data filtering

### Remove inplausible phyla

These phyla are not plausibly found in a rectal swab of a pig, and were
incorrectly annotated, so we are removing them.

``` r
Genus.Counts.Filt <- Genus.AllSamples.Counts %>%
  filter(phylum != "Chordata" , phylum != "Arthropoda" , phylum != "Cnidaria" , 
         phylum != "Porifera" , phylum != "Echinodermata", phylum != "Streptophyta",
         phylum != "Platyhelminthes")
```

Transpose.

``` r
Genus.Counts.Filt.t <- as.tibble(t(Genus.Counts.Filt))
```

    ## Warning: `as.tibble()` was deprecated in tibble 2.0.0.
    ## Please use `as_tibble()` instead.
    ## The signature and semantics have changed, see `?as_tibble`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
    ## Using compatibility `.name_repair`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
# make genus colnames
colnames(Genus.Counts.Filt.t) <- Genus.Counts.Filt.t[6,]

# remove domain, phylum, class, order, family
GenusOnly.Counts.Filt.t <- Genus.Counts.Filt.t[7:66,]

# convert character to numeric
GenusOnly.Counts.Filt.t <- as.data.frame(apply((GenusOnly.Counts.Filt.t), 2, as.numeric))

str(GenusOnly.Counts.Filt.t[,1:5])
```

    ## 'data.frame':    60 obs. of  5 variables:
    ##  $ AHJD-like viruses: num  29 0 153 0 14 1 2 0 0 0 ...
    ##  $ Abiotrophia      : num  5067 5661 4117 1576 3708 ...
    ##  $ Acanthamoeba     : num  0 0 0 1 0 0 0 0 0 1 ...
    ##  $ Acaryochloris    : num  271 416 267 131 230 146 133 141 338 272 ...
    ##  $ Acetivibrio      : num  1988 2981 2071 1012 1991 ...

``` r
# add back sample names as column
GenusOnly.Counts.Filt.t <- GenusOnly.Counts.Filt.t %>%
  mutate(Sample_Name = AllSamples.Metadata$Sample_Name)

# move Sample_Name to first column
GenusOnly.Counts.Filt.t <- GenusOnly.Counts.Filt.t %>%
  relocate(Sample_Name)

kable(head(GenusOnly.Counts.Filt.t))
```

| Sample\_Name                              | AHJD-like viruses | Abiotrophia | Acanthamoeba | Acaryochloris | Acetivibrio | Acetobacter | Acetohalobium | Acholeplasma | Achromobacter | Acidaminococcus | Acidilobus | Acidimicrobium | Acidiphilium | Acidithiobacillus | Acidobacterium | Acidothermus | Acidovorax | Aciduliprofundum | Acinetobacter | Actinobacillus | Actinomyces | Actinosynnema | Aerococcus | Aeromicrobium | Aeromonas | Aeropyrum | Afipia | Aggregatibacter | Agrobacterium | Ahrensia | Ajellomyces | Akkermansia | Albidiferax | Alcanivorax | Algoriphagus | Alicycliphilus | Alicyclobacillus | Aliivibrio | Alistipes | Alkalilimnicola | Alkaliphilus | Allochromatium | Allomyces | Alphabaculovirus | Alphatorquevirus | Alteromonas | Aminobacterium | Aminomonas | Ammonifex | Amycolatopsis | Anabaena | Anaerobaculum | Anaerococcus | Anaerofustis | Anaeromyxobacter | Anaerostipes | Anaerotruncus | Anaplasma | Anoxybacillus | Aquifex | Arcanobacterium | Archaeoglobus | Arcobacter | Aromatoleum | Arthrobacter | Arthroderma | Arthrospira | Ascovirus | Asfivirus | Aspergillus | Asticcacaulis | Atopobium | Aurantimonas | Aureococcus | Avibacterium | Avipoxvirus | Azoarcus | Azorhizobium | Azospirillum | Azotobacter | Babesia | Bacillus | Bacteroides | Bartonella | Basfia | Basidiobolus | Batrachochytrium | Bdellomicrovirus | Bdellovibrio | Beggiatoa | Beijerinckia | Bermanella | Betabaculovirus | Betaentomopoxvirus | Betaretrovirus | Beutenbergia | Bicaudavirus | Bifidobacterium | Bigelowiella | Blastocystis | Blastopirellula | Blattabacterium | Blautia | Bocavirus | Boothiomyces | Bordetella | Borrelia | Botryotinia | Bpp-1-like viruses | Brachybacterium | Brachyspira | Bradyrhizobium | Brevibacillus | Brevibacterium | Brevundimonas | Brucella | Brugia | Bryopsis | Buchnera | Bulleidia | Burkholderia | Butyrivibrio | Caenorhabditis | Cafeteria | Caldanaerobacter | Caldicellulosiruptor | Calditerrivibrio | Caldivirga | Caminibacter | Campylobacter | Candida | Candidatus Accumulibacter | Candidatus Amoebophilus | Candidatus Azobacteroides | Candidatus Blochmannia | Candidatus Carsonella | Candidatus Cloacamonas | Candidatus Desulforudis | Candidatus Hamiltonella | Candidatus Hodgkinia | Candidatus Korarchaeum | Candidatus Koribacter | Candidatus Liberibacter | Candidatus Pelagibacter | Candidatus Phytoplasma | Candidatus Protochlamydia | Candidatus Puniceispirillum | Candidatus Regiella | Candidatus Riesia | Candidatus Solibacter | Candidatus Sulcia | Candidatus Zinderia | Capnocytophaga | Carboxydothermus | Cardiobacterium | Carnobacterium | Catenibacterium | Catenulispora | Catonella | Caulobacter | Cavemovirus | Cellulomonas | Cellulosilyticum | Cellvibrio | Cenarchaeum | Chaetomium | Chelativorans | Chitinophaga | Chlamydia | Chlamydiamicrovirus | Chlamydomonas | Chlamydophila | Chlorella | Chloriridovirus | Chlorobaculum | Chlorobium | Chloroflexus | Chloroherpeton | Chlorovirus | Chondrus | Chromera | Chromobacterium | Chromohalobacter | Chryseobacterium | Chrysodidymus | Chthoniobacter | Citreicella | Citrobacter | Citromicrobium | Cladochytrium | Clavibacter | Clavispora | Clostridium | Coccidioides | Coccolithovirus | Coelomomyces | Collimonas | Collinsella | Colwellia | Comamonas | Conexibacter | Congregibacter | Conidiobolus | Coprinopsis | Coprobacillus | Coprococcus | Coprothermobacter | Coraliomargarita | Corynebacterium | Coxiella | Crinivirus | Croceibacter | Crocosphaera | Cronobacter | Cryptobacterium | Cryptomonas | Cryptosporidium | Cupriavidus | Cyanidioschyzon | Cyanidium | Cyanobium | Cyanophora | Cyanothece | Cylindrospermopsis | Cylindrospermum | Cyprinivirus | Cytomegalovirus | Cytophaga | Debaryomyces | Dechloromonas | Deferribacter | Dehalococcoides | Dehalogenimonas | Deinococcus | Delftia | Denitrovibrio | Dependovirus | Dermacoccus | Desulfarculus | Desulfatibacillum | Desulfitobacterium | Desulfobacterium | Desulfococcus | Desulfohalobium | Desulfomicrobium | Desulfonatronospira | Desulfotalea | Desulfotomaculum | Desulfovibrio | Desulfurispirillum | Desulfurivibrio | Desulfurococcus | Desulfuromonas | Dethiobacter | Dethiosulfovibrio | Dialister | Dichelobacter | Dickeya | Dictyoglomus | Dictyostelium | Dinoroseobacter | Dokdonia | Dorea | Durinskia | Dyadobacter | Ectocarpus | Edwardsiella | Eggerthella | Ehrlichia | Eikenella | Eimeria | Elusimicrobium | Emericella | Emiliania | Encephalitozoon | Endoriftia | Enhydrobacter | Entamoeba | Enterobacter | Enterococcus | Enterocytozoon | Entomophthora | Epsilon15-like viruses | Epulopiscium | Eremococcus | Eremothecium | Erwinia | Erysipelothrix | Erythrobacter | Escherichia | Ethanoligenens | Eubacterium | Euglena | Exiguobacterium | Faecalibacterium | Ferrimonas | Ferroglobus | Ferroplasma | Fervidobacterium | Fibrobacter | Filifactor | Filobasidiella | Finegoldia | Flavobacterium | Floydiella | Fluoribacter | Francisella | Frankia | Fucus | Fulvimarina | Fusobacterium | Gallionella | Gammaretrovirus | Gardnerella | Gemella | Gemmata | Gemmatimonas | Geobacillus | Geobacter | Geodermatophilus | Giardia | Gibberella | Glaciecola | Gloeobacter | Gluconacetobacter | Gluconobacter | Gordonia | Gracilaria | Gracilariopsis | Gramella | Granulibacter | Granulicatella | Guillardia | Haemophilus | Hafnia | Hahella | Halalkalicoccus | Halanaerobium | Haliangium | Haloarcula | Halobacterium | Haloferax | Halogeometricum | Halomicrobium | Halomonas | Haloquadratum | Halorhabdus | Halorhodospira | Halorubrum | Haloterrigena | Halothermothrix | Halothiobacillus | Harpochytrium | Helicobacter | Helicosporidium | Heliobacterium | Hemiselmis | Herbaspirillum | Herminiimonas | Herpetosiphon | Heterosigma | Hirschia | Histophilus | Hoeflea | Holdemania | Hyaloraphidium | Hydrogenivirga | Hydrogenobacter | Hydrogenobaculum | Hyperthermus | Hyphomicrobium | Hyphomonas | Hypocrea | Ichnovirus | Idiomarina | Ignicoccus | Ignisphaera | Ilyobacter | Inovirus | Intrasporangium | Iridovirus | Janibacter | Jannaschia | Janthinobacterium | Jonesia | Jonquetella | Kangiella | Ketogulonicigenium | Kineococcus | Kingella | Klebsiella | Kluyveromyces | Kocuria | Kordia | Kosmotoga | Kribbella | Kryptoperidinium | Ktedonobacter | Kytococcus | L5-like viruses | LUZ24-like viruses | Labrenzia | Laccaria | Lachancea | Lachnum | Lactobacillus | Lactococcus | Lambda-like viruses | Laminaria | Laribacter | Lawsonia | Leadbetterella | Leeuwenhoekiella | Legionella | Leifsonia | Leishmania | Lentisphaera | Leotia | Leptosira | Leptospira | Leptospirillum | Leptothrix | Leptotrichia | Leuconostoc | Limnobacter | Listeria | Listonella | Loa | Lodderomyces | Loktanella | Lutiella | Lymphocryptovirus | Lymphocystivirus | Lyngbya | Lysinibacillus | Macavirus | Macrococcus | Magnaporthe | Magnetococcus | Magnetospirillum | Malassezia | Malawimonas | Mannheimia | Maribacter | Maricaulis | Marinobacter | Marinococcus | Marinomonas | Mariprofundus | Maritimibacter | Marivirga | Megasphaera | Meiothermus | Mesoplasma | Mesorhizobium | Metallosphaera | Methanobrevibacter | Methanocaldococcus | Methanocella | Methanococcoides | Methanococcus | Methanocorpusculum | Methanoculleus | Methanohalobium | Methanohalophilus | Methanoplanus | Methanopyrus | Methanoregula | Methanosaeta | Methanosarcina | Methanosphaera | Methanosphaerula | Methanospirillum | Methanothermobacter | Methanothermococcus | Methanothermus | Methylacidiphilum | Methylibium | Methylobacillus | Methylobacter | Methylobacterium | Methylocella | Methylococcus | Methylophaga | Methylosinus | Methylotenera | Methylovorus | Meyerozyma | Micrococcus | Microcoleus | Microcystis | Micromonas | Micromonospora | Microscilla | Mitsuokella | Mobiluncus | Molluscipoxvirus | Moniliophthora | Monomastix | Monosiga | Moorella | Moraxella | Moritella | Mu-like viruses | Mucilaginibacter | Mycobacterium | Mycoplasma | Myxococcus | N15-like viruses | N4-like viruses | Naegleria | Nakamurella | Nakaseomyces | Nanoarchaeum | Natranaerobius | Natrialba | Natronomonas | Nautilia | Nectria | Neisseria | Neolecta | Neorickettsia | Neosartorya | Nephroselmis | Neptuniibacter | Neurospora | Nitratiruptor | Nitrobacter | Nitrococcus | Nitrosococcus | Nitrosomonas | Nitrosopumilus | Nitrosospira | Nitrospira | Nocardia | Nocardioides | Nocardiopsis | Nodularia | Nosema | Nostoc | Novosphingobium | Oceanibulbus | Oceanicaulis | Oceanicola | Oceanithermus | Oceanobacillus | Ochrobactrum | Ochromonas | Octadecabacter | Odontella | Oedogonium | Oenococcus | Oligotropha | Olsenella | Oltmannsiellopsis | Opitutus | Oribacterium | Orientia | Ornithobacterium | Orthopoxvirus | Oscillatoria | Oscillochloris | Ostreavirus | Ostreococcus | Oxalobacter | P1-like viruses | P2-like viruses | P22-like viruses | Paenibacillus | Paludibacter | Pantoea | Parabacteroides | Parachlamydia | Parachlorella | Paracoccidioides | Paracoccus | Paraglomus | Paramecium | Parascardovia | Parvibaculum | Parvularcula | Pasteurella | Paulinella | Paxillus | Pectobacterium | Pediococcus | Pedobacter | Pelagibaca | Pelobacter | Pelodictyon | Pelotomaculum | Penicillium | Peptoniphilus | Peptostreptococcus | Perkinsus | Persephonella | Petrotoga | Phaeobacter | Phaeodactylum | Phaeosphaeria | Phaeovirus | Phenylobacterium | Phi29-like viruses | PhiC31-like viruses | Phieco32-like viruses | Photobacterium | Photorhabdus | Physoderma | Phytophthora | Pichia | Picrophilus | Pirellula | Planctomyces | Planococcus | Plasmodium | Plectrovirus | Plesiocystis | Podospora | Polaribacter | Polaromonas | Polychytrium | Polynucleobacter | Polysphondylium | Porphyra | Porphyromonas | Postia | Prasinovirus | Prevotella | Prochlorococcus | Propionibacterium | Prosthecochloris | Proteromonas | Proteus | Prototheca | Providencia | Pseudendoclonium | Pseudoalteromonas | Pseudomonas | Pseudoramibacter | Pseudovibrio | Psychrobacter | Psychroflexus | Psychromonas | Pycnococcus | Pylaiella | Pyramidobacter | Pyramimonas | Pyrenophora | Pyrobaculum | Pyrococcus | Pythium | Ralstonia | Ranavirus | Raphidiopsis | Reclinomonas | Reinekea | Renibacterium | Rhadinovirus | Rhizobium | Rhodobacter | Rhodococcus | Rhodomicrobium | Rhodomonas | Rhodopirellula | Rhodopseudomonas | Rhodospirillum | Rhodothermus | Rhopalomyces | Rickettsia | Rickettsiella | Riemerella | Robiginitalea | Roseburia | Roseibium | Roseiflexus | Roseobacter | Roseomonas | Roseovarius | Rothia | Rozella | Rubrobacter | Rudivirus | Ruegeria | Ruminococcus | SP6-like viruses | SPO1-like viruses | SPbeta-like viruses | Saccharomonospora | Saccharomyces | Saccharophagus | Saccharopolyspora | Saccoglossus | Sagittula | Salinibacter | Salinispora | Salmonella | Sanguibacter | Saprolegnia | Scardovia | Scenedesmus | Scheffersomyces | Schizophyllum | Schizosaccharomyces | Sclerotinia | Sebaldella | Segniliparus | Selenomonas | Serratia | Shewanella | Shigella | Shuttleworthia | Sideroxydans | Simonsiella | Simplexvirus | Sinorhizobium | Slackia | Sodalis | Sorangium | Sphaerobacter | Sphingobacterium | Sphingobium | Sphingomonas | Sphingopyxis | Spirochaeta | Spiromicrovirus | Spiroplasma | Spirosoma | Sporosarcina | Stackebrandtia | Staphylococcus | Staphylothermus | Starkeya | Stenotrophomonas | Stigeoclonium | Stigmatella | Streptobacillus | Streptococcus | Streptomyces | Streptosporangium | Subdoligranulum | Sulfitobacter | Sulfolobus | Sulfuricurvum | Sulfurihydrogenibium | Sulfurimonas | Sulfurospirillum | Sulfurovum | Symbiobacterium | Synchytrium | Synechococcus | Synechocystis | Synedra | Syntrophobacter | Syntrophomonas | Syntrophothermus | Syntrophus | T1-like viruses | T4-like viruses | T5-like viruses | T7-like viruses | Talaromyces | Tectivirus | Teredinibacter | Terriglobus | Tetragenococcus | Tetrahymena | Thalassiosira | Thalassobium | Thauera | Theileria | Thermaerobacter | Thermanaerovibrio | Thermincola | Thermoanaerobacter | Thermoanaerobacterium | Thermobaculum | Thermobifida | Thermobispora | Thermococcus | Thermocrinis | Thermodesulfovibrio | Thermofilum | Thermomicrobium | Thermomonospora | Thermoplasma | Thermoproteus | Thermosediminibacter | Thermosinus | Thermosipho | Thermosphaera | Thermosynechococcus | Thermotoga | Thermus | Thioalkalivibrio | Thiobacillus | Thiomicrospira | Thiomonas | Tolumonas | Toxoplasma | Treponema | Trichodesmium | Trichomonas | Trichophyton | Trichoplax | Tropheryma | Truepera | Trypanosoma | Tsukamurella | Tuber | Turicibacter | Uncinocarpus | Ureaplasma | Ustilago | VP2-like phages | Vanderwaltozyma | Varicellovirus | Variovorax | Vaucheria | Veillonella | Verminephrobacter | Verrucomicrobium | Verticillium | Vibrio | Victivallis | Volvox | Vulcanisaeta | Waddlia | Weissella | Whispovirus | Wigglesworthia | Wolbachia | Wolinella | Xanthobacter | Xanthomonas | Xenorhabdus | Xylanimonas | Xylella | Yarrowia | Yatapoxvirus | Yersinia | Zunongwangia | Zygosaccharomyces | Zymomonas | c2-like viruses | phiKMV-like viruses | phiKZ-like viruses | unclassified (derived from Actinobacteria (class)) | unclassified (derived from Alicyclobacillaceae) | unclassified (derived from Alloherpesviridae) | unclassified (derived from Alphaproteobacteria) | unclassified (derived from Alteromonadales) | unclassified (derived from Bacteria) | unclassified (derived from Bacteroidetes) | unclassified (derived from Betaproteobacteria) | unclassified (derived from Burkholderiales) | unclassified (derived from Campylobacterales) | unclassified (derived from Candidatus Poribacteria) | unclassified (derived from Caudovirales) | unclassified (derived from Chromerida) | unclassified (derived from Chroococcales) | unclassified (derived from Clostridiales Family XI. Incertae Sedis) | unclassified (derived from Clostridiales) | unclassified (derived from Deltaproteobacteria) | unclassified (derived from Elusimicrobia) | unclassified (derived from Erysipelotrichaceae) | unclassified (derived from Euryarchaeota) | unclassified (derived from Flavobacteria) | unclassified (derived from Flavobacteriaceae) | unclassified (derived from Flavobacteriales) | unclassified (derived from Fuselloviridae) | unclassified (derived from Gammaproteobacteria) | unclassified (derived from Lachnospiraceae) | unclassified (derived from Marseillevirus family) | unclassified (derived from Methylophilales) | unclassified (derived from Mononegavirales) | unclassified (derived from Myoviridae) | unclassified (derived from Opitutaceae) | unclassified (derived from Pelagophyceae) | unclassified (derived from Phycodnaviridae) | unclassified (derived from Podoviridae) | unclassified (derived from Poxviridae) | unclassified (derived from Proteobacteria) | unclassified (derived from Rhodobacteraceae) | unclassified (derived from Rhodobacterales) | unclassified (derived from Rickettsiales) | unclassified (derived from Ruminococcaceae) | unclassified (derived from Siphoviridae) | unclassified (derived from Thermotogales) | unclassified (derived from Thiotrichales) | unclassified (derived from Verrucomicrobia subdivision 3) | unclassified (derived from Verrucomicrobiales) | unclassified (derived from Vibrionaceae) | unclassified (derived from Vibrionales) | unclassified (derived from Viruses) | unclassified (derived from other sequences) |
|:------------------------------------------|------------------:|------------:|-------------:|--------------:|------------:|------------:|--------------:|-------------:|--------------:|----------------:|-----------:|---------------:|-------------:|------------------:|---------------:|-------------:|-----------:|-----------------:|--------------:|---------------:|------------:|--------------:|-----------:|--------------:|----------:|----------:|-------:|----------------:|--------------:|---------:|------------:|------------:|------------:|------------:|-------------:|---------------:|-----------------:|-----------:|----------:|----------------:|-------------:|---------------:|----------:|-----------------:|-----------------:|------------:|---------------:|-----------:|----------:|--------------:|---------:|--------------:|-------------:|-------------:|-----------------:|-------------:|--------------:|----------:|--------------:|--------:|----------------:|--------------:|-----------:|------------:|-------------:|------------:|------------:|----------:|----------:|------------:|--------------:|----------:|-------------:|------------:|-------------:|------------:|---------:|-------------:|-------------:|------------:|--------:|---------:|------------:|-----------:|-------:|-------------:|-----------------:|-----------------:|-------------:|----------:|-------------:|-----------:|----------------:|-------------------:|---------------:|-------------:|-------------:|----------------:|-------------:|-------------:|----------------:|----------------:|--------:|----------:|-------------:|-----------:|---------:|------------:|-------------------:|----------------:|------------:|---------------:|--------------:|---------------:|--------------:|---------:|-------:|---------:|---------:|----------:|-------------:|-------------:|---------------:|----------:|-----------------:|---------------------:|-----------------:|-----------:|-------------:|--------------:|--------:|--------------------------:|------------------------:|--------------------------:|-----------------------:|----------------------:|-----------------------:|------------------------:|------------------------:|---------------------:|-----------------------:|----------------------:|------------------------:|------------------------:|-----------------------:|--------------------------:|----------------------------:|--------------------:|------------------:|----------------------:|------------------:|--------------------:|---------------:|-----------------:|----------------:|---------------:|----------------:|--------------:|----------:|------------:|------------:|-------------:|-----------------:|-----------:|------------:|-----------:|--------------:|-------------:|----------:|--------------------:|--------------:|--------------:|----------:|----------------:|--------------:|-----------:|-------------:|---------------:|------------:|---------:|---------:|----------------:|-----------------:|-----------------:|--------------:|---------------:|------------:|------------:|---------------:|--------------:|------------:|-----------:|------------:|-------------:|----------------:|-------------:|-----------:|------------:|----------:|----------:|-------------:|---------------:|-------------:|------------:|--------------:|------------:|------------------:|-----------------:|----------------:|---------:|-----------:|-------------:|-------------:|------------:|----------------:|------------:|----------------:|------------:|----------------:|----------:|----------:|-----------:|-----------:|-------------------:|----------------:|-------------:|----------------:|----------:|-------------:|--------------:|--------------:|----------------:|----------------:|------------:|--------:|--------------:|-------------:|------------:|--------------:|------------------:|-------------------:|-----------------:|--------------:|----------------:|-----------------:|--------------------:|-------------:|-----------------:|--------------:|-------------------:|----------------:|----------------:|---------------:|-------------:|------------------:|----------:|--------------:|--------:|-------------:|--------------:|----------------:|---------:|------:|----------:|------------:|-----------:|-------------:|------------:|----------:|----------:|--------:|---------------:|-----------:|----------:|----------------:|-----------:|--------------:|----------:|-------------:|-------------:|---------------:|--------------:|-----------------------:|-------------:|------------:|-------------:|--------:|---------------:|--------------:|------------:|---------------:|------------:|--------:|----------------:|-----------------:|-----------:|------------:|------------:|-----------------:|------------:|-----------:|---------------:|-----------:|---------------:|-----------:|-------------:|------------:|--------:|------:|------------:|--------------:|------------:|----------------:|------------:|--------:|--------:|-------------:|------------:|----------:|-----------------:|--------:|-----------:|-----------:|------------:|------------------:|--------------:|---------:|-----------:|---------------:|---------:|--------------:|---------------:|-----------:|------------:|-------:|--------:|----------------:|--------------:|-----------:|-----------:|--------------:|----------:|----------------:|--------------:|----------:|--------------:|------------:|---------------:|-----------:|--------------:|----------------:|-----------------:|--------------:|-------------:|----------------:|---------------:|-----------:|---------------:|--------------:|--------------:|------------:|---------:|------------:|--------:|-----------:|---------------:|---------------:|----------------:|-----------------:|-------------:|---------------:|-----------:|---------:|-----------:|-----------:|-----------:|------------:|-----------:|---------:|----------------:|-----------:|-----------:|-----------:|------------------:|--------:|------------:|----------:|-------------------:|------------:|---------:|-----------:|--------------:|--------:|-------:|----------:|----------:|-----------------:|--------------:|-----------:|----------------:|-------------------:|----------:|---------:|----------:|--------:|--------------:|------------:|--------------------:|----------:|-----------:|---------:|---------------:|-----------------:|-----------:|----------:|-----------:|-------------:|-------:|----------:|-----------:|---------------:|-----------:|-------------:|------------:|------------:|---------:|-----------:|----:|-------------:|-----------:|---------:|------------------:|-----------------:|--------:|---------------:|----------:|------------:|------------:|--------------:|-----------------:|-----------:|------------:|-----------:|-----------:|-----------:|-------------:|-------------:|------------:|--------------:|---------------:|----------:|------------:|------------:|-----------:|--------------:|---------------:|-------------------:|-------------------:|-------------:|-----------------:|--------------:|-------------------:|---------------:|----------------:|------------------:|--------------:|-------------:|--------------:|-------------:|---------------:|---------------:|-----------------:|-----------------:|--------------------:|--------------------:|---------------:|------------------:|------------:|----------------:|--------------:|-----------------:|-------------:|--------------:|-------------:|-------------:|--------------:|-------------:|-----------:|------------:|------------:|------------:|-----------:|---------------:|------------:|------------:|-----------:|-----------------:|---------------:|-----------:|---------:|---------:|----------:|----------:|----------------:|-----------------:|--------------:|-----------:|-----------:|-----------------:|----------------:|----------:|------------:|-------------:|-------------:|---------------:|----------:|-------------:|---------:|--------:|----------:|---------:|--------------:|------------:|-------------:|---------------:|-----------:|--------------:|------------:|------------:|--------------:|-------------:|---------------:|-------------:|-----------:|---------:|-------------:|-------------:|----------:|-------:|-------:|----------------:|-------------:|-------------:|-----------:|--------------:|---------------:|-------------:|-----------:|---------------:|----------:|-----------:|-----------:|------------:|----------:|------------------:|---------:|-------------:|---------:|-----------------:|--------------:|-------------:|---------------:|------------:|-------------:|------------:|----------------:|----------------:|-----------------:|--------------:|-------------:|--------:|----------------:|--------------:|--------------:|-----------------:|-----------:|-----------:|-----------:|--------------:|-------------:|-------------:|------------:|-----------:|---------:|---------------:|------------:|-----------:|-----------:|-----------:|------------:|--------------:|------------:|--------------:|-------------------:|----------:|--------------:|----------:|------------:|--------------:|--------------:|-----------:|-----------------:|-------------------:|--------------------:|----------------------:|---------------:|-------------:|-----------:|-------------:|-------:|------------:|----------:|-------------:|------------:|-----------:|-------------:|-------------:|----------:|-------------:|------------:|-------------:|-----------------:|----------------:|---------:|--------------:|-------:|-------------:|-----------:|----------------:|------------------:|-----------------:|-------------:|--------:|-----------:|------------:|-----------------:|------------------:|------------:|-----------------:|-------------:|--------------:|--------------:|-------------:|------------:|----------:|---------------:|------------:|------------:|------------:|-----------:|--------:|----------:|----------:|-------------:|-------------:|---------:|--------------:|-------------:|----------:|------------:|------------:|---------------:|-----------:|---------------:|-----------------:|---------------:|-------------:|-------------:|-----------:|--------------:|-----------:|--------------:|----------:|----------:|------------:|------------:|-----------:|------------:|-------:|--------:|------------:|----------:|---------:|-------------:|-----------------:|------------------:|--------------------:|------------------:|--------------:|---------------:|------------------:|-------------:|----------:|-------------:|------------:|-----------:|-------------:|------------:|----------:|------------:|----------------:|--------------:|--------------------:|------------:|-----------:|-------------:|------------:|---------:|-----------:|---------:|---------------:|-------------:|------------:|-------------:|--------------:|--------:|--------:|----------:|--------------:|-----------------:|------------:|-------------:|-------------:|------------:|----------------:|------------:|----------:|-------------:|---------------:|---------------:|----------------:|---------:|-----------------:|--------------:|------------:|----------------:|--------------:|-------------:|------------------:|----------------:|--------------:|-----------:|--------------:|---------------------:|-------------:|-----------------:|-----------:|----------------:|------------:|--------------:|--------------:|--------:|----------------:|---------------:|-----------------:|-----------:|----------------:|----------------:|----------------:|----------------:|------------:|-----------:|---------------:|------------:|----------------:|------------:|--------------:|-------------:|--------:|----------:|----------------:|------------------:|------------:|-------------------:|----------------------:|--------------:|-------------:|--------------:|-------------:|-------------:|--------------------:|------------:|----------------:|----------------:|-------------:|--------------:|---------------------:|------------:|------------:|--------------:|--------------------:|-----------:|--------:|-----------------:|-------------:|---------------:|----------:|----------:|-----------:|----------:|--------------:|------------:|-------------:|-----------:|-----------:|---------:|------------:|-------------:|------:|-------------:|-------------:|-----------:|---------:|----------------:|----------------:|---------------:|-----------:|----------:|------------:|------------------:|-----------------:|-------------:|-------:|------------:|-------:|-------------:|--------:|----------:|------------:|---------------:|----------:|----------:|-------------:|------------:|------------:|------------:|--------:|---------:|-------------:|---------:|-------------:|------------------:|----------:|----------------:|--------------------:|-------------------:|---------------------------------------------------:|------------------------------------------------:|----------------------------------------------:|------------------------------------------------:|--------------------------------------------:|-------------------------------------:|------------------------------------------:|-----------------------------------------------:|--------------------------------------------:|----------------------------------------------:|----------------------------------------------------:|-----------------------------------------:|---------------------------------------:|------------------------------------------:|--------------------------------------------------------------------:|------------------------------------------:|------------------------------------------------:|------------------------------------------:|------------------------------------------------:|------------------------------------------:|------------------------------------------:|----------------------------------------------:|---------------------------------------------:|-------------------------------------------:|------------------------------------------------:|--------------------------------------------:|--------------------------------------------------:|--------------------------------------------:|--------------------------------------------:|---------------------------------------:|----------------------------------------:|------------------------------------------:|--------------------------------------------:|----------------------------------------:|---------------------------------------:|-------------------------------------------:|---------------------------------------------:|--------------------------------------------:|------------------------------------------:|--------------------------------------------:|-----------------------------------------:|------------------------------------------:|------------------------------------------:|----------------------------------------------------------:|-----------------------------------------------:|-----------------------------------------:|----------------------------------------:|------------------------------------:|--------------------------------------------:|
| ShotgunWGS-ControlPig6GutMicrobiome-Day14 |                29 |        5067 |            0 |           271 |        1988 |          66 |          1036 |          779 |           192 |           50181 |         10 |             59 |          244 |               245 |            365 |          389 |        790 |              217 |           794 |           2621 |        1004 |           193 |        341 |            25 |       851 |        31 |     23 |             368 |           517 |       16 |          11 |        1474 |         342 |         255 |          171 |             81 |             1438 |        520 |      1983 |             413 |        13144 |            198 |         0 |                0 |                8 |         148 |            509 |        101 |       915 |           161 |      473 |           130 |         4600 |         1549 |             1311 |         4289 |          6643 |        22 |          1391 |     554 |             531 |           461 |        513 |         412 |          951 |          14 |          92 |         0 |         0 |         167 |           287 |     35624 |          112 |           0 |            0 |           0 |      286 |          190 |          165 |         192 |       8 |    27258 |      318931 |        222 |    744 |            0 |                0 |                1 |          334 |        81 |           96 |         27 |               0 |                  0 |              0 |          406 |            0 |           40993 |            0 |            0 |             295 |              52 |   17707 |         0 |            0 |        947 |      291 |          19 |                 46 |             301 |        5181 |            712 |          1770 |            279 |           126 |      276 |     35 |        1 |       86 |      2321 |         2463 |        40314 |            135 |         0 |             6815 |                 8872 |              295 |         28 |           33 |          4599 |      42 |                       203 |                     196 |                       702 |                     51 |                     0 |                     63 |                    1459 |                      48 |                    0 |                     75 |                   671 |                      23 |                      61 |                    169 |                       177 |                          46 |                  18 |                 7 |                  1612 |                14 |                   0 |           2284 |             4439 |              72 |            232 |            4894 |           250 |       177 |         815 |           0 |          291 |             2656 |        581 |          15 |         30 |           274 |         2203 |       137 |                   1 |           127 |           135 |         0 |               1 |           662 |       2359 |          907 |            468 |           4 |        0 |        0 |             402 |              350 |              328 |             0 |             85 |          13 |         456 |             16 |             0 |         408 |         17 |      297762 |           12 |               0 |            0 |          0 |       19115 |       257 |       119 |          324 |            143 |            0 |          26 |          1405 |       24605 |               336 |              203 |            6812 |      194 |          0 |          475 |          245 |         371 |            3943 |           1 |              26 |         847 |               8 |        25 |        12 |          5 |       1558 |                 55 |               0 |            0 |               0 |      1635 |           32 |           358 |           555 |            3574 |             330 |        1030 |     221 |           644 |            0 |          31 |           309 |               591 |              14579 |              526 |           613 |             264 |              575 |                  70 |          755 |            11294 |          4843 |                331 |             226 |              17 |            644 |          426 |               992 |     24138 |           330 |     453 |         1217 |           175 |             136 |      366 | 20388 |         0 |        1864 |          0 |          296 |       10647 |        64 |        28 |       0 |            386 |         54 |         0 |               2 |         13 |            32 |        99 |          615 |        10711 |             47 |             0 |                      4 |          540 |         116 |           36 |     282 |            287 |           333 |        1977 |          17513 |      276925 |       3 |            2132 |           134155 |        132 |         126 |          95 |              785 |        4768 |        594 |             73 |       2960 |           4955 |          0 |            0 |         419 |     950 |     0 |          40 |          7250 |         131 |               1 |        2547 |     156 |      92 |          211 |        7035 |      5443 |              207 |      24 |         96 |          1 |         470 |               112 |           247 |      125 |          2 |              0 |     1415 |           162 |            478 |         28 |        1195 |      0 |     433 |              43 |           915 |        283 |         87 |            64 |        52 |              61 |            39 |       121 |            59 |          58 |            226 |         41 |            76 |            2226 |               98 |             0 |         2189 |               0 |           5259 |          1 |            147 |           222 |           519 |           4 |      122 |         474 |      31 |      13306 |              0 |             69 |             142 |              263 |           43 |             78 |        167 |        0 |          2 |        288 |         43 |          26 |       2406 |        0 |             149 |          0 |        111 |        128 |               225 |     334 |         271 |       109 |                 68 |         359 |       34 |        785 |            38 |     261 |    115 |       371 |       216 |                0 |           140 |        162 |               1 |                  0 |        61 |       13 |        29 |       0 |        178641 |        2714 |                  65 |         0 |        145 |      335 |           1469 |             1324 |        556 |       216 |         85 |          102 |      0 |         2 |        516 |              1 |        217 |         2436 |        1674 |          31 |     4368 |          0 |   8 |           13 |         42 |       81 |                 0 |                0 |      72 |           1157 |         0 |         503 |          58 |           579 |              640 |         26 |           0 |        381 |        784 |        137 |          482 |            0 |         565 |            47 |             88 |       536 |       39571 |         549 |        194 |           308 |             43 |               1912 |                417 |          179 |              400 |          1200 |               1296 |            383 |             128 |               136 |           276 |          105 |           284 |          269 |           1565 |            499 |              213 |              423 |                 350 |                  10 |             45 |               102 |         181 |             204 |            40 |              588 |          198 |           374 |           13 |           24 |           141 |          173 |         24 |         220 |          73 |         344 |         48 |            219 |         309 |      106534 |       5991 |                0 |             15 |          0 |       57 |     4305 |        92 |        50 |               0 |              657 |          1711 |       1430 |        576 |                0 |               2 |        32 |         220 |           71 |            4 |           1577 |        32 |           62 |      158 |      24 |       899 |        0 |            55 |         155 |            1 |             38 |         96 |           391 |         324 |         121 |           471 |          423 |             41 |          256 |        191 |      365 |          426 |          178 |        69 |      0 |    689 |             271 |           21 |           98 |        114 |           158 |           1728 |          274 |          0 |             53 |         8 |          1 |        776 |         144 |     25686 |                 0 |      955 |        19432 |       45 |                0 |             0 |           42 |             79 |           0 |           92 |         514 |               0 |               1 |                0 |          7823 |         4818 |     526 |           24663 |            32 |             0 |                4 |        337 |          0 |         42 |           846 |          257 |           85 |         410 |          5 |        0 |            694 |        1736 |       3321 |         25 |       2655 |        1255 |          4478 |          50 |          1637 |               1371 |        32 |           222 |       677 |          20 |            36 |            39 |          0 |              188 |                  0 |                   0 |                     0 |            797 |          335 |          0 |           63 |      8 |          57 |       158 |          375 |           1 |        138 |            0 |          118 |        17 |          901 |         425 |            0 |              192 |               0 |       21 |          8873 |      5 |            1 |    1029670 |             660 |               872 |              145 |            0 |     331 |          0 |         118 |                0 |               769 |        2976 |             5626 |           50 |           505 |           190 |          502 |           0 |         0 |           1711 |           0 |          18 |          97 |        443 |       0 |       581 |         0 |           25 |            6 |      203 |           141 |            0 |       814 |         604 |         739 |             83 |          4 |            603 |              937 |            533 |          648 |            0 |        215 |             9 |        491 |           732 |     93165 |        17 |        1503 |         364 |         38 |         136 |    332 |       0 |         817 |         0 |      298 |       145083 |                5 |                35 |                  11 |               133 |            54 |            684 |               421 |           26 |        30 |          456 |         399 |        771 |          550 |           1 |       639 |           0 |              49 |            34 |                  57 |          14 |       2577 |           34 |       36152 |      519 |       2720 |      237 |           5045 |          211 |          50 |            1 |           550 |   13349 |     101 |       434 |           540 |              443 |          82 |          290 |          263 |        2362 |               2 |           0 |      2579 |            0 |            178 |           4604 |              30 |       98 |              247 |             0 |         128 |             681 |        110465 |         1855 |               255 |           82351 |            48 |        158 |           140 |                  608 |          359 |              274 |        232 |            3519 |           0 |          2287 |           497 |       0 |             833 |           3600 |             1217 |        980 |               0 |              22 |               0 |               0 |          46 |          0 |            172 |         226 |              11 |          34 |            61 |           21 |     165 |        16 |             946 |               739 |        2853 |               6464 |                  3550 |           554 |          444 |           223 |          490 |           95 |                 359 |          85 |             250 |             198 |          182 |             6 |                 1924 |        2471 |         905 |            14 |                 463 |       2161 |     570 |              311 |          295 |            238 |        99 |       350 |         43 |      3169 |           408 |         458 |            1 |         68 |         57 |      174 |          46 |          140 |     8 |          774 |            9 |        192 |       48 |               0 |              11 |              0 |        183 |         1 |       23909 |               272 |               65 |            7 |   2167 |         663 |     98 |           17 |      71 |       174 |           0 |             10 |       114 |       430 |          237 |         760 |         138 |         674 |     247 |       52 |            0 |     1199 |         1535 |                 9 |       199 |               2 |                   0 |                  4 |                                                103 |                                            1173 |                                             0 |                                              56 |                                          32 |                                  580 |                                       946 |                                             11 |                                         365 |                                            48 |                                                  26 |                                       64 |                                      0 |                                        33 |                                                                 421 |                                     12636 |                                             186 |                                       246 |                                           17626 |                                       340 |                                       602 |                                           668 |                                          261 |                                          0 |                                             409 |                                       16441 |                                                 0 |                                          11 |                                           0 |                                    206 |                                      86 |                                         2 |                                           0 |                                       7 |                                      0 |                                          3 |                                           39 |                                          39 |                                        16 |                                       15985 |                                      768 |                                       115 |                                         0 |                                                       120 |                                            119 |                                       56 |                                      33 |                                 252 |                                          48 |
| ShotgunWGS-ControlPig8GutMicrobiome-Day0  |                 0 |        5661 |            0 |           416 |        2981 |          86 |          1373 |         1269 |           298 |           39909 |         20 |            126 |          344 |               375 |            539 |          554 |       1502 |              498 |          1045 |           2461 |        1433 |           281 |        354 |            68 |      1124 |        33 |     21 |             416 |           811 |       45 |          12 |        2182 |         633 |         333 |          624 |            149 |             1835 |        704 |      5891 |             337 |        19539 |            220 |         0 |                1 |               18 |         287 |           1036 |        186 |      1199 |           221 |      741 |           249 |         6175 |         2355 |             1943 |         5468 |         13087 |        48 |          1884 |     674 |             974 |           540 |        708 |         525 |         1501 |          18 |         126 |         1 |         0 |         145 |           327 |     42713 |          149 |           2 |            0 |           0 |      470 |          250 |          295 |         300 |      14 |    36494 |      442199 |        338 |    862 |            0 |                0 |                4 |          444 |       114 |          164 |         58 |               1 |                  0 |              0 |          601 |            1 |           95976 |            1 |            5 |             514 |             166 |   23637 |         1 |            0 |       1069 |      509 |          40 |                 14 |             427 |        6708 |           1036 |          2244 |            417 |           171 |      394 |     47 |        0 |       96 |      2327 |         3547 |        39907 |            204 |         0 |             9162 |                11870 |              450 |         53 |           66 |          8514 |      74 |                       256 |                     377 |                      1730 |                     69 |                     2 |                    160 |                    1952 |                      57 |                    0 |                     97 |                   899 |                      18 |                      79 |                    254 |                       284 |                          70 |                  24 |                 1 |                  1865 |                49 |                   0 |           3838 |             5382 |              89 |            338 |            7907 |           369 |       272 |         953 |           0 |          413 |             4072 |        586 |          26 |         17 |           457 |         3454 |       128 |                  18 |           150 |           218 |         6 |               0 |          1106 |       2794 |         1107 |            649 |          14 |        0 |        0 |             599 |              500 |              680 |             1 |            375 |          47 |         688 |             32 |             0 |         616 |         25 |      424056 |           55 |               0 |            0 |          1 |       20485 |       391 |       250 |          412 |            169 |            0 |          43 |          2230 |       26769 |               543 |              473 |            3838 |      228 |          0 |          877 |          240 |         470 |            5494 |           1 |              74 |        1257 |               9 |        18 |       383 |         14 |       2011 |                 56 |               1 |            1 |               0 |      2596 |           35 |           927 |           588 |            3698 |             379 |        1387 |     349 |           956 |            0 |          56 |           436 |               981 |              20232 |              744 |           911 |             396 |              762 |                 124 |         1047 |            13009 |          9335 |                416 |             395 |              35 |            773 |          699 |              1632 |     10396 |           320 |     569 |         1554 |           375 |             193 |      776 | 62948 |         0 |        2525 |          0 |          395 |       12983 |       103 |        56 |       0 |            678 |         62 |         4 |              11 |         42 |            67 |      1267 |          763 |        15015 |             76 |             0 |                      9 |          897 |         230 |           42 |     375 |            436 |           451 |        4404 |          24857 |      252845 |       0 |            2864 |           124365 |        205 |         108 |         114 |              968 |        6963 |       1080 |            117 |       3329 |           7259 |          0 |            1 |         750 |    1431 |     0 |          39 |         10721 |         196 |               1 |        5115 |     243 |     179 |          340 |        8831 |      7012 |              275 |      59 |        174 |          2 |         676 |               221 |           274 |      247 |          3 |              0 |     2455 |           225 |            617 |         31 |        1566 |      0 |     567 |              61 |          1200 |        387 |        114 |            90 |        60 |              96 |            66 |       172 |           109 |          90 |            308 |         50 |            94 |            2722 |              156 |             0 |         3541 |               0 |           5937 |          0 |            181 |           349 |           658 |           1 |      117 |         726 |      45 |      18363 |              0 |             96 |             169 |              325 |           46 |            133 |        255 |        0 |          0 |        486 |         60 |          24 |       2995 |        1 |             149 |          3 |        223 |        202 |               371 |     475 |         374 |       162 |                121 |         625 |       69 |        936 |            72 |     291 |    302 |       644 |       346 |                0 |           242 |        240 |               0 |                  1 |       127 |       42 |        33 |       0 |        153576 |        2761 |                 104 |         0 |        193 |      595 |           1719 |             1775 |        680 |       404 |        128 |          272 |      0 |         1 |        809 |              1 |        418 |         3048 |        1911 |          74 |     5766 |          0 |  14 |           25 |         80 |      108 |                 0 |                0 |     113 |           1562 |         0 |         601 |          91 |           773 |              826 |         35 |           0 |        387 |       1442 |        258 |          703 |            0 |         736 |            85 |            132 |      1023 |       25118 |         864 |        266 |           523 |             57 |               2802 |                696 |          221 |              564 |          1776 |               1748 |            528 |             155 |               232 |           454 |          176 |           350 |          401 |           2403 |            614 |              222 |              509 |                 539 |                  27 |             98 |               207 |         473 |             354 |            84 |              872 |          172 |           534 |           49 |           94 |           231 |          160 |         32 |         285 |         142 |         410 |         74 |            331 |         673 |       91970 |      10473 |                0 |             23 |          0 |       88 |     5010 |       111 |        80 |               4 |              992 |          2469 |       1700 |        910 |                0 |              15 |       112 |         326 |           81 |            5 |           2109 |        73 |           99 |      221 |      35 |      1301 |        0 |            63 |         181 |            0 |             75 |         97 |           503 |         459 |         178 |           765 |          643 |             54 |          349 |        307 |      462 |          851 |          276 |        83 |      1 |    952 |             391 |           30 |          135 |        159 |           243 |           2209 |          488 |          0 |             69 |         3 |          0 |       1108 |         102 |     30551 |                 4 |     1555 |         7440 |       81 |                0 |             0 |           95 |            143 |           0 |          161 |         741 |               2 |              16 |               19 |         11386 |         8899 |     601 |           39080 |            47 |             0 |                7 |        338 |          0 |         71 |          1019 |          366 |          141 |         528 |         15 |        0 |            871 |        2258 |       5681 |         35 |       3491 |        1261 |          5463 |          62 |          2644 |               2008 |        43 |           276 |       995 |          26 |            83 |            54 |          0 |              173 |                  8 |                   0 |                     0 |           1130 |          481 |          0 |          156 |     23 |         109 |       290 |          666 |           1 |        236 |            0 |          162 |        40 |         1645 |        1012 |            0 |              395 |               0 |       23 |         16492 |      7 |            4 |     820657 |            1108 |              1128 |              263 |            0 |     427 |          1 |         240 |                0 |              1015 |        3925 |             2424 |           74 |           754 |           328 |          705 |           0 |         1 |           3061 |           0 |          26 |         134 |        626 |       0 |       761 |         0 |           24 |            2 |      349 |           204 |            1 |      1049 |         795 |         991 |            151 |          3 |            983 |             1365 |            772 |          972 |            1 |        306 |            10 |        840 |          1135 |     72078 |        39 |        2141 |         567 |        106 |         243 |    397 |       0 |        1032 |         0 |      532 |       233438 |                0 |                42 |                  11 |               204 |           108 |            842 |               638 |           18 |        53 |          667 |         550 |       1238 |          699 |           1 |       855 |           5 |              87 |            40 |                 164 |          12 |       3363 |           68 |       33861 |      680 |       3781 |      510 |           4244 |          176 |          92 |            0 |           907 |   17333 |     168 |       566 |           856 |             1017 |         188 |          368 |          326 |        4106 |               7 |           1 |      3749 |            0 |            274 |           5633 |              58 |      178 |              421 |             3 |         206 |            1115 |         32512 |         2943 |               509 |           65446 |           102 |        212 |           206 |                  584 |          432 |              343 |        417 |            4358 |           0 |          4236 |           637 |       0 |            1217 |           4845 |             1665 |       1431 |               0 |              85 |               1 |              35 |          27 |          1 |            249 |         414 |              19 |          56 |            98 |           29 |     243 |        19 |            1323 |               951 |        3534 |               8290 |                  4124 |           739 |          580 |           338 |          814 |          130 |                 489 |         121 |             388 |             323 |          300 |             7 |                 2772 |        5000 |        1238 |            12 |                 599 |       3069 |     812 |              445 |          448 |            358 |       133 |       480 |         41 |      5331 |           516 |        1577 |            7 |         68 |         90 |      280 |         121 |          147 |    27 |         1135 |           23 |        241 |       98 |               0 |              23 |              0 |        323 |         1 |       19797 |               508 |              303 |           24 |   3077 |        1561 |    166 |           30 |     152 |       212 |           0 |             15 |       200 |       646 |          319 |        1179 |         169 |         796 |     327 |      109 |            0 |     1352 |         1959 |                27 |       269 |               3 |                   0 |                  6 |                                                209 |                                            1414 |                                             0 |                                             141 |                                          44 |                                  821 |                                      2204 |                                             17 |                                         815 |                                            72 |                                                  49 |                                      171 |                                      0 |                                        47 |                                                                 520 |                                     16819 |                                             355 |                                       329 |                                           28383 |                                       550 |                                      1064 |                                           881 |                                          725 |                                          0 |                                             739 |                                       21035 |                                                 0 |                                          16 |                                           0 |                                    587 |                                     597 |                                         2 |                                           2 |                                      33 |                                      0 |                                          2 |                                           85 |                                          95 |                                        12 |                                       51452 |                                     1069 |                                       191 |                                         0 |                                                       336 |                                            254 |                                       80 |                                      33 |                                 311 |                                          33 |
| ShotgunWGS-ControlPig3GutMicrobiome-Day14 |               153 |        4117 |            0 |           267 |        2071 |          60 |          1015 |          817 |           197 |           31994 |         25 |             81 |          243 |               222 |            309 |          301 |       1126 |              245 |           653 |           1798 |         698 |           181 |        737 |            18 |       762 |        26 |     12 |             289 |           443 |       13 |          11 |        1419 |         464 |         302 |          272 |            141 |             1369 |        486 |      2917 |             325 |        13386 |            173 |         0 |                0 |                5 |         195 |            607 |        112 |       836 |           195 |      485 |           147 |         4189 |         1512 |             1216 |         4240 |          8660 |        44 |          1416 |     469 |             484 |           407 |        569 |         384 |          778 |          14 |          89 |         0 |         1 |         111 |           213 |     19924 |           84 |           3 |            0 |           0 |      338 |          184 |          238 |         153 |       6 |    26183 |      294856 |        170 |    561 |            0 |                0 |                0 |          315 |        50 |           96 |         38 |               0 |                  2 |              0 |          222 |            0 |           25359 |            0 |            1 |             285 |              90 |   18590 |         0 |            0 |        841 |      246 |          21 |                 12 |             223 |        4840 |            659 |          1647 |            205 |           115 |      225 |     30 |        1 |       76 |      1531 |         2412 |        34373 |            104 |         0 |             6622 |                 9318 |              278 |         21 |           37 |          4769 |      55 |                       210 |                     216 |                       790 |                     31 |                     1 |                     87 |                    1519 |                      39 |                    0 |                     52 |                   694 |                      10 |                      80 |                    150 |                       183 |                          65 |                  11 |                 2 |                  1413 |                31 |                   0 |           2428 |             4140 |              69 |            270 |            3721 |           208 |       246 |         804 |           0 |          203 |             2758 |        492 |          21 |         31 |           309 |         2284 |        75 |                   0 |            73 |           163 |         3 |               1 |           689 |       1986 |          898 |            610 |           7 |        0 |        0 |             374 |              336 |              337 |             0 |            168 |          31 |         470 |             27 |             0 |         335 |         39 |      299010 |           18 |               0 |            0 |          0 |       20465 |       323 |       169 |          293 |            117 |            0 |          14 |          1264 |       24714 |               356 |              247 |            1818 |      145 |          0 |          435 |          188 |         322 |            3273 |           0 |              29 |         873 |               6 |        16 |       462 |         10 |       1563 |                 47 |               0 |            1 |               0 |      1595 |           15 |           369 |           381 |            3124 |             313 |         966 |     246 |           760 |            0 |          28 |           299 |               604 |              14708 |              498 |           585 |             338 |              528 |                  71 |          822 |            10274 |          5420 |                271 |             237 |              22 |            589 |          492 |               895 |      7351 |           244 |     440 |         1109 |           153 |             118 |      414 | 32261 |         2 |        1656 |          2 |          238 |        9588 |        59 |        36 |       0 |            452 |         39 |         4 |               2 |         14 |            53 |       113 |          478 |        13222 |             54 |             0 |                      4 |          584 |         103 |           33 |     305 |            209 |           309 |        3066 |          19759 |      269071 |       0 |            2082 |           175626 |        142 |          80 |          59 |              794 |        4985 |        621 |             74 |       2373 |           4924 |          0 |            1 |         509 |     954 |     0 |          26 |          6742 |         161 |               2 |        1795 |     172 |      80 |          241 |        6468 |      5172 |              177 |      24 |        104 |          0 |         536 |               137 |           170 |      103 |          5 |              0 |     1371 |           176 |            656 |         25 |        1008 |      0 |     369 |              44 |           910 |        274 |         86 |            46 |        37 |              57 |            35 |       103 |            78 |          41 |            202 |         49 |            75 |            2050 |              107 |             0 |         2163 |               0 |           5150 |          0 |            206 |           270 |           522 |           0 |       98 |         571 |      27 |      12582 |              0 |             88 |              99 |              182 |           35 |             76 |        198 |        0 |          0 |        326 |         35 |          16 |       2107 |        1 |              96 |          2 |        117 |        143 |               296 |     243 |         211 |       105 |                 83 |         270 |       44 |        600 |            35 |     155 |    167 |       372 |       185 |                1 |           125 |        122 |               1 |                  1 |        61 |       14 |        16 |       0 |        163871 |        4344 |                  49 |         0 |        160 |      372 |           1507 |             1139 |        484 |       187 |         61 |          154 |      0 |         1 |        531 |              3 |        269 |         2470 |        2143 |          61 |     4387 |          0 |  12 |           17 |         62 |       63 |                 0 |                0 |      94 |           1202 |         0 |         520 |          65 |           531 |              494 |         20 |           0 |        270 |        763 |        151 |          470 |            0 |         463 |            60 |            112 |       542 |       27391 |         521 |        194 |           311 |             34 |               1931 |                418 |          146 |              445 |          1226 |               1154 |            430 |             103 |               117 |           321 |          121 |           247 |          268 |           1527 |            380 |              171 |              373 |                 343 |                  38 |             54 |               138 |         284 |             228 |            61 |              631 |          162 |           376 |           32 |           77 |           136 |          106 |         23 |         177 |          73 |         376 |         56 |            222 |         348 |       57206 |       1189 |                0 |             12 |          1 |       55 |     3974 |        79 |        40 |               0 |              441 |          1460 |       1180 |        576 |                0 |               2 |        37 |         192 |           58 |            4 |           1539 |        51 |           68 |      159 |      21 |       814 |        0 |            32 |         113 |            1 |             36 |         48 |           297 |         355 |         126 |           482 |          406 |             38 |          241 |        212 |      272 |          397 |          158 |        63 |      0 |    701 |             366 |           20 |           68 |        114 |           153 |           1575 |          220 |          0 |             52 |         7 |          0 |        876 |          74 |     13557 |                 2 |     1166 |         7969 |       78 |                0 |             0 |           72 |             96 |           0 |          112 |         582 |               7 |               3 |                1 |          8001 |         5183 |     404 |           23919 |            22 |             0 |                3 |        203 |          0 |         40 |           369 |          258 |          114 |         274 |         30 |        0 |            609 |        1780 |       3350 |         21 |       2774 |        1015 |          4192 |          33 |          1488 |               1213 |        19 |           219 |       660 |          18 |            49 |            22 |          0 |              181 |                  1 |                   1 |                     0 |            850 |          288 |          0 |           55 |     13 |          83 |       229 |          397 |           2 |        111 |            0 |           98 |        20 |          906 |         780 |            0 |              271 |               0 |       32 |          8026 |      5 |            1 |     890023 |             964 |               618 |              169 |            0 |     274 |          2 |         147 |                0 |               709 |        2616 |             2232 |           40 |           465 |           169 |          520 |           0 |         0 |           1543 |           0 |           9 |         119 |        418 |       0 |       576 |         1 |           20 |            2 |      210 |            95 |            0 |       684 |         577 |         650 |             88 |          7 |            579 |              910 |            592 |          631 |            0 |        223 |             7 |        453 |           803 |     83235 |        24 |        1635 |         348 |         62 |         126 |    341 |       0 |         697 |         0 |      377 |       131810 |                5 |                33 |                  13 |               112 |            72 |            871 |               333 |           16 |        28 |          464 |         379 |        760 |          391 |           1 |       376 |           0 |              33 |            31 |                  62 |          13 |       2214 |           42 |       37426 |      513 |       2671 |      415 |           3675 |          139 |         113 |            0 |           586 |   11188 |     122 |       423 |           579 |              515 |         108 |          317 |          249 |        2292 |               0 |           0 |      2461 |            0 |            141 |           4510 |              51 |       99 |              215 |             1 |         136 |             823 |        431655 |         1787 |               294 |           75804 |            76 |        111 |           141 |                  500 |          353 |              258 |        266 |            3115 |           0 |          3911 |           496 |       0 |             750 |           3692 |             1165 |        910 |               0 |              36 |               0 |              20 |           9 |          0 |            153 |         247 |              24 |          32 |            73 |           13 |     162 |        30 |             964 |               633 |        2661 |               5809 |                  2962 |           484 |          355 |           177 |          532 |           85 |                 348 |          93 |             304 |             218 |          175 |             8 |                 1960 |        3508 |         917 |            25 |                 474 |       2175 |     514 |              269 |          324 |            250 |       114 |       285 |         14 |      3013 |           358 |         466 |            0 |         31 |         40 |      167 |          65 |           93 |    19 |          641 |           10 |        152 |       28 |               0 |               9 |              0 |        292 |         0 |       17258 |               319 |              140 |           13 |   2111 |         620 |    109 |           17 |     111 |       219 |           0 |             15 |       129 |       368 |          241 |         770 |         129 |         447 |     213 |       42 |            0 |     1067 |         1432 |                12 |       246 |               4 |                   0 |                  7 |                                                104 |                                            1054 |                                             0 |                                              81 |                                          29 |                                  577 |                                       818 |                                              5 |                                         672 |                                            41 |                                                  21 |                                       63 |                                      0 |                                        37 |                                                                 319 |                                     11368 |                                             207 |                                       247 |                                           19120 |                                       311 |                                       604 |                                           592 |                                          345 |                                          0 |                                             517 |                                       16846 |                                                 0 |                                          16 |                                           0 |                                    214 |                                     163 |                                         1 |                                           0 |                                      19 |                                      3 |                                          0 |                                           63 |                                          48 |                                         9 |                                       21377 |                                     1145 |                                       134 |                                         0 |                                                       159 |                                            138 |                                       34 |                                      16 |                                 182 |                                          62 |
| ShotgunWGS-TomatoPig14GutMicrobiome-Day7  |                 0 |        1576 |            1 |           131 |        1012 |          41 |           399 |          391 |            90 |           11378 |          7 |             27 |          100 |               102 |            135 |          113 |        335 |              195 |           442 |           1000 |         405 |            91 |        135 |             8 |      1263 |        10 |     15 |             263 |           230 |       13 |           6 |         578 |         156 |         183 |          104 |             35 |              555 |        391 |      1410 |             132 |         5297 |            110 |         0 |                0 |               24 |         197 |            242 |         55 |       278 |            52 |      225 |            88 |         1576 |          762 |              457 |         1944 |          3531 |        18 |           532 |     170 |             443 |           217 |        243 |         155 |          406 |           4 |          44 |         0 |         0 |          43 |            62 |      3501 |           38 |           0 |            0 |           0 |      159 |           67 |           79 |         108 |       3 |    10048 |      107303 |         79 |    331 |            0 |                0 |                2 |          148 |        50 |           37 |         61 |               0 |                  0 |              0 |          203 |            1 |            5542 |            0 |            1 |             105 |              40 |    7102 |         0 |            0 |        298 |      213 |          29 |                  4 |             171 |        1844 |            240 |           599 |            121 |            46 |      160 |     54 |        0 |      119 |       694 |         1041 |        10853 |             66 |         0 |             2529 |                 3572 |              122 |          7 |           21 |          2070 |      21 |                        84 |                      84 |                       352 |                     32 |                     0 |                     36 |                     529 |                      43 |                    0 |                     36 |                   240 |                       9 |                      32 |                     97 |                        59 |                          21 |                  34 |                 6 |                   425 |                11 |                   1 |            846 |             1489 |              41 |            101 |            1376 |           104 |        84 |         225 |           1 |          129 |             1136 |        192 |          15 |          6 |            89 |          759 |        28 |                   0 |            36 |            50 |         1 |               1 |           255 |        737 |          308 |            216 |           4 |        2 |        0 |             178 |              173 |              137 |             0 |             33 |          12 |         410 |              5 |             0 |         205 |          7 |      109452 |            9 |               0 |            0 |          6 |        3136 |       214 |        62 |          113 |             55 |            0 |          14 |           471 |        5653 |               142 |               82 |            3171 |      101 |          0 |          184 |           96 |         268 |             921 |           0 |              19 |         320 |               3 |         5 |        13 |          3 |        580 |                 23 |               0 |            0 |               0 |       571 |            6 |           193 |           163 |            1136 |             156 |         371 |      81 |           258 |            0 |          24 |           135 |               286 |               5163 |              200 |           254 |             138 |              207 |                  39 |          304 |             3462 |          3867 |                114 |             112 |              13 |            240 |          219 |               392 |      1396 |           110 |     289 |          424 |            94 |              59 |      217 |  8251 |         0 |         560 |          0 |          257 |        2816 |        38 |        23 |       0 |            247 |         28 |         5 |               6 |         13 |            42 |        48 |          448 |         4365 |             17 |             6 |                      3 |          290 |          71 |           17 |     210 |            122 |           111 |       12096 |           7325 |       55035 |       1 |             693 |            47886 |        196 |          35 |          44 |              255 |        1750 |        305 |             34 |       1326 |           1745 |          2 |            0 |         285 |     330 |     0 |          11 |          2992 |          65 |               0 |         511 |     118 |      31 |           61 |        2475 |      1952 |               95 |       7 |         41 |          1 |         176 |                66 |            58 |       69 |          2 |              0 |      622 |            80 |            190 |          9 |         849 |      0 |     217 |              24 |           304 |        134 |         42 |            26 |        23 |              28 |            14 |        68 |            29 |          32 |            110 |         22 |            35 |             731 |               78 |             0 |         3519 |               0 |           1694 |          1 |             93 |           122 |           200 |           0 |       33 |         357 |      22 |       4418 |              0 |             45 |              57 |               75 |           18 |             25 |         68 |        0 |          3 |        301 |         29 |          12 |        713 |        1 |              52 |          0 |         76 |         59 |               138 |     136 |         158 |        80 |                 29 |         226 |       26 |        610 |            19 |     140 |    103 |       183 |       106 |                0 |            83 |         90 |               0 |                  0 |        60 |       10 |         8 |       0 |        104615 |         809 |                 141 |         0 |         90 |      240 |            435 |              510 |        259 |       104 |         37 |           85 |      0 |         1 |        237 |              0 |        106 |          745 |         545 |          41 |     1651 |          0 |  18 |            6 |         25 |       36 |                 0 |                0 |      24 |            466 |         0 |         206 |          29 |           211 |              225 |         16 |           0 |        187 |        405 |         67 |          298 |            0 |         278 |            36 |             36 |       235 |        6151 |         229 |        129 |           110 |             23 |               2746 |                273 |           87 |              179 |           537 |                481 |            153 |              54 |                81 |           122 |           79 |           120 |          124 |            765 |            196 |               90 |              166 |                 361 |                  23 |             58 |                64 |          98 |             107 |            44 |              234 |           48 |           182 |           28 |           22 |            83 |           68 |         10 |          70 |          27 |         116 |         28 |             84 |         153 |        7754 |       5070 |                0 |              5 |          0 |       46 |     1363 |        58 |       136 |               4 |              222 |           718 |        595 |        234 |                0 |               2 |        27 |          97 |           39 |            6 |            599 |        13 |           21 |       71 |      12 |       511 |        0 |            19 |          55 |            2 |             45 |         30 |           131 |         140 |          49 |           256 |          239 |             22 |          116 |         92 |      141 |          219 |           67 |        24 |      1 |    257 |              80 |           12 |           31 |         38 |            72 |            598 |           94 |          0 |             26 |         4 |          0 |        286 |          38 |      2330 |                 1 |      295 |         3460 |       22 |                0 |             0 |           44 |             29 |           0 |           54 |         378 |             126 |               6 |                5 |          2964 |         1842 |     290 |           10529 |            19 |             0 |                1 |         78 |          0 |          9 |            77 |          108 |           51 |         305 |          7 |        0 |            369 |         615 |       1243 |         10 |        949 |         303 |          1301 |          12 |           900 |                614 |        17 |            68 |       305 |          13 |            32 |            16 |          0 |               53 |                  1 |                   0 |                     0 |            686 |          290 |          0 |           38 |      8 |          50 |        54 |          151 |           4 |        127 |            0 |           57 |         6 |          316 |         210 |            0 |              100 |               0 |        9 |          3428 |      1 |            0 |     277930 |             276 |               501 |               54 |            0 |     292 |          0 |         159 |                0 |               652 |        1471 |              630 |           19 |           263 |            79 |          476 |           0 |         1 |            645 |           0 |          13 |          53 |        217 |       0 |       240 |         0 |            7 |            4 |      113 |            69 |            1 |       283 |         189 |         298 |             45 |          1 |            199 |              438 |            285 |          247 |            0 |        113 |            20 |        175 |           318 |     12198 |         4 |         557 |         189 |         22 |          54 |    145 |       0 |         259 |         0 |      147 |        58127 |                1 |                15 |                   3 |                76 |            41 |            271 |               140 |            5 |        11 |          167 |         157 |       1022 |          198 |           0 |        81 |           0 |              20 |            11 |                  47 |           6 |        823 |           18 |       11115 |      335 |       1912 |     1500 |           1089 |           63 |          54 |            0 |           223 |    3372 |     106 |       122 |           269 |              210 |          62 |          120 |           76 |        1179 |               0 |           0 |       809 |            0 |             73 |           3212 |              22 |       57 |              151 |             0 |          46 |             272 |          9055 |          785 |               139 |           25259 |            27 |         84 |            61 |                  190 |          120 |              112 |        118 |            1170 |           2 |           851 |           164 |       0 |             335 |           1360 |              410 |        352 |               0 |              18 |               0 |               5 |           9 |          1 |            104 |          80 |               9 |          15 |            29 |            6 |      66 |        11 |             357 |               290 |         978 |               2398 |                  1040 |           196 |          194 |            87 |          242 |           42 |                 149 |          33 |             111 |             110 |          101 |             3 |                  819 |        1484 |         357 |             5 |                 158 |        805 |     235 |              163 |          116 |            146 |        64 |       670 |          8 |      2161 |           173 |         251 |            2 |         19 |         26 |       60 |          19 |           48 |     6 |          406 |            8 |         80 |       23 |               0 |              11 |              0 |         67 |         0 |        4735 |                84 |               45 |            7 |   1727 |         407 |     59 |            9 |      67 |        73 |           0 |             17 |        74 |       260 |           86 |         363 |         126 |         207 |     121 |       30 |            0 |      772 |          501 |                 4 |        63 |               2 |                   0 |                  6 |                                                 64 |                                             377 |                                             0 |                                              32 |                                          77 |                                  929 |                                       413 |                                             10 |                                         378 |                                            19 |                                                  17 |                                       59 |                                      0 |                                        17 |                                                                 208 |                                      4424 |                                              93 |                                       112 |                                            5413 |                                       194 |                                       284 |                                           203 |                                          233 |                                          0 |                                             357 |                                        5529 |                                                 0 |                                          13 |                                           0 |                                     83 |                                      60 |                                         1 |                                           1 |                                       3 |                                      0 |                                          1 |                                           41 |                                          80 |                                         5 |                                       11160 |                                      395 |                                        93 |                                         1 |                                                        55 |                                             47 |                                      106 |                                      31 |                                  88 |                                         115 |
| ShotgunWGS-ControlPig5GutMicrobiome-Day7  |                14 |        3708 |            0 |           230 |        1991 |          64 |           876 |          787 |           277 |           29680 |         16 |             90 |          205 |               243 |            341 |          304 |       1189 |              445 |           896 |           1778 |         623 |           179 |        224 |            47 |      1684 |        27 |     31 |             383 |           494 |       30 |          19 |        1214 |         514 |         341 |          302 |            152 |             1170 |        647 |      3503 |             270 |        12270 |            205 |         0 |                0 |                7 |         334 |            632 |        105 |       815 |           168 |      448 |           139 |         3168 |         1475 |             1134 |         3828 |          8846 |        43 |          1142 |     438 |             523 |           438 |        460 |         416 |          735 |           5 |          88 |         0 |         0 |          87 |           125 |     12791 |           99 |           1 |            0 |           0 |      406 |          184 |          220 |         225 |       6 |    22326 |      225359 |        206 |    673 |            0 |                0 |                2 |          294 |        67 |          101 |         64 |               0 |                  0 |              0 |          218 |            1 |           66063 |            0 |            2 |             300 |              96 |   16872 |         2 |            0 |        830 |      262 |          24 |                  6 |             203 |        4201 |            669 |          1387 |            214 |           105 |      258 |     36 |        0 |      202 |      1097 |         2653 |        28995 |            112 |         0 |             5764 |                 8126 |              244 |         30 |           41 |          5906 |      36 |                       308 |                     231 |                       808 |                     57 |                     1 |                     83 |                    1241 |                      57 |                    0 |                     66 |                   569 |                      26 |                      72 |                    175 |                       177 |                          53 |                  27 |                 5 |                   922 |                31 |                   1 |           1881 |             3472 |              77 |            197 |            3799 |           193 |       178 |         546 |           0 |          191 |             2805 |        333 |          24 |         16 |           244 |         1743 |        86 |                   4 |           107 |           113 |         4 |               0 |           656 |       1789 |          706 |            438 |           9 |        0 |        0 |             398 |              364 |              297 |             0 |            230 |          31 |         470 |             29 |             0 |         310 |          8 |      272328 |           31 |               0 |            0 |          0 |       10233 |       357 |       194 |          297 |            150 |            0 |          19 |          1189 |       18751 |               299 |              302 |            1826 |      152 |          0 |          433 |          150 |         340 |            2680 |           0 |              42 |         957 |               6 |        16 |       600 |          4 |       1285 |                 41 |               0 |            0 |               0 |      1338 |           14 |           823 |           372 |            2490 |             226 |         825 |     329 |           550 |            0 |          35 |           279 |               608 |              13231 |              422 |           516 |             261 |              445 |                  56 |          673 |             8336 |          7391 |                273 |             208 |              24 |            478 |          521 |               940 |     12959 |           244 |     480 |          897 |           154 |             122 |      367 | 32973 |         1 |        1235 |          0 |          416 |        7772 |        62 |        48 |       0 |            450 |         45 |         3 |               4 |         37 |            50 |       220 |          505 |         9131 |             76 |             0 |                      6 |          539 |         130 |           39 |     297 |            223 |           296 |        1661 |          16498 |      181032 |       2 |            1686 |           114740 |        270 |          78 |          94 |              626 |        4040 |        672 |             49 |       2351 |           3360 |          0 |            1 |         530 |     890 |     0 |          36 |          6669 |         182 |               2 |        3518 |     173 |     102 |          211 |        5448 |      4373 |              164 |      36 |         76 |          0 |         372 |               133 |           153 |      115 |          0 |              0 |     1229 |           144 |            430 |         15 |        1322 |      0 |     462 |              35 |           747 |        277 |         86 |            82 |        42 |              55 |            39 |       132 |            51 |          80 |            211 |         40 |            69 |            1732 |              111 |             0 |         1893 |               1 |           4174 |          0 |            240 |           301 |           407 |           1 |      106 |         601 |      35 |      10907 |              0 |             68 |             100 |              146 |           41 |             64 |        188 |        0 |          0 |        475 |         42 |          24 |       1924 |        0 |              96 |          0 |        135 |        135 |               347 |     225 |         243 |       191 |                 74 |         314 |       61 |        633 |            23 |     130 |    154 |       353 |       168 |                1 |           142 |        126 |               2 |                  0 |        72 |       19 |        17 |       0 |        549401 |        1876 |                  81 |         1 |        188 |      455 |            975 |             1002 |        548 |       173 |         63 |          174 |      0 |         2 |        485 |              0 |        297 |         1940 |        1203 |          95 |     3963 |          0 |  13 |           24 |         51 |       97 |                 0 |                0 |      67 |            946 |         0 |         422 |          47 |           422 |              597 |         18 |           0 |        266 |        725 |        127 |          619 |            0 |         566 |            75 |             97 |       488 |       16036 |         481 |        172 |           330 |             50 |               1992 |                533 |          211 |              340 |          1187 |               1102 |            361 |              95 |               137 |           298 |          118 |           257 |          294 |           1574 |            339 |              143 |              371 |                 397 |                  25 |             81 |               130 |         370 |             302 |            73 |              615 |          141 |           399 |           60 |          105 |           196 |          147 |         13 |         120 |          68 |         256 |         57 |            216 |         317 |       50915 |       1108 |                0 |              7 |          0 |       56 |     3533 |       107 |       167 |               0 |              500 |          1488 |       1015 |        493 |                0 |               6 |        45 |         159 |           61 |            4 |           1383 |        55 |           75 |      151 |      18 |       918 |        0 |            33 |          95 |            0 |             90 |         49 |           296 |         316 |         116 |           503 |          465 |             31 |          295 |        193 |      207 |          413 |          162 |        62 |      0 |    613 |             290 |           20 |           76 |        111 |           143 |           1336 |          219 |          0 |             60 |         4 |          0 |        693 |          91 |      9197 |                 2 |      858 |         4884 |       34 |                0 |             0 |           73 |             88 |           0 |          112 |         878 |               2 |               3 |                2 |          7065 |         3905 |     494 |           22016 |            30 |             0 |                6 |        226 |          0 |         39 |           694 |          228 |          100 |         545 |         28 |        1 |            690 |        1178 |       2743 |         34 |       2142 |         804 |          3400 |          28 |          1528 |               1306 |        36 |           172 |       604 |          16 |            47 |            30 |          0 |              103 |                  3 |                   1 |                     0 |           1070 |          491 |          1 |           80 |     12 |         102 |       245 |          435 |           0 |         87 |            0 |           95 |        12 |          804 |         838 |            0 |              363 |               0 |       19 |          7825 |      6 |            2 |     507229 |             900 |               659 |              148 |            0 |     430 |          0 |         253 |                1 |               949 |        2830 |             1596 |           53 |           459 |           151 |          730 |           0 |         0 |           1687 |           0 |           9 |          99 |        427 |       1 |       704 |         0 |           21 |            2 |      211 |           114 |            3 |       650 |         672 |         617 |             82 |          3 |            592 |             1020 |            585 |          579 |            0 |        214 |            29 |        406 |           633 |     47030 |        34 |        1358 |         396 |         87 |         159 |    214 |       0 |         588 |         0 |      395 |       129723 |                4 |                39 |                   2 |               112 |            71 |            575 |               351 |            8 |        34 |          379 |         346 |        837 |          346 |           0 |       570 |           1 |              46 |            26 |                  59 |          15 |       2082 |           51 |       20968 |      575 |       3220 |      209 |           2773 |          178 |          68 |            0 |           606 |    9232 |     174 |       340 |           582 |              541 |         159 |          304 |          232 |        2377 |               3 |           0 |      1727 |            0 |            163 |           4225 |              54 |      110 |              266 |             2 |         139 |             645 |         41068 |         1825 |               269 |           43518 |            57 |        125 |           117 |                  408 |          296 |              260 |        286 |            2779 |           0 |          4149 |           393 |       1 |             770 |           2980 |              979 |        907 |               0 |              42 |               1 |               7 |          21 |          0 |            171 |         184 |               6 |          24 |            59 |           19 |     255 |        14 |             840 |               612 |        2151 |               5295 |                  2641 |           400 |          326 |           209 |          537 |           86 |                 290 |          82 |             250 |             185 |          217 |            11 |                 1748 |        2666 |         905 |             5 |                 334 |       1910 |     481 |              348 |          324 |            280 |       111 |       922 |         22 |      3432 |           335 |         467 |            4 |         27 |         59 |      165 |          46 |          104 |    14 |          815 |           13 |        160 |       37 |               0 |              18 |              0 |        322 |         1 |       13193 |               336 |              182 |            4 |   2870 |        1165 |     99 |           22 |      99 |       119 |           0 |             23 |       145 |       296 |          253 |         780 |         193 |         359 |     296 |       52 |            0 |     1077 |          964 |                17 |       192 |               6 |                   0 |                 12 |                                                133 |                                             989 |                                             1 |                                             124 |                                         113 |                                  840 |                                      1026 |                                             30 |                                         600 |                                            43 |                                                  22 |                                      102 |                                      0 |                                        26 |                                                                 369 |                                     10621 |                                             190 |                                       215 |                                           13467 |                                       387 |                                       493 |                                           425 |                                          349 |                                          0 |                                             607 |                                       15029 |                                                 0 |                                          32 |                                           0 |                                    173 |                                     190 |                                         0 |                                           1 |                                       6 |                                      0 |                                          3 |                                           74 |                                          95 |                                        10 |                                       24441 |                                      718 |                                       129 |                                         0 |                                                       197 |                                            139 |                                      122 |                                      41 |                                 200 |                                          13 |
| ShotgunWGS-TomatoPig18GutMicrobiome-Day7  |                 1 |        1159 |            0 |           146 |         585 |          33 |           265 |          338 |           195 |           10604 |          6 |             85 |          139 |               135 |            171 |          261 |       1161 |              145 |           368 |            634 |         527 |           172 |        100 |            49 |       356 |        12 |     13 |             115 |           267 |       17 |           9 |         717 |         512 |         104 |          187 |            138 |              488 |        184 |      1636 |             139 |         4443 |             90 |         0 |                0 |               92 |          74 |            196 |         34 |       225 |           148 |      229 |            53 |         1402 |          450 |              593 |         1437 |          2878 |        27 |           467 |     183 |             244 |           148 |        194 |         183 |          713 |           2 |          37 |         0 |         0 |          79 |           104 |     15164 |           47 |           0 |            0 |           0 |      197 |          110 |          105 |          84 |       5 |     8839 |      168421 |         67 |    214 |            0 |                0 |                0 |          146 |        26 |           73 |         20 |               0 |                  0 |              0 |          203 |            0 |          104397 |            0 |            0 |             176 |              45 |    7266 |         2 |            0 |        447 |      134 |          17 |                 13 |             166 |        2073 |            462 |           502 |            216 |           102 |      134 |     49 |        0 |       32 |       378 |         1349 |        10996 |             58 |         1 |             2222 |                 3243 |              109 |         14 |           11 |          1242 |      10 |                       129 |                      87 |                       418 |                     19 |                     1 |                     25 |                     429 |                       8 |                    0 |                     25 |                   352 |                       7 |                      43 |                     48 |                        69 |                          27 |                   4 |                 0 |                   773 |                10 |                   0 |           1136 |             1109 |              20 |             57 |            1120 |           178 |        58 |         374 |           0 |          201 |              751 |        223 |          11 |          4 |           170 |         1247 |        39 |                   2 |            64 |            51 |         0 |               0 |           333 |       1052 |          372 |            250 |          16 |        0 |        0 |             173 |              132 |              203 |             0 |            168 |          13 |         459 |             19 |             0 |         397 |         15 |      101208 |           12 |               0 |            0 |          0 |       15652 |       149 |       183 |          150 |             67 |            0 |          12 |           389 |        5587 |               101 |              172 |            1853 |       73 |          0 |          204 |           70 |         189 |            1767 |           0 |              10 |         594 |               1 |         6 |       704 |          5 |        614 |                 20 |               0 |            0 |               0 |       789 |           11 |           270 |           104 |            1041 |             230 |         443 |     248 |           232 |            0 |          37 |           106 |               209 |               4918 |              195 |           238 |             109 |              174 |                  33 |          299 |             2951 |          1806 |                 91 |             115 |              14 |            217 |          103 |               314 |      3230 |           101 |     198 |          423 |            62 |              80 |      220 |  8626 |         0 |         796 |          0 |          129 |        4429 |        30 |        20 |       0 |            155 |         28 |         0 |               1 |          6 |            14 |        31 |          352 |         4166 |             30 |             0 |                      0 |          170 |          41 |           28 |     160 |             64 |           214 |       16141 |           8612 |       66371 |       0 |             677 |           136627 |         60 |          35 |          31 |              245 |        2071 |        169 |             40 |        881 |           2888 |          0 |            0 |         223 |     766 |     0 |          18 |          2416 |          66 |               5 |        4843 |      44 |      91 |          208 |        2135 |      2027 |              160 |       4 |         70 |          0 |         231 |                73 |           106 |       92 |          1 |              0 |      678 |           134 |            149 |          6 |         431 |      0 |     169 |              14 |           224 |        144 |         54 |            20 |        14 |              20 |            21 |        60 |            26 |          22 |            143 |         14 |            39 |             643 |               68 |             0 |          927 |               0 |           1360 |          0 |            173 |           142 |           244 |           0 |       40 |         183 |      15 |       4918 |              0 |             19 |              26 |               71 |           14 |             84 |         99 |        0 |          0 |        122 |          9 |           6 |        638 |        1 |             130 |          0 |        167 |         88 |               197 |     238 |         103 |        54 |                 44 |         326 |        9 |        543 |            23 |     114 |     86 |       142 |       175 |                0 |            62 |        107 |               2 |                  1 |        45 |        5 |        16 |       0 |        185934 |         932 |                 255 |         0 |         74 |      181 |            612 |              835 |        237 |       191 |         31 |           51 |      0 |         0 |        235 |              2 |        335 |          796 |         670 |          72 |     1640 |          2 |  27 |           11 |         42 |       30 |                 0 |                0 |      29 |            367 |         0 |         131 |          13 |           201 |              269 |         10 |           0 |        122 |        443 |         82 |          224 |            0 |         208 |            23 |             46 |       242 |       11524 |         250 |         66 |           264 |             20 |                691 |                176 |           70 |              159 |           427 |                485 |            171 |              39 |                37 |           115 |           56 |           122 |           95 |            508 |            171 |               69 |              141 |                 147 |                  19 |             18 |                72 |         412 |             129 |            50 |              441 |          129 |           222 |           15 |          125 |            72 |           71 |         10 |          89 |          41 |         145 |         26 |            151 |         154 |        8831 |       1015 |                0 |              6 |          0 |       24 |     1204 |        29 |        34 |               0 |              318 |          1337 |        396 |        327 |                0 |               0 |        28 |         158 |           38 |            1 |            508 |        19 |           14 |       42 |       9 |       317 |        0 |            14 |          57 |            2 |             19 |         17 |           112 |         209 |          59 |           229 |          261 |             15 |          112 |         81 |      177 |          491 |          122 |        23 |      1 |    347 |             238 |           13 |           54 |         83 |            79 |            493 |          138 |          0 |             22 |         0 |          0 |        344 |          62 |     11037 |                 1 |      718 |         1696 |       13 |                0 |             0 |           39 |             36 |           0 |           52 |         154 |             120 |               2 |               15 |          2815 |         3156 |     176 |           15566 |            15 |             0 |                1 |        111 |          0 |         23 |           886 |          143 |           54 |         114 |         17 |        0 |            246 |         775 |       1779 |         30 |        906 |         503 |          1215 |          13 |           452 |                381 |        25 |            75 |       249 |          12 |            36 |            17 |          0 |               71 |                  2 |                   1 |                     0 |            288 |          155 |          0 |           25 |      7 |          39 |       168 |          303 |           2 |         73 |            0 |           70 |        11 |          515 |         947 |            0 |              332 |               0 |        6 |          4741 |      4 |            0 |     478207 |             723 |               647 |               60 |            0 |     133 |          0 |          60 |                2 |               334 |        1457 |              585 |           15 |           177 |            93 |          177 |           0 |         0 |            559 |           0 |           7 |          42 |        179 |       0 |       380 |         0 |           11 |            2 |       73 |           111 |            1 |       364 |         457 |         492 |             57 |          0 |            358 |              537 |            272 |          315 |            0 |        128 |             4 |        183 |           365 |     15783 |        17 |         734 |         267 |         67 |          92 |    154 |       0 |         361 |         0 |      216 |        47208 |                1 |                 9 |                   7 |                82 |           116 |            419 |               362 |            2 |        19 |          197 |         330 |       1070 |          303 |           0 |       573 |           1 |              24 |             4 |                  40 |           8 |        860 |           23 |        2400 |      247 |       1106 |     2014 |            941 |           84 |          64 |            0 |           354 |    4993 |      46 |       192 |           395 |              264 |         127 |          230 |          169 |         841 |               0 |           0 |      1487 |            0 |            101 |           3009 |              23 |       79 |              145 |             0 |         113 |             358 |         10184 |         1565 |               216 |           38373 |            45 |         53 |            43 |                  253 |          100 |               70 |         75 |            1110 |           0 |          3488 |           206 |       0 |             304 |           1271 |              297 |        318 |               0 |              37 |               0 |               5 |           4 |          0 |             88 |          81 |               4 |          10 |            28 |            6 |     115 |         3 |             343 |               223 |         622 |               2078 |                   913 |           198 |          302 |           160 |          171 |           31 |                 162 |          41 |             141 |             151 |           88 |             4 |                  577 |         674 |         312 |            10 |                 189 |        714 |     219 |              162 |          151 |             84 |        89 |       130 |         10 |      1445 |           156 |         186 |            0 |         29 |         38 |      111 |          24 |           83 |    10 |          210 |            3 |         67 |       38 |               0 |              26 |              0 |        309 |         0 |        5017 |               314 |              151 |            3 |    962 |         237 |     72 |            6 |      48 |        68 |           0 |             11 |        64 |       125 |          154 |         485 |          67 |         328 |     106 |       18 |            0 |      547 |          736 |                12 |       108 |               6 |                   0 |                  5 |                                                139 |                                             314 |                                             0 |                                              53 |                                          10 |                                 1208 |                                       391 |                                              5 |                                         119 |                                            15 |                                                  17 |                                       41 |                                      0 |                                        12 |                                                                 126 |                                      4323 |                                              84 |                                        75 |                                            3971 |                                       141 |                                       272 |                                           292 |                                          164 |                                          0 |                                             273 |                                        5333 |                                                 0 |                                           6 |                                           0 |                                    105 |                                     100 |                                         0 |                                           1 |                                      19 |                                      0 |                                          2 |                                           32 |                                          25 |                                         8 |                                       12458 |                                      413 |                                        48 |                                         0 |                                                       182 |                                             81 |                                       12 |                                       9 |                                  87 |                                          59 |

Calculate relative abundance, and bind back to metadata.

``` r
GenusOnly.Counts.Filt.t.wtotal <- GenusOnly.Counts.Filt.t %>%
  mutate(Total.Counts = rowSums(GenusOnly.Counts.Filt.t[,2:ncol(GenusOnly.Counts.Filt.t)]))

dim(GenusOnly.Counts.Filt.t.wtotal)
```

    ## [1]  60 897

``` r
# create rel abund df
RelAbund.Genus.Filt <- GenusOnly.Counts.Filt.t.wtotal[,2:896]/GenusOnly.Counts.Filt.t.wtotal$Total.Counts

# add back metadata
RelAbund.Genus.Filt <- bind_cols(AllSamples.Metadata, RelAbund.Genus.Filt)
```

### Counting missing data

The goal of these next bits of code are to understand how many missing
values we have in our dataset, to set what parameters we will use for
filtering.

``` r
# how many zeros are in the column AHJD-like viruses?
sum(RelAbund.Genus.Filt$`AHJD-like viruses` == 0)  # this code works
```

    ## [1] 31

``` r
# remove metadata   
# metadata is all character or factor, so can select only numeric columns
RelAbund.Genus.Filt.nometadata <- RelAbund.Genus.Filt %>%
  select_if(is.numeric)

# create a list with the number of zeros for each genus
counting_zeros <- sapply(RelAbund.Genus.Filt.nometadata, 
                         function(x){ (sum(x==0))})

# plot a histogram to look at missing values
counting_zeros_df <- as.data.frame(counting_zeros)

hist(counting_zeros_df$counting_zeros, 
     breaks = 61,
     main = "Histogram of Genera with Zero Relative Intensity",
     sub = "Starting at No Zeros",
     xlab = "Number of zero relative intensity values",
     ylab = "Frequency")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

First column is no missing values, and its so big its hard to see how
many missing values we actually have.

``` r
# plot a histogram to look, but removing genera that are only missing 1 value
counting_zeros_df_missingval <- counting_zeros_df %>%
  rownames_to_column(var = "rowname") %>%
  filter(counting_zeros > 0) %>%
  column_to_rownames(var = "rowname")

# how many genera have at least one missing value?
dim(counting_zeros_df_missingval)
```

    ## [1] 186   1

186 genera have at least one missing value.

``` r
# histogram of number of zeros, starting at 1 zero
hist(counting_zeros_df_missingval$counting_zeros, 
     breaks = 60,
     main = "Histogram of Genera with Zero Relative Intensity",
     sub = "Starting at 1 Zero",
     xlab = "Number of zero relative intensity values",
     ylab = "Frequency")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# plot a histogram to look, but removing genera that have 20 or more zeros
counting_zeros_df_missing20ormore <- counting_zeros_df %>%
  rownames_to_column(var = "rowname") %>%
  filter(counting_zeros >= 20) %>%
  column_to_rownames(var = "rowname")

# histogram of number of zeros, starting at 20 zero
hist(counting_zeros_df_missing20ormore$counting_zeros, 
     breaks = 40,
     main = "Histogram of Genera with Zero Relative Intensity",
     sub = "Starting at 20 Zero",
     xlab = "Number of zero relative intensity values",
     ylab = "Frequency")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
# how many genera have 20 or more missing value?
dim(counting_zeros_df_missing20ormore)
```

    ## [1] 140   1

There are 140 genera that have 20 or more missing values. Because 20
missing values here is 1/3 missing, we decided to use this as our
cutoff.

### Filter for &lt;33% missingness

Our decided criteria:  
Filter out genera from relative abundance table that have &gt; 20 zeros,
or more than 33% missing data.

``` r
# make a character vector of genera names that have > 20 zeros from the rownames in above table
zeros.20 <- c(rownames(counting_zeros_df_missing20ormore))

# filter using this list
RelAbund.Genus.Filt.zerofilt <- RelAbund.Genus.Filt %>%
  rownames_to_column(var = "rowname") %>%
  select(everything(), -all_of(zeros.20)) %>%
  column_to_rownames(var = "rowname")

RelAbund.Genus.Filt.zerofilt[1:3,1:6]
```

    ##                                 Sample_Name Pig    Diet Time_Point
    ## 1 ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## 2  ShotgunWGS-ControlPig8GutMicrobiome-Day0   8 Control      Day 0
    ## 3 ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ##   Diet_By_Time_Point Abiotrophia
    ## 1     Control Day 14 0.001305713
    ## 2      Control Day 0 0.001347804
    ## 3     Control Day 14 0.001066255

``` r
dim(RelAbund.Genus.Filt.zerofilt)
```

    ## [1]  60 760

Our final dataset has 755 genera (because there are 5 columns of
metadata).

Write final dataset genus rel abund to .csv this way we have it.

``` r
write_csv(RelAbund.Genus.Filt.zerofilt,
          file = "Genus_RelAbund_Final_Filtered_WithMetadata.csv")
```

## Microbiome profile

Wrangling to enable collection of some summary statistics about our
microbiome profile, including how many genera belong to different
domains, etc.

### Wrangling

Grab names of final genera.

``` r
# contains inplausible genera removed, but not removed for zeroes
dim(Genus.Counts.Filt)
```

    ## [1] 895  66

``` r
Genus.Counts.Filt[1:5, 1:10]
```

    ## # A tibble: 5 × 10
    ##   domain    phylum    class order family genus `ShotgunWGS-Co…` `ShotgunWGS-Co…`
    ##   <chr>     <chr>     <chr> <chr> <chr>  <chr>            <dbl>            <dbl>
    ## 1 Viruses   unclassi… uncl… Caud… Podov… AHJD…               29                0
    ## 2 Bacteria  Firmicut… Baci… Lact… Aeroc… Abio…             5067             5661
    ## 3 Eukaryota unclassi… uncl… uncl… uncla… Acan…                0                0
    ## 4 Bacteria  Cyanobac… uncl… uncl… uncla… Acar…              271              416
    ## 5 Bacteria  Firmicut… Clos… Clos… Rumin… Acet…             1988             2981
    ## # … with 2 more variables: `ShotgunWGS-ControlPig3GutMicrobiome-Day14` <dbl>,
    ## #   `ShotgunWGS-TomatoPig14GutMicrobiome-Day7` <dbl>

``` r
# final filtered data
RelAbund.Genus.Filt.zerofilt[1:5, 1:10]
```

    ##                                 Sample_Name Pig    Diet Time_Point
    ## 1 ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## 2  ShotgunWGS-ControlPig8GutMicrobiome-Day0   8 Control      Day 0
    ## 3 ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ## 4  ShotgunWGS-TomatoPig14GutMicrobiome-Day7  14  Tomato      Day 7
    ## 5  ShotgunWGS-ControlPig5GutMicrobiome-Day7   5 Control      Day 7
    ##   Diet_By_Time_Point Abiotrophia Acaryochloris  Acetivibrio  Acetobacter
    ## 1     Control Day 14 0.001305713  6.983388e-05 0.0005122869 1.700751e-05
    ## 2      Control Day 0 0.001347804  9.904370e-05 0.0007097339 2.047538e-05
    ## 3     Control Day 14 0.001066255  6.914992e-05 0.0005363651 1.553931e-05
    ## 4       Tomato Day 7 0.001311580  1.090209e-04 0.0008422076 3.412106e-05
    ## 5      Control Day 7 0.001207244  7.488298e-05 0.0006482261 2.083700e-05
    ##   Acetohalobium
    ## 1  0.0002669664
    ## 2  0.0003268918
    ## 3  0.0002628733
    ## 4  0.0003320562
    ## 5  0.0002852065

``` r
# grab colnames which have all the final genera
final_genera <- colnames(RelAbund.Genus.Filt.zerofilt)

# remove metadata colnames
final_genera <- final_genera[6:760]  

final_genera <- as.data.frame(final_genera)

# create a df with the final genera we want to keep for our analysis
final_genera <- final_genera %>%
  rename(genus = final_genera)
```

Get back domain and `inner_join()` with `final_genera` list

``` r
# pull from full dataset the domain and genus columns
Genus.Counts.Filt.Domain.Genera <- Genus.Counts.Filt %>%
  select(domain, genus)

Genus.Counts.Filt.Domain.Genera[1:10,]
```

    ## # A tibble: 10 × 2
    ##    domain    genus            
    ##    <chr>     <chr>            
    ##  1 Viruses   AHJD-like viruses
    ##  2 Bacteria  Abiotrophia      
    ##  3 Eukaryota Acanthamoeba     
    ##  4 Bacteria  Acaryochloris    
    ##  5 Bacteria  Acetivibrio      
    ##  6 Bacteria  Acetobacter      
    ##  7 Bacteria  Acetohalobium    
    ##  8 Bacteria  Acholeplasma     
    ##  9 Bacteria  Achromobacter    
    ## 10 Bacteria  Acidaminococcus

``` r
# want to join Genus.Counts.Filt.Domain.Genera with final_genera
final_genera_withdomain <- inner_join(final_genera, Genus.Counts.Filt.Domain.Genera,
                                      by = "genus")
```

### Count genera

``` r
final_genera_withdomain %>%
  count()
```

    ##     n
    ## 1 755

``` r
final_genera_withdomain %>%
  group_by(domain) %>%
  count()
```

    ## # A tibble: 5 × 2
    ## # Groups:   domain [5]
    ##   domain              n
    ##   <chr>           <int>
    ## 1 Archaea            60
    ## 2 Bacteria          582
    ## 3 Eukaryota          89
    ## 4 other sequences     1
    ## 5 Viruses            23

We have 755 total genera. We have 60 genera from Archaea, 582 from
Bacteria, 89 from Eukaryota, and 23 from Viruses.

### Most prevalent genera

What are the most prevalent genera in our pigs?

``` r
RelAbund.Genus.Filt.zerofilt[1:5, 1:10]
```

    ##                                 Sample_Name Pig    Diet Time_Point
    ## 1 ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## 2  ShotgunWGS-ControlPig8GutMicrobiome-Day0   8 Control      Day 0
    ## 3 ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ## 4  ShotgunWGS-TomatoPig14GutMicrobiome-Day7  14  Tomato      Day 7
    ## 5  ShotgunWGS-ControlPig5GutMicrobiome-Day7   5 Control      Day 7
    ##   Diet_By_Time_Point Abiotrophia Acaryochloris  Acetivibrio  Acetobacter
    ## 1     Control Day 14 0.001305713  6.983388e-05 0.0005122869 1.700751e-05
    ## 2      Control Day 0 0.001347804  9.904370e-05 0.0007097339 2.047538e-05
    ## 3     Control Day 14 0.001066255  6.914992e-05 0.0005363651 1.553931e-05
    ## 4       Tomato Day 7 0.001311580  1.090209e-04 0.0008422076 3.412106e-05
    ## 5      Control Day 7 0.001207244  7.488298e-05 0.0006482261 2.083700e-05
    ##   Acetohalobium
    ## 1  0.0002669664
    ## 2  0.0003268918
    ## 3  0.0002628733
    ## 4  0.0003320562
    ## 5  0.0002852065

``` r
genera_means <- RelAbund.Genus.Filt.zerofilt %>%
  summarize_if(is.numeric, mean)

genera_means_t <- t(genera_means)
genera_means_t <- as.data.frame(genera_means_t)

genera_means_t <- genera_means_t %>%
  rename(rel_abund_genera = V1) %>%
  arrange(-rel_abund_genera)

head(genera_means_t)
```

    ##                  rel_abund_genera
    ## Prevotella             0.22231328
    ## Bacteroides            0.10347888
    ## Clostridium            0.08556113
    ## Lactobacillus          0.06777787
    ## Eubacterium            0.05164571
    ## Faecalibacterium       0.04480044

The most prevalent genera are Prevotella (22.23% average abundance),
Bacteroides (10.34%), Clostridium (8.56%), Lactobacillus (6.78%) and
Eubacterium (5.16%).

## Rarefaction curves

### Create tax and OTU tables

This section uses a different package than the rest of the analysis;
data and metadata need to be uploaded again and made into format
friendly for package.

``` r
# tax table
TAX_tab <- Genus.AllSamples.Counts %>% 
  select(Domain = domain, Phylum = phylum,
         Class = class, Order = order,
         Family = family, Genus = genus)

tax_names <- colnames(TAX_tab)

head(TAX_tab)
```

    ## # A tibble: 6 × 6
    ##   Domain    Phylum                                Class       Order Family Genus
    ##   <chr>     <chr>                                 <chr>       <chr> <chr>  <chr>
    ## 1 Viruses   unclassified (derived from Viruses)   unclassifi… Caud… Podov… AHJD…
    ## 2 Bacteria  Firmicutes                            Bacilli     Lact… Aeroc… Abio…
    ## 3 Eukaryota unclassified (derived from Eukaryota) unclassifi… uncl… uncla… Acan…
    ## 4 Bacteria  Cyanobacteria                         unclassifi… uncl… uncla… Acar…
    ## 5 Bacteria  Firmicutes                            Clostridia  Clos… Rumin… Acet…
    ## 6 Bacteria  Proteobacteria                        Alphaprote… Rhod… Aceto… Acet…

``` r
head(tax_names)
```

    ## [1] "Domain" "Phylum" "Class"  "Order"  "Family" "Genus"

``` r
#  OTU table
OTU_tab <- Genus.AllSamples.Counts[, seq(7, 66)]

head(OTU_tab)
```

    ## # A tibble: 6 × 60
    ##   `ShotgunWGS-ControlPig6Gu…` `ShotgunWGS-Co…` `ShotgunWGS-Co…` `ShotgunWGS-To…`
    ##                         <dbl>            <dbl>            <dbl>            <dbl>
    ## 1                          29                0              153                0
    ## 2                        5067             5661             4117             1576
    ## 3                           0                0                0                1
    ## 4                         271              416              267              131
    ## 5                        1988             2981             2071             1012
    ## 6                          66               86               60               41
    ## # … with 56 more variables: `ShotgunWGS-ControlPig5GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-TomatoPig18GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-TomatoPig16GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-ControlPig10GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-ControlPig2GutMicrobiome-Day0` <dbl>,
    ## #   `ShotgunWGS-TomatoPig18GutMicrobiome-Day0` <dbl>,
    ## #   `ShotgunWGS-ControlPig10GutMicrobiome-Day0` <dbl>, …

### Create metadata

Since metadata is contained in column names, we will parse them from
here.

``` r
raw_names <- colnames(OTU_tab)
names_table <- data.frame(Raw_names = raw_names)
```

First, the string will split by the middle hyphen.

``` r
names_table <- names_table %>% 
  separate(Raw_names, into = c("Shotgun", "Type", "Day")) %>% 
  select(-Shotgun)
```

Now, since the character `GutMicrobiome` is constant over all samples,
it will be removed. In the same manner, the character `Day` will be
removed.

``` r
names_table <- names_table %>% 
  mutate(Type = str_remove(string = Type, pattern = "GutMicrobiome") ) %>% 
  mutate(Type = str_remove(string = Type, pattern = "GutMicrobime") ) %>% 
  mutate(Day = str_remove(string = Day, pattern = "Day"))
head(names_table, 2)
```

    ##          Type Day
    ## 1 ControlPig6  14
    ## 2 ControlPig8   0

Since `Pig` is in the middle of the sample type and the pig number, it
will be used as separator character. And the final result is a tidy
data.

``` r
names_table <- names_table %>% 
  separate(col = Type, into = c("Type", "Pig"), sep = "Pig") %>% 
  mutate(Type = factor(Type), Pig = factor(Pig), Day = as.integer(Day)) %>% 
  select(Type, Day, Pig)
head(names_table)
```

    ##      Type Day Pig
    ## 1 Control  14   6
    ## 2 Control   0   8
    ## 3 Control  14   3
    ## 4  Tomato   7  14
    ## 5 Control   7   5
    ## 6  Tomato   7  18

### Renaming samples

Now that we have tidy data, it’s better to replace long names with
shorter ones. New names will be created as *Type\_Pig\_Day*. We are also
creating names that distinguish the 6 diet by time point groups.

``` r
tmp_names  <- names_table %>%
  mutate(Pig = paste0("P", Pig), Day = paste0("D", Day)) %>% 
  unite("Kronas", Type:Day) %>% unite("Sample", Kronas:Pig, remove = F) %>% 
  select(-Pig)
head(tmp_names)
```

    ##           Sample      Kronas
    ## 1 Control_D14_P6 Control_D14
    ## 2  Control_D0_P8  Control_D0
    ## 3 Control_D14_P3 Control_D14
    ## 4  Tomato_D7_P14   Tomato_D7
    ## 5  Control_D7_P5  Control_D7
    ## 6  Tomato_D7_P18   Tomato_D7

``` r
metadata <- bind_cols(names_table, tmp_names)
head(metadata)
```

    ##      Type Day Pig         Sample      Kronas
    ## 1 Control  14   6 Control_D14_P6 Control_D14
    ## 2 Control   0   8  Control_D0_P8  Control_D0
    ## 3 Control  14   3 Control_D14_P3 Control_D14
    ## 4  Tomato   7  14  Tomato_D7_P14   Tomato_D7
    ## 5 Control   7   5  Control_D7_P5  Control_D7
    ## 6  Tomato   7  18  Tomato_D7_P18   Tomato_D7

Finally, in order to use this metadata with the OTU table, which is
linked by `names`, **the row names of the metadata and the column names
in the OTU table must be the same**.

``` r
rownames(names_table) <-  tmp_names$Sample
colnames(OTU_tab) <-   tmp_names$Sample
```

### Creating a phyloseq object

It’s time create a *phyloseq* object that will allow us to analyze this
data easier.

``` r
rownames(metadata) <- metadata$Sample
gut_microbiome_raw <- phyloseq(otu_table(OTU_tab, taxa_are_rows = T),
                               tax_table(TAX_tab),
                               sample_data(metadata))

colnames(tax_table(gut_microbiome_raw) ) <- tax_names
```

We can see that data at genus level accounts with 1085 taxas in 60
samples. But, we had developed a filtering scheme to remove very low
abundance and inconsistently detected taxa, so let’s merge this full
list

``` r
gut_microbiome_raw
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1085 taxa and 60 samples ]
    ## sample_data() Sample Data:       [ 60 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1085 taxa by 6 taxonomic ranks ]

### Filtering taxa

We want to only include the taxa we ended up using in our final
analysis. We have already created an df `final_genera` which contains
only the genera used in our final analysis.

``` r
final_genera_forphyloseq <- final_genera$genus

# subset to include only final genera
gut_microbiome_clean <- subset_taxa(gut_microbiome_raw, Genus %in% final_genera_forphyloseq)

gut_microbiome_clean
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 755 taxa and 60 samples ]
    ## sample_data() Sample Data:       [ 60 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 755 taxa by 6 taxonomic ranks ]

The final phyloseq object has 871 taxonomic levels, in our cases,
species since it is the lowest taxonomic levels that the sequences were
annotated.

In order to check if we have the 45 phyla, we are gonna to count the
unique phyla names in the dataset.

``` r
length(unique(tax_table(gut_microbiome_clean)[, "Genus"]))
```

    ## [1] 755

We got our 755 genera.

### Creating rarefaction curves

``` r
plot_rarefaction <- ranacapa::ggrare(gut_microbiome_clean, step = 60000,
                                   color = 'Type',  se = F, plot = F) 
```

    ## rarefying sample Control_D14_P6
    ## rarefying sample Control_D0_P8
    ## rarefying sample Control_D14_P3
    ## rarefying sample Tomato_D7_P14
    ## rarefying sample Control_D7_P5
    ## rarefying sample Tomato_D7_P18
    ## rarefying sample Tomato_D7_P16
    ## rarefying sample Control_D7_P10
    ## rarefying sample Control_D0_P2
    ## rarefying sample Tomato_D0_P18
    ## rarefying sample Control_D0_P10
    ## rarefying sample Control_D0_P7
    ## rarefying sample Control_D14_P8
    ## rarefying sample Tomato_D0_P11
    ## rarefying sample Tomato_D0_P19
    ## rarefying sample Tomato_D14_P17
    ## rarefying sample Control_D14_P9
    ## rarefying sample Control_D14_P10
    ## rarefying sample Tomato_D7_P19
    ## rarefying sample Control_D14_P5
    ## rarefying sample Control_D7_P2
    ## rarefying sample Control_D7_P6
    ## rarefying sample Tomato_D0_P12
    ## rarefying sample Tomato_D0_P14
    ## rarefying sample Control_D14_P7
    ## rarefying sample Tomato_D14_P11
    ## rarefying sample Tomato_D0_P20
    ## rarefying sample Control_D0_P9
    ## rarefying sample Tomato_D7_P11
    ## rarefying sample Tomato_D7_P13
    ## rarefying sample Tomato_D0_P17
    ## rarefying sample Tomato_D14_P19
    ## rarefying sample Tomato_D0_P13
    ## rarefying sample Control_D14_P2
    ## rarefying sample Control_D7_P1
    ## rarefying sample Tomato_D7_P15
    ## rarefying sample Tomato_D0_P15
    ## rarefying sample Tomato_D7_P12
    ## rarefying sample Tomato_D14_P14
    ## rarefying sample Tomato_D14_P20
    ## rarefying sample Control_D0_P1
    ## rarefying sample Control_D14_P4
    ## rarefying sample Control_D0_P6
    ## rarefying sample Tomato_D0_P16
    ## rarefying sample Tomato_D14_P16
    ## rarefying sample Tomato_D14_P18
    ## rarefying sample Control_D7_P7
    ## rarefying sample Control_D7_P4
    ## rarefying sample Tomato_D14_P13
    ## rarefying sample Control_D7_P8
    ## rarefying sample Tomato_D14_P15
    ## rarefying sample Tomato_D14_P12
    ## rarefying sample Tomato_D7_P20
    ## rarefying sample Control_D14_P1
    ## rarefying sample Control_D0_P3
    ## rarefying sample Control_D0_P5
    ## rarefying sample Control_D0_P4
    ## rarefying sample Control_D7_P9
    ## rarefying sample Control_D7_P3
    ## rarefying sample Tomato_D7_P17

``` r
plot_rarefaction <-  plot_rarefaction + theme_test() +
  facet_wrap("Day", scales = "free_x", ncol = 1, 
             labeller = labeller(Day = c( `0` = "Day 0" ,
                                          `7` = "Day 7",
                                          `14`= "Day 14")) ) + 
  labs(color = "Diet",
       title = "Rarefaction curves") +
  scale_color_manual(values = c( "steelblue2", "tomato2"))
  
plot_rarefaction
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

## Krona plots for exploratory analysis

The `psadd` package is able to create Krona plots with an phyloseq
object. Two Kronas will be created, per sample and per category
`Day + Type`. The Krona plots only include the final filtered taxa we
used in our analysis.

``` r
# Write kronas per samples
plot_krona(physeq = gut_microbiome, output = "./plots/kronas/per_sample", 
           variable = "Sample")

# Write kronas per category (Sample type + Day) i.e. Tomato_D7
plot_krona(physeq = gut_microbiome, output = "./plots/kronas/per_category/", 
           variable = "Kronas")
```

## PERMANOVA

Use PERMANOVA to conduct statistical analysis of overall microbial
profile differences among groups.

### All samples, full model

Test the overall effect of `Diet`, `Time_Point` and their interaction of
the overall microbiome. ORIGINAL

``` r
# create factors
factors_time_diet_pig <- RelAbund.Genus.Filt.zerofilt %>% select(Time_Point, Diet, Pig)

# create permutations
perm_time_diet_pig <- how(nperm = 9999)
setBlocks(perm_time_diet_pig) <- with(factors_time_diet_pig, Pig)

# run permanova
AllData.Genus.Filt.permanova <- adonis2(RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Diet*Time_Point,
                                        data = factors_time_diet_pig,
                                        permutations = perm_time_diet_pig,
                                        method = "bray")

AllData.Genus.Filt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  with(factors_time_diet_pig, Pig) 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Diet * Time_Point, data = factors_time_diet_pig, permutations = perm_time_diet_pig, method = "bray")
    ##                 Df SumOfSqs      R2      F Pr(>F)    
    ## Diet             1  0.05879 0.04954 3.4114 0.0001 ***
    ## Time_Point       2  0.16612 0.13999 4.8196 0.0001 ***
    ## Diet:Time_Point  2  0.03113 0.02623 0.9031 0.3831    
    ## Residual        54  0.93061 0.78424                  
    ## Total           59  1.18665 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

-   Diet: p = 0.0001, significant  
-   Time\_Point: p = 0.0001, significant  
-   Diet\*Time\_Point: p = 0.3831, non-significant

Comparison when you don’t filter out for missing values

``` r
AllData.Genus.Filt.permanova.no0filt <- adonis2(RelAbund.Genus.Filt[,-c(1:5)]~Diet*Time_Point,
                                        data = factors_time_diet_pig,
                                        permutations = perm_time_diet_pig,
                                        method = "bray")
AllData.Genus.Filt.permanova.no0filt
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  with(factors_time_diet_pig, Pig) 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = RelAbund.Genus.Filt[, -c(1:5)] ~ Diet * Time_Point, data = factors_time_diet_pig, permutations = perm_time_diet_pig, method = "bray")
    ##                 Df SumOfSqs      R2      F Pr(>F)    
    ## Diet             1  0.05880 0.04954 3.4114 0.0001 ***
    ## Time_Point       2  0.16616 0.14000 4.8200 0.0001 ***
    ## Diet:Time_Point  2  0.03113 0.02623 0.9031 0.3750    
    ## Residual        54  0.93079 0.78423                  
    ## Total           59  1.18689 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

-   Diet: p = 0.0001, significant  
-   Time\_Point: p = 0.0001, significant  
-   Diet\*Time\_Point: p = 0.3750, non-significant

Significance is the same whether you filter for missing data or not.

NEW

``` r
set.seed(2021)
# create factors
factors_time_diet_pig_genus <- RelAbund.Genus.Filt.zerofilt %>% select(Time_Point, Diet, Pig)

# create permutations
perm_time_diet_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                plots = Plots(strata=factors_time_diet_pig_genus$Pig,
                                              type="free",))
# run permanova
AllData.Genus.Filt.permanova <- adonis2(RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Diet*Time_Point,
                                        data = factors_time_diet_pig_genus,
                                        permutations = perm_time_diet_pig_genus,
                                        method = "bray",
                                        by = "margin")

AllData.Genus.Filt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Plots: factors_time_diet_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Diet * Time_Point, data = factors_time_diet_pig_genus, permutations = perm_time_diet_pig_genus, method = "bray", by = "margin")
    ##                 Df SumOfSqs      R2      F Pr(>F)
    ## Diet:Time_Point  2  0.03113 0.02623 0.9031   0.36
    ## Residual        54  0.93061 0.78424              
    ## Total           59  1.18665 1.00000

Interaction not significant (p=.355) so remove from model

``` r
set.seed(2021)
# create factors
factors_time_diet_pig_genus <- RelAbund.Genus.Filt.zerofilt %>% select(Time_Point, Diet, Pig)

# create permutations
perm_time_diet_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                plots = Plots(strata=factors_time_diet_pig_genus$Pig, type="free",))
# run permanova
AllData.Genus.Filt.permanova <- adonis2(RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Diet+Time_Point,
                                        data = factors_time_diet_pig_genus,
                                        permutations = perm_time_diet_pig_genus,
                                        method = "bray",
                                        by = "margin")

AllData.Genus.Filt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Plots: factors_time_diet_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Diet + Time_Point, data = factors_time_diet_pig_genus, permutations = perm_time_diet_pig_genus, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## Diet        1  0.05879 0.04954 3.4232  0.060 . 
    ## Time_Point  2  0.16612 0.13999 4.8364  0.005 **
    ## Residual   56  0.96174 0.81047                 
    ## Total      59  1.18665 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Diet not significant p=.060 but close Time significant p=.005

Test for homogeneity of multivariate dispersions

``` r
dis <- vegdist(RelAbund.Genus.Filt.zerofilt[,-c(1:5)], method = "bray")
mod <- betadisper(dis, RelAbund.Genus.Filt.zerofilt$Diet)
permutest(mod)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.001513 0.0015128 1.0653    999  0.332
    ## Residuals 58 0.082364 0.0014201

``` r
dis <- vegdist(RelAbund.Genus.Filt.zerofilt[,-c(1:5)], method = "bray")
mod <- betadisper(dis, RelAbund.Genus.Filt.zerofilt$Time)
permutest(mod)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups     2 0.008867 0.0044335 2.6538    999  0.081 .
    ## Residuals 57 0.095227 0.0016707                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

MANOVA TRIAL

``` r
a <- do.call(rbind, lapply(RelAbund.Genus.Filt.zerofilt, as.data.frame))

dep_vars <- cbind(RelAbund.Genus.Filt.zerofilt[-c(1:5)])
fit <- manova(cbind(RelAbund.Genus.Filt.zerofilt$Abiotrophia,RelAbund.Genus.Filt.zerofilt$Acaryochloris)~Diet*Time_Point + (1|Pig), data=RelAbund.Genus.Filt.zerofilt)

tidy(fit)
```

### Post Hoc PERMANOVA within Each Diet

#### Within Control Diet Only

Effect in control diet of time.

``` r
set.seed(2021)
# filter for only control
control.RelAbund.Genus.Filt.zerofilt <- subset(RelAbund.Genus.Filt.zerofilt, Diet == "Control")

# create factors
factors_control_genera <- droplevels(control.RelAbund.Genus.Filt.zerofilt %>% select(Time_Point, Pig))

# create permutations
perm_control_genera <- how(within = Within(type="series", constant=TRUE),
                                plots = Plots(strata=factors_control_genera$Pig, type="none",))
# run permanova
Control.ByTime.Genus.zerofilt <- adonis2(control.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                        data = factors_control_genera,
                                        permutations = perm_control_genera,
                                        method = "bray",
                                        by = "margin")
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

``` r
Control.ByTime.Genus.zerofilt
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_genera$Pig, plot permutation: none
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 2
    ## 
    ## adonis2(formula = control.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_genera, permutations = perm_control_genera, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2     F Pr(>F)
    ## Time_Point  2  0.13507 0.22578 3.937 0.3333
    ## Residual   27  0.46317 0.77422             
    ## Total      29  0.59824 1.00000

Significant effect of time (p = 0.005) within control samples.

Now do pairwise comparisons to see where the significance is coming from

##### Control T1 vs Control T2

``` r
set.seed(2021)
# filter data set for only samples at T1 and T2
control.T1T2.RelAbund.Genus.Filt.zerofilt <- subset(control.RelAbund.Genus.Filt.zerofilt,
                                               Time_Point != "Day 14")

# create factors
factors_control_T1T2_pig_genus <- droplevels(control.T1T2.RelAbund.Genus.Filt.zerofilt %>%
                                               select(Time_Point, Pig))

# create permutations
perm_control_T1T2_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_control_T1T2_pig_genus$Pig,
                                                 type = "free"))

# run PERMANOVA
Control.T1T2.Genus.zerofilt.permanova <- adonis2(control.T1T2.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                                 data = factors_control_T1T2_pig_genus,
                                                 permutations = perm_control_T1T2_pig_genus, 
                                                 method = "bray",
                                                 by = "margin")

Control.T1T2.Genus.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_T1T2_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = control.T1T2.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_T1T2_pig_genus, permutations = perm_control_T1T2_pig_genus, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)  
    ## Time_Point  1  0.03492 0.08631 1.7003   0.03 *
    ## Residual   18  0.36971 0.91369                
    ## Total      19  0.40464 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Significant p = .030

##### Control T1 vs T3

``` r
set.seed(2021)
# filter data set for only samples at T1 and T3
control.T1T3.RelAbund.Genus.Filt.zerofilt <- subset(control.RelAbund.Genus.Filt.zerofilt,
                                               Time_Point != "Day 7")

# create factors
factors_control_T1T3_pig_genus <- droplevels(control.T1T3.RelAbund.Genus.Filt.zerofilt %>%
                                               select(Time_Point, Pig))

# create permutations
perm_control_T1T3_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_control_T1T3_pig_genus$Pig,
                                                 type = "free"))

# run PERMANOVA
Control.T1T3.Genus.zerofilt.permanova <- adonis2(control.T1T3.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                                 data = factors_control_T1T3_pig_genus,
                                                 permutations = perm_control_T1T3_pig_genus, 
                                                 method = "bray",
                                                 by = "margin")

Control.T1T3.Genus.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_T1T3_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = control.T1T3.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_T1T3_pig_genus, permutations = perm_control_T1T3_pig_genus, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## Time_Point  1  0.09340 0.28463 7.1617   0.01 **
    ## Residual   18  0.23474 0.71537                 
    ## Total      19  0.32814 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Significant p = .010

##### Control T2 vs T3

``` r
set.seed(2021)
# filter data set for only samples at T1 and T3
control.T2T3.RelAbund.Genus.Filt.zerofilt <- subset(control.RelAbund.Genus.Filt.zerofilt,
                                               Time_Point != "Day 0")

# create factors
factors_control_T2T3_pig_genus <- droplevels(control.T2T3.RelAbund.Genus.Filt.zerofilt %>%
                                               select(Time_Point, Pig))

# create permutations
perm_control_T2T3_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_control_T2T3_pig_genus$Pig,
                                                 type = "free"))

# run PERMANOVA
Control.T2T3.Genus.zerofilt.permanova <- adonis2(control.T2T3.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                                 data = factors_control_T2T3_pig_genus,
                                                 permutations = perm_control_T2T3_pig_genus, 
                                                 method = "bray",
                                                 by = "margin")

Control.T2T3.Genus.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_T2T3_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = control.T2T3.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_T2T3_pig_genus, permutations = perm_control_T2T3_pig_genus, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)  
    ## Time_Point  1  0.07429 0.18752 4.1544  0.015 *
    ## Residual   18  0.32189 0.81248                
    ## Total      19  0.39618 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Sig: p=.015

#### Tomato

Effect of tomato diet over time.

``` r
set.seed(2021)
# filter for only tomato
tomato.RelAbund.Genus.Filt.zerofilt <- subset(RelAbund.Genus.Filt.zerofilt, Diet == "Tomato")

# create factors
factors_tomato_genera <- droplevels(tomato.RelAbund.Genus.Filt.zerofilt %>% select(Time_Point, Pig))

# create permutations
perm_tomato_genera <- how(within = Within(type="series", constant=TRUE),
                          plots = Plots(strata=factors_tomato_genera$Pig, type="none",))
# run permanova
Tomato.ByTime.Genus.zerofilt <- adonis2(tomato.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                        data = factors_tomato_genera,
                                        permutations = perm_tomato_genera,
                                        method = "bray",
                                        by = "margin")
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

``` r
Tomato.ByTime.Genus.zerofilt
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_tomato_genera$Pig, plot permutation: none
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 2
    ## 
    ## adonis2(formula = tomato.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_tomato_genera, permutations = perm_tomato_genera, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)
    ## Time_Point  2  0.06217 0.11739 1.7955 0.3333
    ## Residual   27  0.46745 0.88261              
    ## Total      29  0.52962 1.00000

Significant effect of time (p = 0.01) within tomato samples

Now do pairwise comparisons to see where the significance is coming from

##### Tomato T1 vs Tomato T2

``` r
set.seed(2021)
# filter data set for only samples at T1 and T2
tomato.T1T2.RelAbund.Genus.Filt.zerofilt <- subset(tomato.RelAbund.Genus.Filt.zerofilt,
                                               Time_Point != "Day 14")

# create factors
factors_tomato_T1T2_pig_genus <- droplevels(tomato.T1T2.RelAbund.Genus.Filt.zerofilt %>%
                                               select(Time_Point, Pig))

# create permutations
perm_tomato_T1T2_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_tomato_T1T2_pig_genus$Pig,
                                                 type = "free"))

# run PERMANOVA
tomato.T1T2.Genus.zerofilt.permanova <- adonis2(tomato.T1T2.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                                 data = factors_tomato_T1T2_pig_genus,
                                                 permutations = perm_tomato_T1T2_pig_genus, 
                                                 method = "bray",
                                                 by = "margin")

tomato.T1T2.Genus.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_tomato_T1T2_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = tomato.T1T2.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_tomato_T1T2_pig_genus, permutations = perm_tomato_T1T2_pig_genus, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)  
    ## Time_Point  1  0.02556 0.07483 1.4559   0.09 .
    ## Residual   18  0.31598 0.92517                
    ## Total      19  0.34154 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Not significant .090

##### tomato T1 vs T3

``` r
set.seed(2021)
# filter data set for only samples at T1 and T3
tomato.T1T3.RelAbund.Genus.Filt.zerofilt <- subset(tomato.RelAbund.Genus.Filt.zerofilt,
                                               Time_Point != "Day 7")

# create factors
factors_tomato_T1T3_pig_genus <- droplevels(tomato.T1T3.RelAbund.Genus.Filt.zerofilt %>%
                                               select(Time_Point, Pig))

# create permutations
perm_tomato_T1T3_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_tomato_T1T3_pig_genus$Pig,
                                                 type = "free"))

# run PERMANOVA
tomato.T1T3.Genus.zerofilt.permanova <- adonis2(tomato.T1T3.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                                 data = factors_tomato_T1T3_pig_genus,
                                                 permutations = perm_tomato_T1T3_pig_genus, 
                                                 method = "bray",
                                                 by = "margin")

tomato.T1T3.Genus.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_tomato_T1T3_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = tomato.T1T3.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_tomato_T1T3_pig_genus, permutations = perm_tomato_T1T3_pig_genus, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)
    ## Time_Point  1  0.04969 0.14718 3.1064   0.15
    ## Residual   18  0.28793 0.85282              
    ## Total      19  0.33762 1.00000

Significant p = .150

##### tomato T2 vs T3

``` r
set.seed(2021)
# filter data set for only samples at T1 and T3
tomato.T2T3.RelAbund.Genus.Filt.zerofilt <- subset(tomato.RelAbund.Genus.Filt.zerofilt,
                                               Time_Point != "Day 0")

# create factors
factors_tomato_T2T3_pig_genus <- droplevels(tomato.T2T3.RelAbund.Genus.Filt.zerofilt %>%
                                               select(Time_Point, Pig))

# create permutations
perm_tomato_T2T3_pig_genus <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_tomato_T2T3_pig_genus$Pig,
                                                 type = "free"))

# run PERMANOVA
tomato.T2T3.Genus.zerofilt.permanova <- adonis2(tomato.T2T3.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Time_Point,
                                                 data = factors_tomato_T2T3_pig_genus,
                                                 permutations = perm_tomato_T2T3_pig_genus, 
                                                 method = "bray",
                                                 by = "margin")

tomato.T2T3.Genus.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_tomato_T2T3_pig_genus$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = tomato.T2T3.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_tomato_T2T3_pig_genus, permutations = perm_tomato_T2T3_pig_genus, method = "bray", by = "margin")
    ##            Df SumOfSqs     R2      F Pr(>F)
    ## Time_Point  1  0.01801 0.0516 0.9794   0.12
    ## Residual   18  0.33098 0.9484              
    ## Total      19  0.34899 1.0000

p is non significant =.120

### Subset by time

#### Day 0

Effect of diet at day 0.

``` r
# filter for day 0 only
d0.RelAbund.Genus.Filt.zerofilt <- subset(RelAbund.Genus.Filt.zerofilt, Time_Point == "Day 0")

# create factors
# don't need to include pig, since no repeated measures here 
# only testing Diet within a time point
factors_day0_genera <- d0.RelAbund.Genus.Filt.zerofilt %>% 
  select(Diet)

# create permutations
perm_day0_genera <- how(nperm = 9999)

# run PERMANOVA
d0.ByTime.Genus.zerofilt <- adonis2(d0.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Diet,
                                         data = factors_day0_genera,
                                         permutations = perm_day0_genera,
                                         method = "bray")
d0.ByTime.Genus.zerofilt
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = d0.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Diet, data = factors_day0_genera, permutations = perm_day0_genera, method = "bray")
    ##          Df SumOfSqs     R2      F Pr(>F)
    ## Diet      1 0.020138 0.0676 1.3051 0.2462
    ## Residual 18 0.277747 0.9324              
    ## Total    19 0.297885 1.0000

Non-significant effect of Diet (p = 0.2402) at day 0.

#### Day 7

Effect of diet at day 7.

``` r
# filter for day 7 only
d7.RelAbund.Genus.Filt.zerofilt <- subset(RelAbund.Genus.Filt.zerofilt, Time_Point == "Day 7")

# create factors
# don't need to include pig, since no repeated measures here 
# only testing Diet within a time point
factors_day7_genera <- d7.RelAbund.Genus.Filt.zerofilt %>% 
  select(Diet)

# create permutations
perm_day7_genera <- how(nperm = 9999)

# run PERMANOVA
d7.ByTime.Genus.zerofilt <- adonis2(d7.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Diet,
                                         data = factors_day7_genera,
                                         permutations = perm_day7_genera,
                                         method = "bray")
d7.ByTime.Genus.zerofilt
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = d7.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Diet, data = factors_day7_genera, permutations = perm_day7_genera, method = "bray")
    ##          Df SumOfSqs     R2      F Pr(>F)
    ## Diet      1  0.02780 0.0638 1.2267 0.2762
    ## Residual 18  0.40795 0.9362              
    ## Total    19  0.43575 1.0000

Non-significant effect of Diet (p = 0.2836) at day 7.

#### Day 14

Effect of diet at day 14.

``` r
# filter for day 14 only
d14.RelAbund.Genus.Filt.zerofilt <- subset(RelAbund.Genus.Filt.zerofilt, Time_Point == "Day 14")

# create factors
# don't need to include pig, since no repeated measures here 
# only testing Diet within a time point
factors_day14_genera <- d14.RelAbund.Genus.Filt.zerofilt %>% 
  select(Diet)

# create permutations
perm_day14_genera <- how(nperm = 9999)

# run PERMANOVA
d14.ByTime.Genus.zerofilt <- adonis2(d14.RelAbund.Genus.Filt.zerofilt[,-c(1:5)]~Diet,
                                     data = factors_day14_genera,
                                     permutations = perm_day14_genera,
                                     method = "bray")
d14.ByTime.Genus.zerofilt
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = d14.RelAbund.Genus.Filt.zerofilt[, -c(1:5)] ~ Diet, data = factors_day14_genera, permutations = perm_day14_genera, method = "bray")
    ##          Df SumOfSqs      R2     F Pr(>F)   
    ## Diet      1 0.041978 0.14631 3.085 0.0062 **
    ## Residual 18 0.244923 0.85369                
    ## Total    19 0.286900 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Significant effect of Diet (p = 0.005) at day 14.

## PCoA Beta Diversity

#### All samples

``` r
# calculate distances
genus.filt.dist.20zeros <- vegdist(RelAbund.Genus.Filt.zerofilt[6:ncol(RelAbund.Genus.Filt.zerofilt)], 
                                   method = "bray")

# do multi-dimensional scaling (the PCoA calculations) on those distances
scale.genus.filt.20zeros <- cmdscale(genus.filt.dist.20zeros, k=2)

# make into data frame
scale.genus.filt.df.20zeros <- as.data.frame(cbind(scale.genus.filt.20zeros, 
                                                   AllSamples.Metadata))

# do PCoA again, but get eigen values
scale.genus.filt.20zeros.eig <- cmdscale(genus.filt.dist.20zeros, k=2, eig = TRUE)

# convert eigenvalues to percentages and assign to a variable
eigs.genus.filt.20zeros <- (100* ((scale.genus.filt.20zeros.eig$eig)/(sum(scale.genus.filt.20zeros.eig$eig))))

# round the converted eigenvalues
round.eigs.genus.20zeros <- round(eigs.genus.filt.20zeros, 3)
```

All samples, one PCoA

``` r
PCoA_genera_20zeros_allsamples <- scale.genus.filt.df.20zeros %>%
ggplot(aes(x = `1`, y = `2`, fill = Diet_By_Time_Point)) +
  geom_point(size=3, color = "black", shape = 21, alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4", "sienna1","firebrick3","tomato4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))+
  labs(x=paste("PC1: ", round.eigs.genus.20zeros[1], "%"), 
       y=paste("PC2: ", round.eigs.genus.20zeros[2], "%"), 
       fill="Diet & Time Point",
       title = "Beta Diversity",
       subtitle = "Genus Level") 

PCoA_genera_20zeros_allsamples
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

``` r
ggsave("Figures/BetaDiversity_PCoA_Genera_allsamples.png", 
       plot = PCoA_genera_20zeros_allsamples, 
       dpi = 800, 
       width = 10, 
       height = 8)
```

Re-level factors

``` r
scale.genus.filt.df.20zeros <- scale.genus.filt.df.20zeros %>% 
  mutate(Time_Point = fct_relevel(Time_Point, c("Day 0", "Day 7", "Day 14")))
```

##### Facet by time point

``` r
PCoA_genera_20zeros_facetbytime <- scale.genus.filt.df.20zeros %>%
ggplot(aes(x = `1`, y = `2`, fill = Diet_By_Time_Point)) +
  geom_hline(yintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_point(size=3, color = "black", shape = 21, alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4", "sienna1","firebrick3","tomato4")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x=paste("PC1: ", round.eigs.genus.20zeros[1], "%"), 
       y=paste("PC2: ", round.eigs.genus.20zeros[2], "%"), 
       fill="Diet & Time Point",
       title = "Beta Diversity",
       subtitle = "Genera Level, Subset by Time Point") +
  facet_wrap(~Time_Point)

PCoA_genera_20zeros_facetbytime
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

``` r
ggsave("Figures/BetaDiversity_PCoA_Genera_FacetByTimePoint.png", 
       plot = PCoA_genera_20zeros_facetbytime, 
       dpi = 800, 
       width = 10, 
       height = 6)
```

##### Facet by diet

``` r
PCoA_genera_20zeros_facetbydiet <- scale.genus.filt.df.20zeros %>%
ggplot(aes(x = `1`, y = `2`, fill = Diet_By_Time_Point)) +
  geom_hline(yintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_point(size=3, color = "black", shape = 21, alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4", "sienna1","firebrick3","tomato4")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x=paste("PC1: ", round.eigs.genus.20zeros[1], "%"), 
       y=paste("PC2: ", round.eigs.genus.20zeros[2], "%"), 
       fill="Diet & Time Point",
       title = "Beta Diversity",
       subtitle = "Genera Level, Subset by Diet") +
  facet_wrap(~Diet)

PCoA_genera_20zeros_facetbydiet
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->

``` r
ggsave("Figures/BetaDiversity_PCoA_Genera_FacetByDiet.png", 
       plot = PCoA_genera_20zeros_facetbydiet, 
       dpi = 800, 
       width = 10, 
       height = 6)
```

### Subset

Ended up not using this as part of the paper. Since the input is
different here (i.e., the PCoA only has the subset data as an input) the
output looks slightly different and we didn’t feel this was the most
accurate depction of the data.

#### Control only

``` r
# calculate distances
control.RelAbund.Genus.Filt.zerofilt.dist <- vegdist(control.RelAbund.Genus.Filt.zerofilt[,-c(1:5)], 
                                                     method = "bray")

# calculate to make PCoA
control.scale.genus.filt.20zeros <- cmdscale(control.RelAbund.Genus.Filt.zerofilt.dist, k=2)

# filter metadata
meta.control <- subset(AllSamples.Metadata, Diet == "Control")

# make into data frame and add metadata
control.scale.genus.filt.20zeros.df <- as.data.frame(cbind(meta.control, control.scale.genus.filt.20zeros))

# get eigenvalues
control.scale.genus.filt.20zeros.eig <- cmdscale(control.RelAbund.Genus.Filt.zerofilt.dist, k=2, eig = TRUE)
control.eigs.genus.filt.20zeros <- (100*((control.scale.genus.filt.20zeros.eig$eig)/(sum(control.scale.genus.filt.20zeros.eig$eig))))
control.round.eigs.genus.20zeros <- round(control.eigs.genus.filt.20zeros, 3)
```

Reset factor levels

``` r
control.scale.genus.filt.20zeros.df$Time_Point <- factor(control.scale.genus.filt.20zeros.df$Time_Point, levels = c("Day 0", "Day 7", "Day 14"))
```

Plot

``` r
control.scale.genus.filt.20zeros.df %>%
ggplot(aes(x = `1`, y = `2`, fill = Time_Point)) +
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4"))+
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=(paste(control.round.eigs.genus.20zeros[1], "%")), 
       y=(paste(control.round.eigs.genus.20zeros[2], "%")), 
       fill = "Time Point",
       title = "Beta Diversity",
       subtitle = "Genera Level, Control Samples Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

#### Tomato only

``` r
# calculate distances
tomato.RelAbund.Genus.Filt.zerofilt.dist <- vegdist(tomato.RelAbund.Genus.Filt.zerofilt[,-c(1:5)], method = "bray")

# calculate to make PCoA
tomato.scale.genus.filt.20zeros <- cmdscale(tomato.RelAbund.Genus.Filt.zerofilt.dist, k=2)

# filter metadata
meta.tomato <- subset(AllSamples.Metadata, Diet == "Tomato")

# make into data frame and add metadata
tomato.scale.genus.filt.20zeros.df <- as.data.frame(cbind(meta.tomato, tomato.scale.genus.filt.20zeros))

# get eigenvalues
tomato.scale.genus.filt.20zeros.eig <- cmdscale(tomato.RelAbund.Genus.Filt.zerofilt.dist, k=2, eig = TRUE)
tomato.eigs.genus.filt.20zeros <- (100*((tomato.scale.genus.filt.20zeros.eig$eig)/(sum(tomato.scale.genus.filt.20zeros.eig$eig))))
tomato.round.eigs.genus.20zeros <- round(tomato.eigs.genus.filt.20zeros, 3)
```

Reset factor levels

``` r
tomato.scale.genus.filt.20zeros.df$Time_Point <- factor(tomato.scale.genus.filt.20zeros.df$Time_Point, levels = c("Day 0", "Day 7", "Day 14"))
```

Plot

``` r
tomato.scale.genus.filt.20zeros.df %>%
ggplot(aes(x = `1`, y = `2`, fill = Time_Point))+
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values = c("sienna1","firebrick3","tomato4"))+
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=(paste(tomato.round.eigs.genus.20zeros[1], "%")), 
       y=(paste(tomato.round.eigs.genus.20zeros[2], "%")), 
       fill = "Time Point",
       title = "Beta Diversity",
       subtitle = "Genera Level, Tomato Samples Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

#### Day 0 Only

``` r
# calculate distances
d0.RelAbund.Genus.Filt.zerofilt.dist <- vegdist(d0.RelAbund.Genus.Filt.zerofilt[,-c(1:5)], method = "bray")

# calculate to make PCoA
d0.scale.genus.filt.20zeros <- cmdscale(d0.RelAbund.Genus.Filt.zerofilt.dist, k=2)

# filter metadata
meta.day0 <- subset(AllSamples.Metadata, Time_Point == "Day 0")

# make into data frame and add metadata
d0.scale.genus.filt.20zeros.df <- as.data.frame(cbind(meta.day0, d0.scale.genus.filt.20zeros))

# get eigenvalues
d0.scale.genus.filt.20zeros.eig <- cmdscale(d0.RelAbund.Genus.Filt.zerofilt.dist, k=2, eig = TRUE)
d0.eigs.genus.filt.20zeros <- (100*((d0.scale.genus.filt.20zeros.eig$eig)/(sum(d0.scale.genus.filt.20zeros.eig$eig))))
d0.round.eigs.genus.20zeros <- round(d0.eigs.genus.filt.20zeros, 3)
```

Plot

``` r
d0.scale.genus.filt.20zeros.df %>%
ggplot(aes( x= `1`, y = `2`, fill = Diet))+
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values = c("steelblue2", "tomato2")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=(paste(d0.round.eigs.genus.20zeros[1], "%")), 
       y=(paste(d0.round.eigs.genus.20zeros[2], "%")),
       title = "Beta Diversity",
       subtitle = "Genera Level, Day 0 Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->

#### Day 7 Only

``` r
# calculate distances
d7.RelAbund.Genus.Filt.zerofilt.dist <- vegdist(d7.RelAbund.Genus.Filt.zerofilt[,-c(1:5)], method = "bray")

# calculate to make PCoA
d7.scale.genus.filt.20zeros <- cmdscale(d7.RelAbund.Genus.Filt.zerofilt.dist, k=2)

# filter metadata
meta.day7 <- subset(AllSamples.Metadata, Time_Point == "Day 7")

# make into data frame and add metadata
d7.scale.genus.filt.20zeros.df <- as.data.frame(cbind(meta.day7, d7.scale.genus.filt.20zeros))

# get eigenvalues
d7.scale.genus.filt.20zeros.eig <- cmdscale(d7.RelAbund.Genus.Filt.zerofilt.dist, k=2, eig = TRUE)
d7.eigs.genus.filt.20zeros <- (100*((d7.scale.genus.filt.20zeros.eig$eig)/(sum(d7.scale.genus.filt.20zeros.eig$eig))))
d7.round.eigs.genus.20zeros <- round(d7.eigs.genus.filt.20zeros, 3)
```

Plot

``` r
d7.scale.genus.filt.20zeros.df %>%
ggplot(aes( x= `1`, y = `2`, fill = Diet))+
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_color_manual(values = c("steelblue2", "tomato2")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=(paste(d7.round.eigs.genus.20zeros[1], "%")), 
       y=(paste(d7.round.eigs.genus.20zeros[2], "%")),
       title = "Beta Diversity",
       subtitle = "Genera Level, Day 7 Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-68-1.png)<!-- -->

#### Day 14 Only

``` r
# calculate distances
d14.RelAbund.Genus.Filt.zerofilt.dist <- vegdist(d14.RelAbund.Genus.Filt.zerofilt[,-c(1:5)], method = "bray")
# calculate to make PCoA
d14.scale.genus.filt.20zeros <- cmdscale(d14.RelAbund.Genus.Filt.zerofilt.dist, k=2)

# filter metadata
meta.day14 <- subset(AllSamples.Metadata, Time_Point == "Day 14")

# make into data frame and add metadata
d14.scale.genus.filt.20zeros.df <- as.data.frame(cbind(meta.day14, d14.scale.genus.filt.20zeros))

# get eigenvalues
d14.scale.genus.filt.20zeros.eig <- cmdscale(d14.RelAbund.Genus.Filt.zerofilt.dist, k=2, eig = TRUE)
d14.eigs.genus.filt.20zeros <- (100*((d14.scale.genus.filt.20zeros.eig$eig)/(sum(d14.scale.genus.filt.20zeros.eig$eig))))
d14.round.eigs.genus.20zeros <- round(d14.eigs.genus.filt.20zeros, 3)
```

Plot

``` r
d14.scale.genus.filt.20zeros.df %>%
ggplot(aes( x= `1`, y = `2`, fill = Diet))+
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_color_manual(values = c("steelblue2", "tomato2")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=(paste(d14.round.eigs.genus.20zeros[1], "%")), 
       y=(paste(d0.round.eigs.genus.20zeros[2], "%")),
       title = "Beta Diversity",
       subtitle = "Genera Level, Day 14 Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->

## Alpha Diversity

### Wrangling

``` r
kable(head(RelAbund.Genus.Filt.zerofilt))
```

| Sample\_Name                              | Pig | Diet    | Time\_Point | Diet\_By\_Time\_Point | Abiotrophia | Acaryochloris | Acetivibrio | Acetobacter | Acetohalobium | Acholeplasma | Achromobacter | Acidaminococcus | Acidilobus | Acidimicrobium | Acidiphilium | Acidithiobacillus | Acidobacterium | Acidothermus | Acidovorax | Aciduliprofundum | Acinetobacter | Actinobacillus | Actinomyces | Actinosynnema | Aerococcus | Aeromicrobium | Aeromonas | Aeropyrum |   Afipia | Aggregatibacter | Agrobacterium | Ahrensia | Ajellomyces | Akkermansia | Albidiferax | Alcanivorax | Algoriphagus | Alicycliphilus | Alicyclobacillus | Aliivibrio | Alistipes | Alkalilimnicola | Alkaliphilus | Allochromatium | Alphatorquevirus | Alteromonas | Aminobacterium | Aminomonas | Ammonifex | Amycolatopsis |  Anabaena | Anaerobaculum | Anaerococcus | Anaerofustis | Anaeromyxobacter | Anaerostipes | Anaerotruncus | Anaplasma | Anoxybacillus |   Aquifex | Arcanobacterium | Archaeoglobus | Arcobacter | Aromatoleum | Arthrobacter | Arthroderma | Arthrospira | Aspergillus | Asticcacaulis | Atopobium | Aurantimonas |  Azoarcus | Azorhizobium | Azospirillum | Azotobacter | Babesia |  Bacillus | Bacteroides | Bartonella |    Basfia | Bdellovibrio | Beggiatoa | Beijerinckia | Bermanella | Beutenbergia | Bifidobacterium | Blastopirellula | Blattabacterium |   Blautia | Bordetella |  Borrelia | Botryotinia | Bpp-1-like viruses | Brachybacterium | Brachyspira | Bradyrhizobium | Brevibacillus | Brevibacterium | Brevundimonas |  Brucella |   Brugia | Buchnera | Bulleidia | Burkholderia | Butyrivibrio | Caenorhabditis | Caldanaerobacter | Caldicellulosiruptor | Calditerrivibrio | Caldivirga | Caminibacter | Campylobacter |  Candida | Candidatus Accumulibacter | Candidatus Amoebophilus | Candidatus Azobacteroides | Candidatus Blochmannia | Candidatus Cloacamonas | Candidatus Desulforudis | Candidatus Hamiltonella | Candidatus Korarchaeum | Candidatus Koribacter | Candidatus Liberibacter | Candidatus Pelagibacter | Candidatus Phytoplasma | Candidatus Protochlamydia | Candidatus Puniceispirillum | Candidatus Regiella | Candidatus Riesia | Candidatus Solibacter | Candidatus Sulcia | Capnocytophaga | Carboxydothermus | Cardiobacterium | Carnobacterium | Catenibacterium | Catenulispora | Catonella | Caulobacter | Cellulomonas | Cellulosilyticum | Cellvibrio | Cenarchaeum | Chaetomium | Chelativorans | Chitinophaga | Chlamydia | Chlamydomonas | Chlamydophila | Chlorella | Chlorobaculum | Chlorobium | Chloroflexus | Chloroherpeton | Chlorovirus | Chromobacterium | Chromohalobacter | Chryseobacterium | Chthoniobacter | Citreicella | Citrobacter | Citromicrobium | Clavibacter | Clavispora | Clostridium | Coccidioides | Collinsella | Colwellia | Comamonas | Conexibacter | Congregibacter | Coprinopsis | Coprobacillus | Coprococcus | Coprothermobacter | Coraliomargarita | Corynebacterium | Coxiella | Croceibacter | Crocosphaera | Cronobacter | Cryptobacterium | Cryptosporidium | Cupriavidus | Cyanidioschyzon | Cyanidium | Cyanobium | Cyanophora | Cyanothece | Cylindrospermopsis | Cytophaga | Debaryomyces | Dechloromonas | Deferribacter | Dehalococcoides | Dehalogenimonas | Deinococcus |   Delftia | Denitrovibrio | Dermacoccus | Desulfarculus | Desulfatibacillum | Desulfitobacterium | Desulfobacterium | Desulfococcus | Desulfohalobium | Desulfomicrobium | Desulfonatronospira | Desulfotalea | Desulfotomaculum | Desulfovibrio | Desulfurispirillum | Desulfurivibrio | Desulfurococcus | Desulfuromonas | Dethiobacter | Dethiosulfovibrio | Dialister | Dichelobacter |   Dickeya | Dictyoglomus | Dictyostelium | Dinoroseobacter |  Dokdonia |     Dorea | Dyadobacter | Edwardsiella | Eggerthella | Ehrlichia | Eikenella | Elusimicrobium | Emericella | Emiliania | Encephalitozoon | Endoriftia | Enhydrobacter | Entamoeba | Enterobacter | Enterococcus | Enterocytozoon | Epsilon15-like viruses | Epulopiscium | Eremococcus | Eremothecium |   Erwinia | Erysipelothrix | Erythrobacter | Escherichia | Ethanoligenens | Eubacterium | Exiguobacterium | Faecalibacterium | Ferrimonas | Ferroglobus | Ferroplasma | Fervidobacterium | Fibrobacter | Filifactor | Filobasidiella | Finegoldia | Flavobacterium | Francisella |   Frankia | Fulvimarina | Fusobacterium | Gallionella | Gammaretrovirus | Gardnerella |  Gemella |  Gemmata | Gemmatimonas | Geobacillus | Geobacter | Geodermatophilus |  Giardia | Gibberella | Gloeobacter | Gluconacetobacter | Gluconobacter | Gordonia | Gracilaria |  Gramella | Granulibacter | Granulicatella | Guillardia | Haemophilus |   Hahella | Halalkalicoccus | Halanaerobium | Haliangium | Haloarcula | Halobacterium | Haloferax | Halogeometricum | Halomicrobium | Halomonas | Haloquadratum | Halorhabdus | Halorhodospira | Halorubrum | Haloterrigena | Halothermothrix | Halothiobacillus | Helicobacter | Heliobacterium | Herbaspirillum | Herminiimonas | Herpetosiphon | Hirschia | Histophilus |  Hoeflea | Holdemania | Hydrogenivirga | Hydrogenobacter | Hydrogenobaculum | Hyperthermus | Hyphomicrobium | Hyphomonas | Idiomarina | Ignicoccus | Ignisphaera | Ilyobacter | Intrasporangium | Janibacter | Jannaschia | Janthinobacterium |   Jonesia | Jonquetella | Kangiella | Ketogulonicigenium | Kineococcus | Kingella | Klebsiella | Kluyveromyces |   Kocuria |   Kordia | Kosmotoga | Kribbella | Ktedonobacter | Kytococcus | L5-like viruses | Labrenzia | Laccaria | Lachancea | Lactobacillus | Lactococcus | Lambda-like viruses | Laribacter |  Lawsonia | Leadbetterella | Leeuwenhoekiella | Legionella | Leifsonia | Leishmania | Lentisphaera | Leptosira | Leptospira | Leptospirillum | Leptothrix | Leptotrichia | Leuconostoc | Limnobacter |  Listeria |      Loa | Lodderomyces | Loktanella | Lutiella |  Lyngbya | Lysinibacillus | Macrococcus | Magnaporthe | Magnetococcus | Magnetospirillum | Malassezia | Mannheimia | Maribacter | Maricaulis | Marinobacter | Marinomonas | Mariprofundus | Maritimibacter | Marivirga | Megasphaera | Meiothermus | Mesoplasma | Mesorhizobium | Metallosphaera | Methanobrevibacter | Methanocaldococcus | Methanocella | Methanococcoides | Methanococcus | Methanocorpusculum | Methanoculleus | Methanohalobium | Methanohalophilus | Methanoplanus | Methanopyrus | Methanoregula | Methanosaeta | Methanosarcina | Methanosphaera | Methanosphaerula | Methanospirillum | Methanothermobacter | Methanothermococcus | Methanothermus | Methylacidiphilum | Methylibium | Methylobacillus | Methylobacter | Methylobacterium | Methylocella | Methylococcus | Methylophaga | Methylosinus | Methylotenera | Methylovorus | Meyerozyma | Micrococcus | Microcoleus | Microcystis | Micromonas | Micromonospora | Microscilla | Mitsuokella | Mobiluncus | Moniliophthora | Monosiga |  Moorella | Moraxella | Moritella | Mucilaginibacter | Mycobacterium | Mycoplasma | Myxococcus | N4-like viruses | Naegleria | Nakamurella | Nakaseomyces | Nanoarchaeum | Natranaerobius | Natrialba | Natronomonas | Nautilia | Nectria | Neisseria | Neorickettsia | Neosartorya | Neptuniibacter | Neurospora | Nitratiruptor | Nitrobacter | Nitrococcus | Nitrosococcus | Nitrosomonas | Nitrosopumilus | Nitrosospira | Nitrospira |  Nocardia | Nocardioides | Nocardiopsis | Nodularia |    Nostoc | Novosphingobium | Oceanibulbus | Oceanicaulis | Oceanicola | Oceanithermus | Oceanobacillus | Ochrobactrum | Octadecabacter | Odontella | Oenococcus | Oligotropha | Olsenella |  Opitutus | Oribacterium | Orientia | Oscillatoria | Oscillochloris | Ostreococcus | Oxalobacter | P1-like viruses | P2-like viruses | P22-like viruses | Paenibacillus | Paludibacter |   Pantoea | Parabacteroides | Parachlamydia | Paracoccidioides | Paracoccus | Paramecium | Parascardovia | Parvibaculum | Parvularcula | Pasteurella | Paulinella | Pectobacterium | Pediococcus | Pedobacter | Pelagibaca | Pelobacter | Pelodictyon | Pelotomaculum | Penicillium | Peptoniphilus | Peptostreptococcus | Perkinsus | Persephonella | Petrotoga | Phaeobacter | Phaeodactylum | Phaeosphaeria | Phenylobacterium | Phi29-like viruses | Photobacterium | Photorhabdus | Phytophthora |  Pichia | Picrophilus | Pirellula | Planctomyces | Plasmodium | Plesiocystis | Podospora | Polaribacter | Polaromonas | Polynucleobacter | Porphyra | Porphyromonas |  Postia | Prevotella | Prochlorococcus | Propionibacterium | Prosthecochloris |   Proteus | Providencia | Pseudoalteromonas | Pseudomonas | Pseudoramibacter | Pseudovibrio | Psychrobacter | Psychroflexus | Psychromonas | Pyramidobacter | Pyrenophora | Pyrobaculum | Pyrococcus | Ralstonia | Raphidiopsis | Reclinomonas | Reinekea | Renibacterium | Rhizobium | Rhodobacter | Rhodococcus | Rhodomicrobium | Rhodomonas | Rhodopirellula | Rhodopseudomonas | Rhodospirillum | Rhodothermus | Rickettsia | Rickettsiella | Riemerella | Robiginitalea | Roseburia | Roseibium | Roseiflexus | Roseobacter | Roseomonas | Roseovarius |    Rothia | Rubrobacter |  Ruegeria | Ruminococcus | SP6-like viruses | SPO1-like viruses | SPbeta-like viruses | Saccharomonospora | Saccharomyces | Saccharophagus | Saccharopolyspora | Saccoglossus | Sagittula | Salinibacter | Salinispora | Salmonella | Sanguibacter | Scardovia | Scheffersomyces | Schizophyllum | Schizosaccharomyces | Sclerotinia | Sebaldella | Segniliparus | Selenomonas |  Serratia | Shewanella |  Shigella | Shuttleworthia | Sideroxydans | Simonsiella | Sinorhizobium |   Slackia |  Sodalis | Sorangium | Sphaerobacter | Sphingobacterium | Sphingobium | Sphingomonas | Sphingopyxis | Spirochaeta | Spirosoma | Stackebrandtia | Staphylococcus | Staphylothermus | Starkeya | Stenotrophomonas | Stigmatella | Streptobacillus | Streptococcus | Streptomyces | Streptosporangium | Subdoligranulum | Sulfitobacter | Sulfolobus | Sulfuricurvum | Sulfurihydrogenibium | Sulfurimonas | Sulfurospirillum | Sulfurovum | Symbiobacterium | Synechococcus | Synechocystis | Syntrophobacter | Syntrophomonas | Syntrophothermus | Syntrophus | T4-like viruses | T7-like viruses | Talaromyces | Teredinibacter | Terriglobus | Tetragenococcus | Tetrahymena | Thalassiosira | Thalassobium |  Thauera | Theileria | Thermaerobacter | Thermanaerovibrio | Thermincola | Thermoanaerobacter | Thermoanaerobacterium | Thermobaculum | Thermobifida | Thermobispora | Thermococcus | Thermocrinis | Thermodesulfovibrio | Thermofilum | Thermomicrobium | Thermomonospora | Thermoplasma | Thermoproteus | Thermosediminibacter | Thermosinus | Thermosipho | Thermosphaera | Thermosynechococcus | Thermotoga |   Thermus | Thioalkalivibrio | Thiobacillus | Thiomicrospira | Thiomonas | Tolumonas | Toxoplasma | Treponema | Trichodesmium | Trichomonas | Trichoplax | Tropheryma | Truepera | Trypanosoma | Tsukamurella |   Tuber | Turicibacter | Uncinocarpus | Ureaplasma | Ustilago | Vanderwaltozyma | Variovorax | Veillonella | Verminephrobacter | Verrucomicrobium | Verticillium |    Vibrio | Victivallis |   Volvox | Vulcanisaeta |  Waddlia | Weissella | Wigglesworthia | Wolbachia | Wolinella | Xanthobacter | Xanthomonas | Xenorhabdus | Xylanimonas |   Xylella | Yarrowia |  Yersinia | Zunongwangia | Zygosaccharomyces | Zymomonas | phiKZ-like viruses | unclassified (derived from Actinobacteria (class)) | unclassified (derived from Alicyclobacillaceae) | unclassified (derived from Alphaproteobacteria) | unclassified (derived from Alteromonadales) | unclassified (derived from Bacteria) | unclassified (derived from Bacteroidetes) | unclassified (derived from Betaproteobacteria) | unclassified (derived from Burkholderiales) | unclassified (derived from Campylobacterales) | unclassified (derived from Candidatus Poribacteria) | unclassified (derived from Caudovirales) | unclassified (derived from Chroococcales) | unclassified (derived from Clostridiales Family XI. Incertae Sedis) | unclassified (derived from Clostridiales) | unclassified (derived from Deltaproteobacteria) | unclassified (derived from Elusimicrobia) | unclassified (derived from Erysipelotrichaceae) | unclassified (derived from Euryarchaeota) | unclassified (derived from Flavobacteria) | unclassified (derived from Flavobacteriaceae) | unclassified (derived from Flavobacteriales) | unclassified (derived from Gammaproteobacteria) | unclassified (derived from Lachnospiraceae) | unclassified (derived from Methylophilales) | unclassified (derived from Myoviridae) | unclassified (derived from Opitutaceae) | unclassified (derived from Podoviridae) | unclassified (derived from Rhodobacteraceae) | unclassified (derived from Rhodobacterales) | unclassified (derived from Rickettsiales) | unclassified (derived from Ruminococcaceae) | unclassified (derived from Siphoviridae) | unclassified (derived from Thermotogales) | unclassified (derived from Verrucomicrobia subdivision 3) | unclassified (derived from Verrucomicrobiales) | unclassified (derived from Vibrionaceae) | unclassified (derived from Vibrionales) | unclassified (derived from Viruses) | unclassified (derived from other sequences) |
|:------------------------------------------|:----|:--------|:------------|:----------------------|------------:|--------------:|------------:|------------:|--------------:|-------------:|--------------:|----------------:|-----------:|---------------:|-------------:|------------------:|---------------:|-------------:|-----------:|-----------------:|--------------:|---------------:|------------:|--------------:|-----------:|--------------:|----------:|----------:|---------:|----------------:|--------------:|---------:|------------:|------------:|------------:|------------:|-------------:|---------------:|-----------------:|-----------:|----------:|----------------:|-------------:|---------------:|-----------------:|------------:|---------------:|-----------:|----------:|--------------:|----------:|--------------:|-------------:|-------------:|-----------------:|-------------:|--------------:|----------:|--------------:|----------:|----------------:|--------------:|-----------:|------------:|-------------:|------------:|------------:|------------:|--------------:|----------:|-------------:|----------:|-------------:|-------------:|------------:|--------:|----------:|------------:|-----------:|----------:|-------------:|----------:|-------------:|-----------:|-------------:|----------------:|----------------:|----------------:|----------:|-----------:|----------:|------------:|-------------------:|----------------:|------------:|---------------:|--------------:|---------------:|--------------:|----------:|---------:|---------:|----------:|-------------:|-------------:|---------------:|-----------------:|---------------------:|-----------------:|-----------:|-------------:|--------------:|---------:|--------------------------:|------------------------:|--------------------------:|-----------------------:|-----------------------:|------------------------:|------------------------:|-----------------------:|----------------------:|------------------------:|------------------------:|-----------------------:|--------------------------:|----------------------------:|--------------------:|------------------:|----------------------:|------------------:|---------------:|-----------------:|----------------:|---------------:|----------------:|--------------:|----------:|------------:|-------------:|-----------------:|-----------:|------------:|-----------:|--------------:|-------------:|----------:|--------------:|--------------:|----------:|--------------:|-----------:|-------------:|---------------:|------------:|----------------:|-----------------:|-----------------:|---------------:|------------:|------------:|---------------:|------------:|-----------:|------------:|-------------:|------------:|----------:|----------:|-------------:|---------------:|------------:|--------------:|------------:|------------------:|-----------------:|----------------:|---------:|-------------:|-------------:|------------:|----------------:|----------------:|------------:|----------------:|----------:|----------:|-----------:|-----------:|-------------------:|----------:|-------------:|--------------:|--------------:|----------------:|----------------:|------------:|----------:|--------------:|------------:|--------------:|------------------:|-------------------:|-----------------:|--------------:|----------------:|-----------------:|--------------------:|-------------:|-----------------:|--------------:|-------------------:|----------------:|----------------:|---------------:|-------------:|------------------:|----------:|--------------:|----------:|-------------:|--------------:|----------------:|----------:|----------:|------------:|-------------:|------------:|----------:|----------:|---------------:|-----------:|----------:|----------------:|-----------:|--------------:|----------:|-------------:|-------------:|---------------:|-----------------------:|-------------:|------------:|-------------:|----------:|---------------:|--------------:|------------:|---------------:|------------:|----------------:|-----------------:|-----------:|------------:|------------:|-----------------:|------------:|-----------:|---------------:|-----------:|---------------:|------------:|----------:|------------:|--------------:|------------:|----------------:|------------:|---------:|---------:|-------------:|------------:|----------:|-----------------:|---------:|-----------:|------------:|------------------:|--------------:|---------:|-----------:|----------:|--------------:|---------------:|-----------:|------------:|----------:|----------------:|--------------:|-----------:|-----------:|--------------:|----------:|----------------:|--------------:|----------:|--------------:|------------:|---------------:|-----------:|--------------:|----------------:|-----------------:|-------------:|---------------:|---------------:|--------------:|--------------:|---------:|------------:|---------:|-----------:|---------------:|----------------:|-----------------:|-------------:|---------------:|-----------:|-----------:|-----------:|------------:|-----------:|----------------:|-----------:|-----------:|------------------:|----------:|------------:|----------:|-------------------:|------------:|---------:|-----------:|--------------:|----------:|---------:|----------:|----------:|--------------:|-----------:|----------------:|----------:|---------:|----------:|--------------:|------------:|--------------------:|-----------:|----------:|---------------:|-----------------:|-----------:|----------:|-----------:|-------------:|----------:|-----------:|---------------:|-----------:|-------------:|------------:|------------:|----------:|---------:|-------------:|-----------:|---------:|---------:|---------------:|------------:|------------:|--------------:|-----------------:|-----------:|-----------:|-----------:|-----------:|-------------:|------------:|--------------:|---------------:|----------:|------------:|------------:|-----------:|--------------:|---------------:|-------------------:|-------------------:|-------------:|-----------------:|--------------:|-------------------:|---------------:|----------------:|------------------:|--------------:|-------------:|--------------:|-------------:|---------------:|---------------:|-----------------:|-----------------:|--------------------:|--------------------:|---------------:|------------------:|------------:|----------------:|--------------:|-----------------:|-------------:|--------------:|-------------:|-------------:|--------------:|-------------:|-----------:|------------:|------------:|------------:|-----------:|---------------:|------------:|------------:|-----------:|---------------:|---------:|----------:|----------:|----------:|-----------------:|--------------:|-----------:|-----------:|----------------:|----------:|------------:|-------------:|-------------:|---------------:|----------:|-------------:|---------:|--------:|----------:|--------------:|------------:|---------------:|-----------:|--------------:|------------:|------------:|--------------:|-------------:|---------------:|-------------:|-----------:|----------:|-------------:|-------------:|----------:|----------:|----------------:|-------------:|-------------:|-----------:|--------------:|---------------:|-------------:|---------------:|----------:|-----------:|------------:|----------:|----------:|-------------:|---------:|-------------:|---------------:|-------------:|------------:|----------------:|----------------:|-----------------:|--------------:|-------------:|----------:|----------------:|--------------:|-----------------:|-----------:|-----------:|--------------:|-------------:|-------------:|------------:|-----------:|---------------:|------------:|-----------:|-----------:|-----------:|------------:|--------------:|------------:|--------------:|-------------------:|----------:|--------------:|----------:|------------:|--------------:|--------------:|-----------------:|-------------------:|---------------:|-------------:|-------------:|--------:|------------:|----------:|-------------:|-----------:|-------------:|----------:|-------------:|------------:|-----------------:|---------:|--------------:|--------:|-----------:|----------------:|------------------:|-----------------:|----------:|------------:|------------------:|------------:|-----------------:|-------------:|--------------:|--------------:|-------------:|---------------:|------------:|------------:|-----------:|----------:|-------------:|-------------:|---------:|--------------:|----------:|------------:|------------:|---------------:|-----------:|---------------:|-----------------:|---------------:|-------------:|-----------:|--------------:|-----------:|--------------:|----------:|----------:|------------:|------------:|-----------:|------------:|----------:|------------:|----------:|-------------:|-----------------:|------------------:|--------------------:|------------------:|--------------:|---------------:|------------------:|-------------:|----------:|-------------:|------------:|-----------:|-------------:|----------:|----------------:|--------------:|--------------------:|------------:|-----------:|-------------:|------------:|----------:|-----------:|----------:|---------------:|-------------:|------------:|--------------:|----------:|---------:|----------:|--------------:|-----------------:|------------:|-------------:|-------------:|------------:|----------:|---------------:|---------------:|----------------:|---------:|-----------------:|------------:|----------------:|--------------:|-------------:|------------------:|----------------:|--------------:|-----------:|--------------:|---------------------:|-------------:|-----------------:|-----------:|----------------:|--------------:|--------------:|----------------:|---------------:|-----------------:|-----------:|----------------:|----------------:|------------:|---------------:|------------:|----------------:|------------:|--------------:|-------------:|---------:|----------:|----------------:|------------------:|------------:|-------------------:|----------------------:|--------------:|-------------:|--------------:|-------------:|-------------:|--------------------:|------------:|----------------:|----------------:|-------------:|--------------:|---------------------:|------------:|------------:|--------------:|--------------------:|-----------:|----------:|-----------------:|-------------:|---------------:|----------:|----------:|-----------:|----------:|--------------:|------------:|-----------:|-----------:|---------:|------------:|-------------:|--------:|-------------:|-------------:|-----------:|---------:|----------------:|-----------:|------------:|------------------:|-----------------:|-------------:|----------:|------------:|---------:|-------------:|---------:|----------:|---------------:|----------:|----------:|-------------:|------------:|------------:|------------:|----------:|---------:|----------:|-------------:|------------------:|----------:|-------------------:|---------------------------------------------------:|------------------------------------------------:|------------------------------------------------:|--------------------------------------------:|-------------------------------------:|------------------------------------------:|-----------------------------------------------:|--------------------------------------------:|----------------------------------------------:|----------------------------------------------------:|-----------------------------------------:|------------------------------------------:|--------------------------------------------------------------------:|------------------------------------------:|------------------------------------------------:|------------------------------------------:|------------------------------------------------:|------------------------------------------:|------------------------------------------:|----------------------------------------------:|---------------------------------------------:|------------------------------------------------:|--------------------------------------------:|--------------------------------------------:|---------------------------------------:|----------------------------------------:|----------------------------------------:|---------------------------------------------:|--------------------------------------------:|------------------------------------------:|--------------------------------------------:|-----------------------------------------:|------------------------------------------:|----------------------------------------------------------:|-----------------------------------------------:|-----------------------------------------:|----------------------------------------:|------------------------------------:|--------------------------------------------:|
| ShotgunWGS-ControlPig6GutMicrobiome-Day14 | 6   | Control | Day 14      | Control Day 14        |   0.0013057 |      6.98e-05 |   0.0005123 |    1.70e-05 |     0.0002670 |    0.0002007 |     0.0000495 |       0.0129311 |    2.6e-06 |       1.52e-05 |     6.29e-05 |          6.31e-05 |      0.0000941 |    0.0001002 |  0.0002036 |        0.0000559 |     0.0002046 |      0.0006754 |   0.0002587 |      4.97e-05 |  0.0000879 |      6.40e-06 | 0.0002193 |   8.0e-06 | 5.90e-06 |       0.0000948 |     0.0001332 | 4.10e-06 |     2.8e-06 |   0.0003798 |   0.0000881 |   0.0000657 |    0.0000441 |       2.09e-05 |        0.0003706 |  0.0001340 | 0.0005110 |       0.0001064 |    0.0033871 |       5.10e-05 |         2.10e-06 |   0.0000381 |      0.0001312 |   2.60e-05 | 0.0002358 |      4.15e-05 | 0.0001219 |      3.35e-05 |    0.0011854 |    0.0003992 |        0.0003378 |    0.0011052 |     0.0017118 |  5.70e-06 |     0.0003584 | 0.0001428 |       0.0001368 |     0.0001188 |  0.0001322 |   0.0001062 |    0.0002451 |     3.6e-06 |    2.37e-05 |    4.30e-05 |      7.40e-05 | 0.0091799 |     2.89e-05 | 0.0000737 |     4.90e-05 |     4.25e-05 |    4.95e-05 | 2.1e-06 | 0.0070241 |   0.0821852 |   5.72e-05 | 0.0001917 |    0.0000861 |  2.09e-05 |     2.47e-05 |   7.00e-06 |    0.0001046 |       0.0105635 |       0.0000760 |        1.34e-05 | 0.0045629 |  0.0002440 | 0.0000750 |    4.90e-06 |           1.19e-05 |       0.0000776 |   0.0013351 |      0.0001835 |     0.0004561 |      0.0000719 |      3.25e-05 | 0.0000711 | 9.00e-06 | 2.22e-05 | 0.0005981 |    0.0006347 |    0.0103885 |       3.48e-05 |        0.0017562 |            0.0022862 |        0.0000760 |   7.20e-06 |     8.50e-06 |     0.0011851 | 1.08e-05 |                 0.0000523 |                5.05e-05 |                 0.0001809 |               1.31e-05 |               1.62e-05 |               0.0003760 |                1.24e-05 |               1.93e-05 |             0.0001729 |                 5.9e-06 |                1.57e-05 |               4.35e-05 |                  4.56e-05 |                    1.19e-05 |            4.60e-06 |           1.8e-06 |             0.0004154 |          3.60e-06 |      0.0005886 |        0.0011439 |        1.86e-05 |       5.98e-05 |       0.0012611 |      6.44e-05 |  4.56e-05 |   0.0002100 |    0.0000750 |        0.0006844 |  0.0001497 |    3.90e-06 |    7.7e-06 |     0.0000706 |    0.0005677 |  3.53e-05 |      3.27e-05 |      3.48e-05 |   0.0e+00 |     0.0001706 |  0.0006079 |    0.0002337 |      0.0001206 |     1.0e-06 |       0.0001036 |        0.0000902 |        0.0000845 |       2.19e-05 |    3.30e-06 |   0.0001175 |       4.10e-06 |   0.0001051 |   4.40e-06 |   0.0767302 |     3.10e-06 |   0.0049257 | 0.0000662 | 0.0000307 |     8.35e-05 |       3.68e-05 |    6.70e-06 |     0.0003621 |   0.0063405 |         0.0000866 |        0.0000523 |       0.0017554 | 5.00e-05 |    0.0001224 |     6.31e-05 |   0.0000956 |       0.0010161 |        6.70e-06 |   0.0002183 |         2.1e-06 |   6.4e-06 | 0.0000031 |    1.3e-06 |  0.0004015 |           1.42e-05 | 0.0004213 |      8.2e-06 |     0.0000923 |     0.0001430 |       0.0009210 |       0.0000850 |   0.0002654 | 0.0000569 |     0.0001660 |    8.00e-06 |     0.0000796 |         0.0001523 |          0.0037569 |        0.0001355 |     0.0001580 |       0.0000680 |        0.0001482 |            1.80e-05 |    0.0001946 |        0.0029103 |     0.0012480 |           8.53e-05 |        5.82e-05 |        4.40e-06 |      0.0001660 |    0.0001098 |         0.0002556 | 0.0062201 |      8.50e-05 | 0.0001167 |    0.0003136 |      4.51e-05 |        3.50e-05 | 0.0000943 | 0.0052538 |   0.0004803 |    0.0000763 |   0.0027436 |  1.65e-05 |  7.20e-06 |      0.0000995 |   1.39e-05 |   0.0e+00 |         5.0e-07 |   3.30e-06 |      8.20e-06 | 0.0000255 |    0.0001585 |    0.0027601 |       1.21e-05 |                1.0e-06 |    0.0001392 |    2.99e-05 |     9.30e-06 | 0.0000727 |      0.0000740 |     0.0000858 |   0.0005095 |      0.0045129 |   0.0713607 |       0.0005494 |        0.0345703 |  0.0000340 |    3.25e-05 |    2.45e-05 |        0.0002023 |   0.0012287 |  0.0001531 |       1.88e-05 |  0.0007628 |      0.0012769 |   0.0001080 | 0.0002448 |    1.03e-05 |     0.0018682 |    3.38e-05 |         3.0e-07 |   0.0006563 | 4.02e-05 | 2.37e-05 |    0.0000544 |   0.0018128 | 0.0014026 |         5.33e-05 | 6.20e-06 |   2.47e-05 |   0.0001211 |          2.89e-05 |      6.36e-05 | 3.22e-05 |    5.0e-07 | 0.0003646 |      4.17e-05 |      0.0001232 |    7.2e-06 |   0.0003079 | 0.0001116 |        1.11e-05 |     0.0002358 |  0.0000729 |   2.24e-05 |      1.65e-05 |  1.34e-05 |        1.57e-05 |      1.00e-05 |  3.12e-05 |      1.52e-05 |    1.49e-05 |       5.82e-05 |   1.06e-05 |      1.96e-05 |       0.0005736 |         2.53e-05 |    0.0005641 |      0.0013552 |       3.79e-05 |     0.0000572 |     0.0001337 | 3.14e-05 |   0.0001221 | 8.00e-06 |  0.0034288 |       1.78e-05 |        3.66e-05 |         6.78e-05 |     1.11e-05 |       2.01e-05 |   4.30e-05 |  0.0000742 |   1.11e-05 |     6.7e-06 |  0.0006200 |        3.84e-05 |   2.86e-05 |   3.30e-05 |         0.0000580 | 0.0000861 |   0.0000698 |  2.81e-05 |           1.75e-05 |   0.0000925 | 8.80e-06 |  0.0002023 |      9.80e-06 | 0.0000673 | 2.96e-05 | 0.0000956 |  5.57e-05 |      3.61e-05 |   4.17e-05 |         3.0e-07 |  1.57e-05 |  3.3e-06 |   7.5e-06 |     0.0460339 |   0.0006994 |           0.0000167 |   3.74e-05 | 0.0000863 |      0.0003785 |        0.0003412 |  0.0001433 | 0.0000557 |   2.19e-05 |     2.63e-05 |     5e-07 |  0.0001330 |        3.0e-07 |  0.0000559 |    0.0006277 |   0.0004314 |    8.00e-06 | 0.0011256 | 2.10e-06 |      3.3e-06 |   1.08e-05 | 2.09e-05 | 1.86e-05 |      0.0002981 |   0.0001296 |    1.49e-05 |     0.0001492 |        0.0001649 |   6.70e-06 |  0.0000982 |  0.0002020 |   3.53e-05 |    0.0001242 |   0.0001456 |      1.21e-05 |       2.27e-05 | 0.0001381 |   0.0101970 |   0.0001415 |  0.0000500 |     0.0000794 |       1.11e-05 |          0.0004927 |          0.0001075 |     4.61e-05 |        0.0001031 |     0.0003092 |          0.0003340 |      0.0000987 |        3.30e-05 |          3.50e-05 |     0.0000711 |     2.71e-05 |      7.32e-05 |    0.0000693 |      0.0004033 |      0.0001286 |         5.49e-05 |        0.0001090 |           0.0000902 |            2.60e-06 |       1.16e-05 |          2.63e-05 |   0.0000466 |        5.26e-05 |      1.03e-05 |        0.0001515 |     5.10e-05 |     0.0000964 |     3.30e-06 |     6.20e-06 |      3.63e-05 |     4.46e-05 |    6.2e-06 |    5.67e-05 |    1.88e-05 |    8.86e-05 |   1.24e-05 |       5.64e-05 |   0.0000796 |   0.0274527 |  0.0015438 |        3.9e-06 | 1.47e-05 | 0.0011094 |  2.37e-05 | 0.0000129 |        0.0001693 |     0.0004409 |  0.0003685 |  0.0001484 |         5.0e-07 |  8.20e-06 |    5.67e-05 |     1.83e-05 |      1.0e-06 |      0.0004064 |  8.20e-06 |     1.60e-05 | 4.07e-05 | 6.2e-06 | 0.0002317 |      1.42e-05 |    3.99e-05 |       9.80e-06 |   2.47e-05 |     0.0001008 |   0.0000835 |    3.12e-05 |     0.0001214 |    0.0001090 |       1.06e-05 |     6.60e-05 |   4.92e-05 | 0.0000941 |    0.0001098 |     4.59e-05 |  1.78e-05 | 0.0001775 |       0.0000698 |      5.4e-06 |     2.53e-05 |   2.94e-05 |      4.07e-05 |      0.0004453 |    0.0000706 |       1.37e-05 |   2.1e-06 |  0.0002000 |    3.71e-05 | 0.0066190 | 0.0002461 |    0.0050074 | 1.16e-05 |     1.08e-05 |       2.04e-05 |     2.37e-05 |   0.0001325 |       0.0000000 |         3.0e-07 |          0.0e+00 |     0.0020159 |    0.0012415 | 0.0001355 |       0.0063554 |      8.20e-06 |          1.0e-06 |   8.68e-05 |   1.08e-05 |     0.0002180 |     6.62e-05 |     2.19e-05 |   0.0001057 |    1.3e-06 |      0.0001788 |   0.0004473 |  0.0008558 |   6.40e-06 |  0.0006842 |   0.0003234 |     0.0011539 |    1.29e-05 |     0.0004218 |          0.0003533 |  8.20e-06 |      5.72e-05 | 0.0001745 |    5.20e-06 |      9.30e-06 |      1.00e-05 |         4.84e-05 |            0.0e+00 |      0.0002054 |    0.0000863 |     1.62e-05 | 2.1e-06 |    1.47e-05 |  4.07e-05 |    0.0000966 |  0.0000356 |     3.04e-05 |   4.4e-06 |    0.0002322 |   0.0001095 |        0.0000495 |  5.4e-06 |     0.0022865 | 1.3e-06 |  0.2653352 |       0.0001701 |         0.0002247 |         3.74e-05 | 0.0000853 |   0.0000304 |         0.0001982 |   0.0007669 |        0.0014498 |     1.29e-05 |     0.0001301 |      4.90e-05 |    0.0001294 |      0.0004409 |    4.60e-06 |    2.50e-05 |  0.0001142 | 0.0001497 |      6.4e-06 |      1.5e-06 | 5.23e-05 |      3.63e-05 | 0.0002098 |   0.0001556 |   0.0001904 |       2.14e-05 |    1.0e-06 |      0.0001554 |        0.0002415 |      0.0001373 |    0.0001670 |   5.54e-05 |      2.30e-06 |  0.0001265 |     0.0001886 | 0.0240077 |  4.40e-06 |   0.0003873 |   0.0000938 |   9.80e-06 |    3.50e-05 | 0.0000856 |   0.0002105 | 0.0000768 |    0.0373864 |          1.3e-06 |          9.00e-06 |             2.8e-06 |          3.43e-05 |      1.39e-05 |      0.0001763 |         0.0001085 |      6.7e-06 |  7.70e-06 |    0.0001175 |   0.0001028 |  0.0001987 |    0.0001417 | 0.0001647 |        1.26e-05 |       8.8e-06 |            1.47e-05 |     3.6e-06 |  0.0006641 |     8.80e-06 |   0.0093160 | 0.0001337 |  0.0007009 | 0.0000611 |      0.0013000 |     5.44e-05 |    1.29e-05 |     0.0001417 | 0.0034399 | 2.60e-05 | 0.0001118 |     0.0001392 |        0.0001142 |    2.11e-05 |     7.47e-05 |     6.78e-05 |   0.0006087 | 0.0006646 |       4.59e-05 |      0.0011864 |        7.70e-06 | 2.53e-05 |        0.0000636 |    3.30e-05 |       0.0001755 |     0.0284657 |    0.0004780 |         0.0000657 |       0.0212210 |      1.24e-05 |   4.07e-05 |      3.61e-05 |            0.0001567 |    0.0000925 |         7.06e-05 |   5.98e-05 |       0.0009068 |     0.0005893 |     0.0001281 |       0.0002147 |      0.0009277 |        0.0003136 |  0.0002525 |        5.70e-06 |         0.0e+00 |    1.19e-05 |       4.43e-05 |    5.82e-05 |         2.8e-06 |    8.80e-06 |      1.57e-05 |      5.4e-06 | 4.25e-05 |   4.1e-06 |       0.0002438 |         0.0001904 |   0.0007352 |          0.0016657 |             0.0009148 |     0.0001428 |    0.0001144 |      5.75e-05 |    0.0001263 |     2.45e-05 |           0.0000925 |    2.19e-05 |        6.44e-05 |        5.10e-05 |     4.69e-05 |       1.5e-06 |            0.0004958 |   0.0006368 |   0.0002332 |       3.6e-06 |           0.0001193 |  0.0005569 | 0.0001469 |        0.0000801 |    0.0000760 |      0.0000613 |  2.55e-05 | 0.0000902 |   1.11e-05 | 0.0008166 |     0.0001051 |   0.0001180 |   1.75e-05 |   1.47e-05 | 4.48e-05 |    1.19e-05 |     3.61e-05 | 2.1e-06 |    0.0001995 |      2.3e-06 |   4.95e-05 | 1.24e-05 |        2.80e-06 |  0.0000472 |   0.0061611 |         0.0000701 |         1.67e-05 |      1.8e-06 | 0.0005584 |   0.0001708 | 2.53e-05 |      4.4e-06 | 1.83e-05 |  4.48e-05 |       2.60e-06 |  2.94e-05 | 0.0001108 |     6.11e-05 |   0.0001958 |   0.0000356 |   0.0001737 | 0.0000636 | 1.34e-05 | 0.0003090 |    0.0003956 |           2.3e-06 |  5.13e-05 |            1.0e-06 |                                           2.65e-05 |                                       0.0003023 |                                        1.44e-05 |                                    8.20e-06 |                            0.0001495 |                                 0.0002438 |                                        2.8e-06 |                                   0.0000941 |                                      1.24e-05 |                                            6.70e-06 |                                 1.65e-05 |                                  8.50e-06 |                                                           0.0001085 |                                 0.0032562 |                                        4.79e-05 |                                  6.34e-05 |                                       0.0045420 |                                 0.0000876 |                                 0.0001551 |                                     0.0001721 |                                    0.0000673 |                                       0.0001054 |                                   0.0042367 |                                    2.80e-06 |                              0.0000531 |                               0.0000222 |                                1.80e-06 |                                     1.00e-05 |                                    1.00e-05 |                                   4.1e-06 |                                   0.0041192 |                                0.0001979 |                                  2.96e-05 |                                                 0.0000309 |                                       3.07e-05 |                                 1.44e-05 |                                8.50e-06 |                            6.49e-05 |                                    1.24e-05 |
| ShotgunWGS-ControlPig8GutMicrobiome-Day0  | 8   | Control | Day 0       | Control Day 0         |   0.0013478 |      9.90e-05 |   0.0007097 |    2.05e-05 |     0.0003269 |    0.0003021 |     0.0000709 |       0.0095018 |    4.8e-06 |       3.00e-05 |     8.19e-05 |          8.93e-05 |      0.0001283 |    0.0001319 |  0.0003576 |        0.0001186 |     0.0002488 |      0.0005859 |   0.0003412 |      6.69e-05 |  0.0000843 |      1.62e-05 | 0.0002676 |   7.9e-06 | 5.00e-06 |       0.0000990 |     0.0001931 | 1.07e-05 |     2.9e-06 |   0.0005195 |   0.0001507 |   0.0000793 |    0.0001486 |       3.55e-05 |        0.0004369 |  0.0001676 | 0.0014026 |       0.0000802 |    0.0046520 |       5.24e-05 |         4.30e-06 |   0.0000683 |      0.0002467 |   4.43e-05 | 0.0002855 |      5.26e-05 | 0.0001764 |      5.93e-05 |    0.0014702 |    0.0005607 |        0.0004626 |    0.0013019 |     0.0031158 |  1.14e-05 |     0.0004486 | 0.0001605 |       0.0002319 |     0.0001286 |  0.0001686 |   0.0001250 |    0.0003574 |     4.3e-06 |    3.00e-05 |    3.45e-05 |      7.79e-05 | 0.0101694 |     3.55e-05 | 0.0001119 |     5.95e-05 |     7.02e-05 |    7.14e-05 | 3.3e-06 | 0.0086887 |   0.1052813 |   8.05e-05 | 0.0002052 |    0.0001057 |  2.71e-05 |     3.90e-05 |   1.38e-05 |    0.0001431 |       0.0228505 |       0.0001224 |        3.95e-05 | 0.0056276 |  0.0002545 | 0.0001212 |    9.50e-06 |           3.30e-06 |       0.0001017 |   0.0015971 |      0.0002467 |     0.0005343 |      0.0000993 |      4.07e-05 | 0.0000938 | 1.12e-05 | 2.29e-05 | 0.0005540 |    0.0008445 |    0.0095013 |       4.86e-05 |        0.0021813 |            0.0028261 |        0.0001071 |   1.26e-05 |     1.57e-05 |     0.0020271 | 1.76e-05 |                 0.0000609 |                8.98e-05 |                 0.0004119 |               1.64e-05 |               3.81e-05 |               0.0004647 |                1.36e-05 |               2.31e-05 |             0.0002140 |                 4.3e-06 |                1.88e-05 |               6.05e-05 |                  6.76e-05 |                    1.67e-05 |            5.70e-06 |           2.0e-07 |             0.0004440 |          1.17e-05 |      0.0009138 |        0.0012814 |        2.12e-05 |       8.05e-05 |       0.0018825 |      8.79e-05 |  6.48e-05 |   0.0002269 |    0.0000983 |        0.0009695 |  0.0001395 |    6.20e-06 |    4.0e-06 |     0.0001088 |    0.0008223 |  3.05e-05 |      3.57e-05 |      5.19e-05 |   1.4e-06 |     0.0002633 |  0.0006652 |    0.0002636 |      0.0001545 |     3.3e-06 |       0.0001426 |        0.0001190 |        0.0001619 |       8.93e-05 |    1.12e-05 |   0.0001638 |       7.60e-06 |   0.0001467 |   6.00e-06 |   0.1009617 |     1.31e-05 |   0.0048772 | 0.0000931 | 0.0000595 |     9.81e-05 |       4.02e-05 |    1.02e-05 |     0.0005309 |   0.0063733 |         0.0001293 |        0.0001126 |       0.0009138 | 5.43e-05 |    0.0002088 |     5.71e-05 |   0.0001119 |       0.0013080 |        1.76e-05 |   0.0002993 |         2.1e-06 |   4.3e-06 | 0.0000912 |    3.3e-06 |  0.0004788 |           1.33e-05 | 0.0006181 |      8.3e-06 |     0.0002207 |     0.0001400 |       0.0008804 |       0.0000902 |   0.0003302 | 0.0000831 |     0.0002276 |    1.33e-05 |     0.0001038 |         0.0002336 |          0.0048170 |        0.0001771 |     0.0002169 |       0.0000943 |        0.0001814 |            2.95e-05 |    0.0002493 |        0.0030973 |     0.0022225 |           9.90e-05 |        9.40e-05 |        8.30e-06 |      0.0001840 |    0.0001664 |         0.0003886 | 0.0024751 |      7.62e-05 | 0.0001355 |    0.0003700 |      8.93e-05 |        4.60e-05 | 0.0001848 | 0.0149870 |   0.0006012 |    0.0000940 |   0.0030911 |  2.45e-05 |  1.33e-05 |      0.0001614 |   1.48e-05 |   1.0e-06 |         2.6e-06 |   1.00e-05 |      1.60e-05 | 0.0003017 |    0.0001817 |    0.0035749 |       1.81e-05 |                2.1e-06 |    0.0002136 |    5.48e-05 |     1.00e-05 | 0.0000893 |      0.0001038 |     0.0001074 |   0.0010485 |      0.0059181 |   0.0601988 |       0.0006819 |        0.0296095 |  0.0000488 |    2.57e-05 |    2.71e-05 |        0.0002305 |   0.0016578 |  0.0002571 |       2.79e-05 |  0.0007926 |      0.0017283 |   0.0001786 | 0.0003407 |    9.30e-06 |     0.0025525 |    4.67e-05 |         2.0e-07 |   0.0012178 | 5.79e-05 | 4.26e-05 |    0.0000809 |   0.0021025 | 0.0016695 |         6.55e-05 | 1.40e-05 |   4.14e-05 |   0.0001609 |          5.26e-05 |      6.52e-05 | 5.88e-05 |    7.0e-07 | 0.0005845 |      5.36e-05 |      0.0001469 |    7.4e-06 |   0.0003728 | 0.0001350 |        1.45e-05 |     0.0002857 |  0.0000921 |   2.71e-05 |      2.14e-05 |  1.43e-05 |        2.29e-05 |      1.57e-05 |  4.10e-05 |      2.60e-05 |    2.14e-05 |       7.33e-05 |   1.19e-05 |      2.24e-05 |       0.0006481 |         3.71e-05 |    0.0008431 |      0.0014135 |       4.31e-05 |     0.0000831 |     0.0001567 | 2.79e-05 |   0.0001729 | 1.07e-05 |  0.0043720 |       2.29e-05 |        4.02e-05 |         7.74e-05 |     1.10e-05 |       3.17e-05 |   6.07e-05 |  0.0001157 |   1.43e-05 |     5.7e-06 |  0.0007131 |        3.55e-05 |   5.31e-05 |   4.81e-05 |         0.0000883 | 0.0001131 |   0.0000890 |  3.86e-05 |           2.88e-05 |   0.0001488 | 1.64e-05 |  0.0002228 |      1.71e-05 | 0.0000693 | 7.19e-05 | 0.0001533 |  8.24e-05 |      5.76e-05 |   5.71e-05 |         0.0e+00 |  3.02e-05 |  1.0e-05 |   7.9e-06 |     0.0365643 |   0.0006574 |           0.0000248 |   4.60e-05 | 0.0001417 |      0.0004093 |        0.0004226 |  0.0001619 | 0.0000962 |   3.05e-05 |     6.48e-05 |     2e-07 |  0.0001926 |        2.0e-07 |  0.0000995 |    0.0007257 |   0.0004550 |    1.76e-05 | 0.0013728 | 3.30e-06 |      6.0e-06 |   1.90e-05 | 2.57e-05 | 2.69e-05 |      0.0003719 |   0.0001431 |    2.17e-05 |     0.0001840 |        0.0001967 |   8.30e-06 |  0.0000921 |  0.0003433 |   6.14e-05 |    0.0001674 |   0.0001752 |      2.02e-05 |       3.14e-05 | 0.0002436 |   0.0059802 |   0.0002057 |  0.0000633 |     0.0001245 |       1.36e-05 |          0.0006671 |          0.0001657 |     5.26e-05 |        0.0001343 |     0.0004228 |          0.0004162 |      0.0001257 |        3.69e-05 |          5.52e-05 |     0.0001081 |     4.19e-05 |      8.33e-05 |    0.0000955 |      0.0005721 |      0.0001462 |         5.29e-05 |        0.0001212 |           0.0001283 |            6.40e-06 |       2.33e-05 |          4.93e-05 |   0.0001126 |        8.43e-05 |      2.00e-05 |        0.0002076 |     4.10e-05 |     0.0001271 |     1.17e-05 |     2.24e-05 |      5.50e-05 |     3.81e-05 |    7.6e-06 |    6.79e-05 |    3.38e-05 |    9.76e-05 |   1.76e-05 |       7.88e-05 |   0.0001602 |   0.0218968 |  0.0024935 |        5.5e-06 | 2.10e-05 | 0.0011928 |  2.64e-05 | 0.0000190 |        0.0002362 |     0.0005878 |  0.0004047 |  0.0002167 |         3.6e-06 |  2.67e-05 |    7.76e-05 |     1.93e-05 |      1.2e-06 |      0.0005021 |  1.74e-05 |     2.36e-05 | 5.26e-05 | 8.3e-06 | 0.0003097 |      1.50e-05 |    4.31e-05 |       1.79e-05 |   2.31e-05 |     0.0001198 |   0.0001093 |    4.24e-05 |     0.0001821 |    0.0001531 |       1.29e-05 |     8.31e-05 |   7.31e-05 | 0.0001100 |    0.0002026 |     6.57e-05 |  1.98e-05 | 0.0002267 |       0.0000931 |      7.1e-06 |     3.21e-05 |   3.79e-05 |      5.79e-05 |      0.0005259 |    0.0001162 |       1.64e-05 |   7.0e-07 |  0.0002638 |    2.43e-05 | 0.0072738 | 0.0003702 |    0.0017714 | 1.93e-05 |     2.26e-05 |       3.40e-05 |     3.83e-05 |   0.0001764 |       0.0000005 |         3.8e-06 |          4.5e-06 |     0.0027108 |    0.0021187 | 0.0001431 |       0.0093044 |      1.12e-05 |          1.7e-06 |   8.05e-05 |   1.69e-05 |     0.0002426 |     8.71e-05 |     3.36e-05 |   0.0001257 |    3.6e-06 |      0.0002074 |   0.0005376 |  0.0013526 |   8.30e-06 |  0.0008312 |   0.0003002 |     0.0013007 |    1.48e-05 |     0.0006295 |          0.0004781 |  1.02e-05 |      6.57e-05 | 0.0002369 |    6.20e-06 |      1.98e-05 |      1.29e-05 |         4.12e-05 |            1.9e-06 |      0.0002690 |    0.0001145 |     3.71e-05 | 5.5e-06 |    2.60e-05 |  6.90e-05 |    0.0001586 |  0.0000562 |     3.86e-05 |   9.5e-06 |    0.0003917 |   0.0002409 |        0.0000940 |  5.5e-06 |     0.0039265 | 1.7e-06 |  0.1953868 |       0.0002638 |         0.0002686 |         6.26e-05 | 0.0001017 |   0.0000571 |         0.0002417 |   0.0009345 |        0.0005771 |     1.76e-05 |     0.0001795 |      7.81e-05 |    0.0001679 |      0.0007288 |    6.20e-06 |    3.19e-05 |  0.0001490 | 0.0001812 |      5.7e-06 |      5.0e-07 | 8.31e-05 |      4.86e-05 | 0.0002498 |   0.0001893 |   0.0002359 |       3.60e-05 |    7.0e-07 |      0.0002340 |        0.0003250 |      0.0001838 |    0.0002314 |   7.29e-05 |      2.40e-06 |  0.0002000 |     0.0002702 | 0.0171608 |  9.30e-06 |   0.0005097 |   0.0001350 |   2.52e-05 |    5.79e-05 | 0.0000945 |   0.0002457 | 0.0001267 |    0.0555783 |          0.0e+00 |          1.00e-05 |             2.6e-06 |          4.86e-05 |      2.57e-05 |      0.0002005 |         0.0001519 |      4.3e-06 |  1.26e-05 |    0.0001588 |   0.0001309 |  0.0002948 |    0.0001664 | 0.0002036 |        2.07e-05 |       9.5e-06 |            3.90e-05 |     2.9e-06 |  0.0008007 |     1.62e-05 |   0.0080618 | 0.0001619 |  0.0009002 | 0.0001214 |      0.0010104 |     4.19e-05 |    2.19e-05 |     0.0002159 | 0.0041267 | 4.00e-05 | 0.0001348 |     0.0002038 |        0.0002421 |    4.48e-05 |     8.76e-05 |     7.76e-05 |   0.0009776 | 0.0008926 |       6.52e-05 |      0.0013411 |        1.38e-05 | 4.24e-05 |        0.0001002 |    4.90e-05 |       0.0002655 |     0.0077406 |    0.0007007 |         0.0001212 |       0.0155818 |      2.43e-05 |   5.05e-05 |      4.90e-05 |            0.0001390 |    0.0001029 |         8.17e-05 |   9.93e-05 |       0.0010376 |     0.0010085 |     0.0001517 |       0.0002898 |      0.0011535 |        0.0003964 |  0.0003407 |        2.02e-05 |         8.3e-06 |    6.40e-06 |       5.93e-05 |    9.86e-05 |         4.5e-06 |    1.33e-05 |      2.33e-05 |      6.9e-06 | 5.79e-05 |   4.5e-06 |       0.0003150 |         0.0002264 |   0.0008414 |          0.0019737 |             0.0009819 |     0.0001759 |    0.0001381 |      8.05e-05 |    0.0001938 |     3.10e-05 |           0.0001164 |    2.88e-05 |        9.24e-05 |        7.69e-05 |     7.14e-05 |       1.7e-06 |            0.0006600 |   0.0011904 |   0.0002948 |       2.9e-06 |           0.0001426 |  0.0007307 | 0.0001933 |        0.0001059 |    0.0001067 |      0.0000852 |  3.17e-05 | 0.0001143 |   9.80e-06 | 0.0012692 |     0.0001229 |   0.0003755 |   1.62e-05 |   2.14e-05 | 6.67e-05 |    2.88e-05 |     3.50e-05 | 6.4e-06 |    0.0002702 |      5.5e-06 |   5.74e-05 | 2.33e-05 |        5.50e-06 |  0.0000769 |   0.0047134 |         0.0001209 |         7.21e-05 |      5.7e-06 | 0.0007326 |   0.0003717 | 3.95e-05 |      7.1e-06 | 3.62e-05 |  5.05e-05 |       3.60e-06 |  4.76e-05 | 0.0001538 |     7.59e-05 |   0.0002807 |   0.0000402 |   0.0001895 | 0.0000779 | 2.60e-05 | 0.0003219 |    0.0004664 |           6.4e-06 |  6.40e-05 |            1.4e-06 |                                           4.98e-05 |                                       0.0003367 |                                        3.36e-05 |                                    1.05e-05 |                            0.0001955 |                                 0.0005247 |                                        4.0e-06 |                                   0.0001940 |                                      1.71e-05 |                                            1.17e-05 |                                 4.07e-05 |                                  1.12e-05 |                                                           0.0001238 |                                 0.0040044 |                                        8.45e-05 |                                  7.83e-05 |                                       0.0067576 |                                 0.0001309 |                                 0.0002533 |                                     0.0002098 |                                    0.0001726 |                                       0.0001759 |                                   0.0050081 |                                    3.80e-06 |                              0.0001398 |                               0.0001421 |                                7.90e-06 |                                     2.02e-05 |                                    2.26e-05 |                                   2.9e-06 |                                   0.0122500 |                                0.0002545 |                                  4.55e-05 |                                                 0.0000800 |                                       6.05e-05 |                                 1.90e-05 |                                7.90e-06 |                            7.40e-05 |                                    7.90e-06 |
| ShotgunWGS-ControlPig3GutMicrobiome-Day14 | 3   | Control | Day 14      | Control Day 14        |   0.0010663 |      6.91e-05 |   0.0005364 |    1.55e-05 |     0.0002629 |    0.0002116 |     0.0000510 |       0.0082861 |    6.5e-06 |       2.10e-05 |     6.29e-05 |          5.75e-05 |      0.0000800 |    0.0000780 |  0.0002916 |        0.0000635 |     0.0001691 |      0.0004657 |   0.0001808 |      4.69e-05 |  0.0001909 |      4.70e-06 | 0.0001973 |   6.7e-06 | 3.10e-06 |       0.0000748 |     0.0001147 | 3.40e-06 |     2.8e-06 |   0.0003675 |   0.0001202 |   0.0000782 |    0.0000704 |       3.65e-05 |        0.0003546 |  0.0001259 | 0.0007555 |       0.0000842 |    0.0034668 |       4.48e-05 |         1.30e-06 |   0.0000505 |      0.0001572 |   2.90e-05 | 0.0002165 |      5.05e-05 | 0.0001256 |      3.81e-05 |    0.0010849 |    0.0003916 |        0.0003149 |    0.0010981 |     0.0022428 |  1.14e-05 |     0.0003667 | 0.0001215 |       0.0001254 |     0.0001054 |  0.0001474 |   0.0000995 |    0.0002015 |     3.6e-06 |    2.30e-05 |    2.87e-05 |      5.52e-05 | 0.0051601 |     2.18e-05 | 0.0000875 |     4.77e-05 |     6.16e-05 |    3.96e-05 | 1.6e-06 | 0.0067811 |   0.0763643 |   4.40e-05 | 0.0001453 |    0.0000816 |  1.29e-05 |     2.49e-05 |   9.80e-06 |    0.0000575 |       0.0065677 |       0.0000738 |        2.33e-05 | 0.0048146 |  0.0002178 | 0.0000637 |    5.40e-06 |           3.10e-06 |       0.0000578 |   0.0012535 |      0.0001707 |     0.0004266 |      0.0000531 |      2.98e-05 | 0.0000583 | 7.80e-06 | 1.97e-05 | 0.0003965 |    0.0006247 |    0.0089022 |       2.69e-05 |        0.0017150 |            0.0024133 |        0.0000720 |   5.40e-06 |     9.60e-06 |     0.0012351 | 1.42e-05 |                 0.0000544 |                5.59e-05 |                 0.0002046 |               8.00e-06 |               2.25e-05 |               0.0003934 |                1.01e-05 |               1.35e-05 |             0.0001797 |                 2.6e-06 |                2.07e-05 |               3.88e-05 |                  4.74e-05 |                    1.68e-05 |            2.80e-06 |           5.0e-07 |             0.0003660 |          8.00e-06 |      0.0006288 |        0.0010722 |        1.79e-05 |       6.99e-05 |       0.0009637 |      5.39e-05 |  6.37e-05 |   0.0002082 |    0.0000526 |        0.0007143 |  0.0001274 |    5.40e-06 |    8.0e-06 |     0.0000800 |    0.0005915 |  1.94e-05 |      1.89e-05 |      4.22e-05 |   8.0e-07 |     0.0001784 |  0.0005144 |    0.0002326 |      0.0001580 |     1.8e-06 |       0.0000969 |        0.0000870 |        0.0000873 |       4.35e-05 |    8.00e-06 |   0.0001217 |       7.00e-06 |   0.0000868 |   1.01e-05 |   0.0774401 |     4.70e-06 |   0.0053002 | 0.0000837 | 0.0000438 |     7.59e-05 |       3.03e-05 |    3.60e-06 |     0.0003274 |   0.0064006 |         0.0000922 |        0.0000640 |       0.0004708 | 3.76e-05 |    0.0001127 |     4.87e-05 |   0.0000834 |       0.0008477 |        7.50e-06 |   0.0002261 |         1.6e-06 |   4.1e-06 | 0.0001197 |    2.6e-06 |  0.0004048 |           1.22e-05 | 0.0004131 |      3.9e-06 |     0.0000956 |     0.0000987 |       0.0008091 |       0.0000811 |   0.0002502 | 0.0000637 |     0.0001968 |    7.30e-06 |     0.0000774 |         0.0001564 |          0.0038092 |        0.0001290 |     0.0001515 |       0.0000875 |        0.0001367 |            1.84e-05 |    0.0002129 |        0.0026608 |     0.0014037 |           7.02e-05 |        6.14e-05 |        5.70e-06 |      0.0001525 |    0.0001274 |         0.0002318 | 0.0019038 |      6.32e-05 | 0.0001140 |    0.0002872 |      3.96e-05 |        3.06e-05 | 0.0001072 | 0.0083552 |   0.0004289 |    0.0000616 |   0.0024832 |  1.53e-05 |  9.30e-06 |      0.0001171 |   1.01e-05 |   1.0e-06 |         5.0e-07 |   3.60e-06 |      1.37e-05 | 0.0000293 |    0.0001238 |    0.0034243 |       1.40e-05 |                1.0e-06 |    0.0001512 |    2.67e-05 |     8.50e-06 | 0.0000790 |      0.0000541 |     0.0000800 |   0.0007941 |      0.0051174 |   0.0696863 |       0.0005392 |        0.0454851 |  0.0000368 |    2.07e-05 |    1.53e-05 |        0.0002056 |   0.0012911 |  0.0001608 |       1.92e-05 |  0.0006146 |      0.0012753 |   0.0001318 | 0.0002471 |    6.70e-06 |     0.0017461 |    4.17e-05 |         5.0e-07 |   0.0004649 | 4.45e-05 | 2.07e-05 |    0.0000624 |   0.0016751 | 0.0013395 |         4.58e-05 | 6.20e-06 |   2.69e-05 |   0.0001388 |          3.55e-05 |      4.40e-05 | 2.67e-05 |    1.3e-06 | 0.0003551 |      4.56e-05 |      0.0001699 |    6.5e-06 |   0.0002611 | 0.0000956 |        1.14e-05 |     0.0002357 |  0.0000710 |   2.23e-05 |      1.19e-05 |  9.60e-06 |        1.48e-05 |      9.10e-06 |  2.67e-05 |      2.02e-05 |    1.06e-05 |       5.23e-05 |   1.27e-05 |      1.94e-05 |       0.0005309 |         2.77e-05 |    0.0005602 |      0.0013338 |       5.34e-05 |     0.0000699 |     0.0001352 | 2.54e-05 |   0.0001479 | 7.00e-06 |  0.0032586 |       2.28e-05 |        2.56e-05 |         4.71e-05 |     9.10e-06 |       1.97e-05 |   5.13e-05 |  0.0000844 |   9.10e-06 |     4.1e-06 |  0.0005457 |        2.49e-05 |   3.03e-05 |   3.70e-05 |         0.0000767 | 0.0000629 |   0.0000546 |  2.72e-05 |           2.15e-05 |   0.0000699 | 1.14e-05 |  0.0001554 |      9.10e-06 | 0.0000401 | 4.33e-05 | 0.0000963 |  4.79e-05 |      3.24e-05 |   3.16e-05 |         3.0e-07 |  1.58e-05 |  3.6e-06 |   4.1e-06 |     0.0424407 |   0.0011250 |           0.0000127 |   4.14e-05 | 0.0000963 |      0.0003903 |        0.0002950 |  0.0001254 | 0.0000484 |   1.58e-05 |     3.99e-05 |     3e-07 |  0.0001375 |        8.0e-07 |  0.0000697 |    0.0006397 |   0.0005550 |    1.58e-05 | 0.0011362 | 3.10e-06 |      4.4e-06 |   1.61e-05 | 1.63e-05 | 2.43e-05 |      0.0003113 |   0.0001347 |    1.68e-05 |     0.0001375 |        0.0001279 |   5.20e-06 |  0.0000699 |  0.0001976 |   3.91e-05 |    0.0001217 |   0.0001199 |      1.55e-05 |       2.90e-05 | 0.0001404 |   0.0070940 |   0.0001349 |  0.0000502 |     0.0000805 |       8.80e-06 |          0.0005001 |          0.0001083 |     3.78e-05 |        0.0001152 |     0.0003175 |          0.0002989 |      0.0001114 |        2.67e-05 |          3.03e-05 |     0.0000831 |     3.13e-05 |      6.40e-05 |    0.0000694 |      0.0003955 |      0.0000984 |         4.43e-05 |        0.0000966 |           0.0000888 |            9.80e-06 |       1.40e-05 |          3.57e-05 |   0.0000736 |        5.90e-05 |      1.58e-05 |        0.0001634 |     4.20e-05 |     0.0000974 |     8.30e-06 |     1.99e-05 |      3.52e-05 |     2.75e-05 |    6.0e-06 |    4.58e-05 |    1.89e-05 |    9.74e-05 |   1.45e-05 |       5.75e-05 |   0.0000901 |   0.0148157 |  0.0003079 |        3.1e-06 | 1.42e-05 | 0.0010292 |  2.05e-05 | 0.0000104 |        0.0001142 |     0.0003781 |  0.0003056 |  0.0001492 |         5.0e-07 |  9.60e-06 |    4.97e-05 |     1.50e-05 |      1.0e-06 |      0.0003986 |  1.32e-05 |     1.76e-05 | 4.12e-05 | 5.4e-06 | 0.0002108 |      8.30e-06 |    2.93e-05 |       9.30e-06 |   1.24e-05 |     0.0000769 |   0.0000919 |    3.26e-05 |     0.0001248 |    0.0001051 |       9.80e-06 |     6.24e-05 |   5.49e-05 | 0.0000704 |    0.0001028 |     4.09e-05 |  1.63e-05 | 0.0001816 |       0.0000948 |      5.2e-06 |     1.76e-05 |   2.95e-05 |      3.96e-05 |      0.0004079 |    0.0000570 |       1.35e-05 |   1.8e-06 |  0.0002269 |    1.92e-05 | 0.0035111 | 0.0003020 |    0.0020639 | 2.02e-05 |     1.86e-05 |       2.49e-05 |     2.90e-05 |   0.0001507 |       0.0000018 |         8.0e-07 |          3.0e-07 |     0.0020722 |    0.0013423 | 0.0001046 |       0.0061947 |      5.70e-06 |          8.0e-07 |   5.26e-05 |   1.04e-05 |     0.0000956 |     6.68e-05 |     2.95e-05 |   0.0000710 |    7.8e-06 |      0.0001577 |   0.0004610 |  0.0008676 |   5.40e-06 |  0.0007184 |   0.0002629 |     0.0010857 |    8.50e-06 |     0.0003854 |          0.0003142 |  4.90e-06 |      5.67e-05 | 0.0001709 |    4.70e-06 |      1.27e-05 |      5.70e-06 |         4.69e-05 |            3.0e-07 |      0.0002201 |    0.0000746 |     1.42e-05 | 3.4e-06 |    2.15e-05 |  5.93e-05 |    0.0001028 |  0.0000287 |     2.54e-05 |   5.2e-06 |    0.0002346 |   0.0002020 |        0.0000702 |  8.3e-06 |     0.0020786 | 1.3e-06 |  0.2305057 |       0.0002497 |         0.0001601 |         4.38e-05 | 0.0000710 |   0.0000381 |         0.0001836 |   0.0006775 |        0.0005781 |     1.04e-05 |     0.0001204 |      4.38e-05 |    0.0001347 |      0.0003996 |    2.30e-06 |    3.08e-05 |  0.0001083 | 0.0001492 |      5.2e-06 |      5.0e-07 | 5.44e-05 |      2.46e-05 | 0.0001771 |   0.0001494 |   0.0001683 |       2.28e-05 |    1.8e-06 |      0.0001500 |        0.0002357 |      0.0001533 |    0.0001634 |   5.78e-05 |      1.80e-06 |  0.0001173 |     0.0002080 | 0.0215569 |  6.20e-06 |   0.0004234 |   0.0000901 |   1.61e-05 |    3.26e-05 | 0.0000883 |   0.0001805 | 0.0000976 |    0.0341373 |          1.3e-06 |          8.50e-06 |             3.4e-06 |          2.90e-05 |      1.86e-05 |      0.0002256 |         0.0000862 |      4.1e-06 |  7.30e-06 |    0.0001202 |   0.0000982 |  0.0001968 |    0.0001013 | 0.0000974 |        8.50e-06 |       8.0e-06 |            1.61e-05 |     3.4e-06 |  0.0005734 |     1.09e-05 |   0.0096929 | 0.0001329 |  0.0006918 | 0.0001075 |      0.0009518 |     3.60e-05 |    2.93e-05 |     0.0001518 | 0.0028976 | 3.16e-05 | 0.0001096 |     0.0001500 |        0.0001334 |    2.80e-05 |     8.21e-05 |     6.45e-05 |   0.0005936 | 0.0006374 |       3.65e-05 |      0.0011680 |        1.32e-05 | 2.56e-05 |        0.0000557 |    3.52e-05 |       0.0002131 |     0.1117937 |    0.0004628 |         0.0000761 |       0.0196324 |      1.97e-05 |   2.87e-05 |      3.65e-05 |            0.0001295 |    0.0000914 |         6.68e-05 |   6.89e-05 |       0.0008067 |     0.0010129 |     0.0001285 |       0.0001942 |      0.0009562 |        0.0003017 |  0.0002357 |        9.30e-06 |         5.2e-06 |    2.30e-06 |       3.96e-05 |    6.40e-05 |         6.2e-06 |    8.30e-06 |      1.89e-05 |      3.4e-06 | 4.20e-05 |   7.8e-06 |       0.0002497 |         0.0001639 |   0.0006892 |          0.0015045 |             0.0007671 |     0.0001254 |    0.0000919 |      4.58e-05 |    0.0001378 |     2.20e-05 |           0.0000901 |    2.41e-05 |        7.87e-05 |        5.65e-05 |     4.53e-05 |       2.1e-06 |            0.0005076 |   0.0009085 |   0.0002375 |       6.5e-06 |           0.0001228 |  0.0005633 | 0.0001331 |        0.0000697 |    0.0000839 |      0.0000647 |  2.95e-05 | 0.0000738 |   3.60e-06 | 0.0007803 |     0.0000927 |   0.0001207 |   8.00e-06 |   1.04e-05 | 4.33e-05 |    1.68e-05 |     2.41e-05 | 4.9e-06 |    0.0001660 |      2.6e-06 |   3.94e-05 | 7.30e-06 |        2.30e-06 |  0.0000756 |   0.0044696 |         0.0000826 |         3.63e-05 |      3.4e-06 | 0.0005467 |   0.0001606 | 2.82e-05 |      4.4e-06 | 2.87e-05 |  5.67e-05 |       3.90e-06 |  3.34e-05 | 0.0000953 |     6.24e-05 |   0.0001994 |   0.0000334 |   0.0001158 | 0.0000552 | 1.09e-05 | 0.0002763 |    0.0003709 |           3.1e-06 |  6.37e-05 |            1.8e-06 |                                           2.69e-05 |                                       0.0002730 |                                        2.10e-05 |                                    7.50e-06 |                            0.0001494 |                                 0.0002119 |                                        1.3e-06 |                                   0.0001740 |                                      1.06e-05 |                                            5.40e-06 |                                 1.63e-05 |                                  9.60e-06 |                                                           0.0000826 |                                 0.0029442 |                                        5.36e-05 |                                  6.40e-05 |                                       0.0049519 |                                 0.0000805 |                                 0.0001564 |                                     0.0001533 |                                    0.0000894 |                                       0.0001339 |                                   0.0043629 |                                    4.10e-06 |                              0.0000554 |                               0.0000422 |                                4.90e-06 |                                     1.63e-05 |                                    1.24e-05 |                                   2.3e-06 |                                   0.0055364 |                                0.0002965 |                                  3.47e-05 |                                                 0.0000412 |                                       3.57e-05 |                                 8.80e-06 |                                4.10e-06 |                            4.71e-05 |                                    1.61e-05 |
| ShotgunWGS-TomatoPig14GutMicrobiome-Day7  | 14  | Tomato  | Day 7       | Tomato Day 7          |   0.0013116 |      1.09e-04 |   0.0008422 |    3.41e-05 |     0.0003321 |    0.0003254 |     0.0000749 |       0.0094690 |    5.8e-06 |       2.25e-05 |     8.32e-05 |          8.49e-05 |      0.0001123 |    0.0000940 |  0.0002788 |        0.0001623 |     0.0003678 |      0.0008322 |   0.0003370 |      7.57e-05 |  0.0001123 |      6.70e-06 | 0.0010511 |   8.3e-06 | 1.25e-05 |       0.0002189 |     0.0001914 | 1.08e-05 |     5.0e-06 |   0.0004810 |   0.0001298 |   0.0001523 |    0.0000866 |       2.91e-05 |        0.0004619 |  0.0003254 | 0.0011734 |       0.0001099 |    0.0044083 |       9.15e-05 |         2.00e-05 |   0.0001639 |      0.0002014 |   4.58e-05 | 0.0002314 |      4.33e-05 | 0.0001872 |      7.32e-05 |    0.0013116 |    0.0006342 |        0.0003803 |    0.0016178 |     0.0029386 |  1.50e-05 |     0.0004427 | 0.0001415 |       0.0003687 |     0.0001806 |  0.0002022 |   0.0001290 |    0.0003379 |     3.3e-06 |    3.66e-05 |    3.58e-05 |      5.16e-05 | 0.0029136 |     3.16e-05 | 0.0001323 |     5.58e-05 |     6.57e-05 |    8.99e-05 | 2.5e-06 | 0.0083622 |   0.0892998 |   6.57e-05 | 0.0002755 |    0.0001232 |  4.16e-05 |     3.08e-05 |   5.08e-05 |    0.0001689 |       0.0046122 |       0.0000874 |        3.33e-05 | 0.0059104 |  0.0002480 | 0.0001773 |    2.41e-05 |           3.30e-06 |       0.0001423 |   0.0015346 |      0.0001997 |     0.0004985 |      0.0001007 |      3.83e-05 | 0.0001332 | 4.49e-05 | 9.90e-05 | 0.0005776 |    0.0008663 |    0.0090321 |       5.49e-05 |        0.0021047 |            0.0029727 |        0.0001015 |   5.80e-06 |     1.75e-05 |     0.0017227 | 1.75e-05 |                 0.0000699 |                6.99e-05 |                 0.0002929 |               2.66e-05 |               3.00e-05 |               0.0004402 |                3.58e-05 |               3.00e-05 |             0.0001997 |                 7.5e-06 |                2.66e-05 |               8.07e-05 |                  4.91e-05 |                    1.75e-05 |            2.83e-05 |           5.0e-06 |             0.0003537 |          9.20e-06 |      0.0007041 |        0.0012392 |        3.41e-05 |       8.41e-05 |       0.0011451 |      8.66e-05 |  6.99e-05 |   0.0001872 |    0.0001074 |        0.0009454 |  0.0001598 |    1.25e-05 |    5.0e-06 |     0.0000741 |    0.0006317 |  2.33e-05 |      3.00e-05 |      4.16e-05 |   8.0e-07 |     0.0002122 |  0.0006133 |    0.0002563 |      0.0001798 |     3.3e-06 |       0.0001481 |        0.0001440 |        0.0001140 |       2.75e-05 |    1.00e-05 |   0.0003412 |       4.20e-06 |   0.0001706 |   5.80e-06 |   0.0910882 |     7.50e-06 |   0.0026098 | 0.0001781 | 0.0000516 |     9.40e-05 |       4.58e-05 |    1.17e-05 |     0.0003920 |   0.0047045 |         0.0001182 |        0.0000682 |       0.0026390 | 8.41e-05 |    0.0001531 |     7.99e-05 |   0.0002230 |       0.0007665 |        1.58e-05 |   0.0002663 |         2.5e-06 |   4.2e-06 | 0.0000108 |    2.5e-06 |  0.0004827 |           1.91e-05 | 0.0004752 |      5.0e-06 |     0.0001606 |     0.0001357 |       0.0009454 |       0.0001298 |   0.0003088 | 0.0000674 |     0.0002147 |    2.00e-05 |     0.0001123 |         0.0002380 |          0.0042968 |        0.0001664 |     0.0002114 |       0.0001148 |        0.0001723 |            3.25e-05 |    0.0002530 |        0.0028811 |     0.0032182 |           9.49e-05 |        9.32e-05 |        1.08e-05 |      0.0001997 |    0.0001823 |         0.0003262 | 0.0011618 |      9.15e-05 | 0.0002405 |    0.0003529 |      7.82e-05 |        4.91e-05 | 0.0001806 | 0.0068667 |   0.0004660 |    0.0002139 |   0.0023435 |  3.16e-05 |  1.91e-05 |      0.0002056 |   2.33e-05 |   4.2e-06 |         5.0e-06 |   1.08e-05 |      3.50e-05 | 0.0000399 |    0.0003728 |    0.0036326 |       1.41e-05 |                2.5e-06 |    0.0002413 |    5.91e-05 |     1.41e-05 | 0.0001748 |      0.0001015 |     0.0000924 |   0.0100665 |      0.0060960 |   0.0458013 |       0.0005767 |        0.0398517 |  0.0001631 |    2.91e-05 |    3.66e-05 |        0.0002122 |   0.0014564 |  0.0002538 |       2.83e-05 |  0.0011035 |      0.0014522 |   0.0002372 | 0.0002746 |    9.20e-06 |     0.0024900 |    5.41e-05 |         0.0e+00 |   0.0004253 | 9.82e-05 | 2.58e-05 |    0.0000508 |   0.0020597 | 0.0016245 |         7.91e-05 | 5.80e-06 |   3.41e-05 |   0.0001465 |          5.49e-05 |      4.83e-05 | 5.74e-05 |    1.7e-06 | 0.0005176 |      6.66e-05 |      0.0001581 |    7.5e-06 |   0.0007066 | 0.0001806 |        2.00e-05 |     0.0002530 |  0.0001115 |   3.50e-05 |      2.16e-05 |  1.91e-05 |        2.33e-05 |      1.17e-05 |  5.66e-05 |      2.41e-05 |    2.66e-05 |       9.15e-05 |   1.83e-05 |      2.91e-05 |       0.0006084 |         6.49e-05 |    0.0029286 |      0.0014098 |       7.74e-05 |     0.0001015 |     0.0001664 | 2.75e-05 |   0.0002971 | 1.83e-05 |  0.0036768 |       3.74e-05 |        4.74e-05 |         6.24e-05 |     1.50e-05 |       2.08e-05 |   5.66e-05 |  0.0002505 |   2.41e-05 |     1.0e-05 |  0.0005934 |        4.33e-05 |   6.32e-05 |   4.91e-05 |         0.0001148 | 0.0001132 |   0.0001315 |  6.66e-05 |           2.41e-05 |   0.0001881 | 2.16e-05 |  0.0005077 |      1.58e-05 | 0.0001165 | 8.57e-05 | 0.0001523 |  8.82e-05 |      6.91e-05 |   7.49e-05 |         0.0e+00 |  4.99e-05 |  8.3e-06 |   6.7e-06 |     0.0870628 |   0.0006733 |           0.0001173 |   7.49e-05 | 0.0001997 |      0.0003620 |        0.0004244 |  0.0002155 | 0.0000866 |   3.08e-05 |     7.07e-05 |     8e-07 |  0.0001972 |        0.0e+00 |  0.0000882 |    0.0006200 |   0.0004536 |    3.41e-05 | 0.0013740 | 1.50e-05 |      5.0e-06 |   2.08e-05 | 3.00e-05 | 2.00e-05 |      0.0003878 |   0.0001714 |    2.41e-05 |     0.0001756 |        0.0001872 |   1.33e-05 |  0.0001556 |  0.0003370 |   5.58e-05 |    0.0002480 |   0.0002314 |      3.00e-05 |       3.00e-05 | 0.0001956 |   0.0051190 |   0.0001906 |  0.0001074 |     0.0000915 |       1.91e-05 |          0.0022853 |          0.0002272 |     7.24e-05 |        0.0001490 |     0.0004469 |          0.0004003 |      0.0001273 |        4.49e-05 |          6.74e-05 |     0.0001015 |     6.57e-05 |      9.99e-05 |    0.0001032 |      0.0006366 |      0.0001631 |         7.49e-05 |        0.0001381 |           0.0003004 |            1.91e-05 |       4.83e-05 |          5.33e-05 |   0.0000816 |        8.90e-05 |      3.66e-05 |        0.0001947 |     3.99e-05 |     0.0001515 |     2.33e-05 |     1.83e-05 |      6.91e-05 |     5.66e-05 |    8.3e-06 |    5.83e-05 |    2.25e-05 |    9.65e-05 |   2.33e-05 |       6.99e-05 |   0.0001273 |   0.0064530 |  0.0042194 |        4.2e-06 | 3.83e-05 | 0.0011343 |  4.83e-05 | 0.0001132 |        0.0001848 |     0.0005975 |  0.0004952 |  0.0001947 |         1.7e-06 |  2.25e-05 |    8.07e-05 |     3.25e-05 |      5.0e-06 |      0.0004985 |  1.08e-05 |     1.75e-05 | 5.91e-05 | 1.0e-05 | 0.0004253 |      1.58e-05 |    4.58e-05 |       3.74e-05 |   2.50e-05 |     0.0001090 |   0.0001165 |    4.08e-05 |     0.0002130 |    0.0001989 |       1.83e-05 |     9.65e-05 |   7.66e-05 | 0.0001173 |    0.0001823 |     5.58e-05 |  2.00e-05 | 0.0002139 |       0.0000666 |      1.0e-05 |     2.58e-05 |   3.16e-05 |      5.99e-05 |      0.0004977 |    0.0000782 |       2.16e-05 |   3.3e-06 |  0.0002380 |    3.16e-05 | 0.0019391 | 0.0002455 |    0.0028795 | 1.83e-05 |     3.66e-05 |       2.41e-05 |     4.49e-05 |   0.0003146 |       0.0001049 |         5.0e-06 |          4.2e-06 |     0.0024667 |    0.0015330 | 0.0002413 |       0.0087625 |      1.58e-05 |          8.0e-07 |   6.49e-05 |   7.50e-06 |     0.0000641 |     8.99e-05 |     4.24e-05 |   0.0002538 |    5.8e-06 |      0.0003071 |   0.0005118 |  0.0010345 |   8.30e-06 |  0.0007898 |   0.0002522 |     0.0010827 |    1.00e-05 |     0.0007490 |          0.0005110 |  1.41e-05 |      5.66e-05 | 0.0002538 |    1.08e-05 |      2.66e-05 |      1.33e-05 |         4.41e-05 |            8.0e-07 |      0.0005709 |    0.0002413 |     3.16e-05 | 6.7e-06 |    4.16e-05 |  4.49e-05 |    0.0001257 |  0.0001057 |     4.74e-05 |   5.0e-06 |    0.0002630 |   0.0001748 |        0.0000832 |  7.5e-06 |     0.0028529 | 8.0e-07 |  0.2312992 |       0.0002297 |         0.0004169 |         4.49e-05 | 0.0002430 |   0.0001323 |         0.0005426 |   0.0012242 |        0.0005243 |     1.58e-05 |     0.0002189 |      6.57e-05 |    0.0003961 |      0.0005368 |    1.08e-05 |    4.41e-05 |  0.0001806 | 0.0001997 |      5.8e-06 |      3.3e-06 | 9.40e-05 |      5.74e-05 | 0.0002355 |   0.0001573 |   0.0002480 |       3.74e-05 |    8.0e-07 |      0.0001656 |        0.0003645 |      0.0002372 |    0.0002056 |   9.40e-05 |      1.66e-05 |  0.0001456 |     0.0002646 | 0.0101514 |  3.30e-06 |   0.0004635 |   0.0001573 |   1.83e-05 |    4.49e-05 | 0.0001207 |   0.0002155 | 0.0001223 |    0.0483745 |          8.0e-07 |          1.25e-05 |             2.5e-06 |          6.32e-05 |      3.41e-05 |      0.0002255 |         0.0001165 |      4.2e-06 |  9.20e-06 |    0.0001390 |   0.0001307 |  0.0008505 |    0.0001648 | 0.0000674 |        1.66e-05 |       9.2e-06 |            3.91e-05 |     5.0e-06 |  0.0006849 |     1.50e-05 |   0.0092501 | 0.0002788 |  0.0015912 | 0.0012483 |      0.0009063 |     5.24e-05 |    4.49e-05 |     0.0001856 | 0.0028062 | 8.82e-05 | 0.0001015 |     0.0002239 |        0.0001748 |    5.16e-05 |     9.99e-05 |     6.32e-05 |   0.0009812 | 0.0006733 |       6.08e-05 |      0.0026731 |        1.83e-05 | 4.74e-05 |        0.0001257 |    3.83e-05 |       0.0002264 |     0.0075358 |    0.0006533 |         0.0001157 |       0.0210211 |      2.25e-05 |   6.99e-05 |      5.08e-05 |            0.0001581 |    0.0000999 |         9.32e-05 |   9.82e-05 |       0.0009737 |     0.0007082 |     0.0001365 |       0.0002788 |      0.0011318 |        0.0003412 |  0.0002929 |        1.50e-05 |         4.2e-06 |    7.50e-06 |       8.66e-05 |    6.66e-05 |         7.5e-06 |    1.25e-05 |      2.41e-05 |      5.0e-06 | 5.49e-05 |   9.2e-06 |       0.0002971 |         0.0002413 |   0.0008139 |          0.0019957 |             0.0008655 |     0.0001631 |    0.0001615 |      7.24e-05 |    0.0002014 |     3.50e-05 |           0.0001240 |    2.75e-05 |        9.24e-05 |        9.15e-05 |     8.41e-05 |       2.5e-06 |            0.0006816 |   0.0012350 |   0.0002971 |       4.2e-06 |           0.0001315 |  0.0006699 | 0.0001956 |        0.0001357 |    0.0000965 |      0.0001215 |  5.33e-05 | 0.0005576 |   6.70e-06 | 0.0017984 |     0.0001440 |   0.0002089 |   1.58e-05 |   2.16e-05 | 4.99e-05 |    1.58e-05 |     3.99e-05 | 5.0e-06 |    0.0003379 |      6.7e-06 |   6.66e-05 | 1.91e-05 |        9.20e-06 |  0.0000558 |   0.0039406 |         0.0000699 |         3.74e-05 |      5.8e-06 | 0.0014372 |   0.0003387 | 4.91e-05 |      7.5e-06 | 5.58e-05 |  6.08e-05 |       1.41e-05 |  6.16e-05 | 0.0002164 |     7.16e-05 |   0.0003021 |   0.0001049 |   0.0001723 | 0.0001007 | 2.50e-05 | 0.0006425 |    0.0004169 |           3.3e-06 |  5.24e-05 |            5.0e-06 |                                           5.33e-05 |                                       0.0003137 |                                        2.66e-05 |                                    6.41e-05 |                            0.0007731 |                                 0.0003437 |                                        8.3e-06 |                                   0.0003146 |                                      1.58e-05 |                                            1.41e-05 |                                 4.91e-05 |                                  1.41e-05 |                                                           0.0001731 |                                 0.0036817 |                                        7.74e-05 |                                  9.32e-05 |                                       0.0045048 |                                 0.0001615 |                                 0.0002364 |                                     0.0001689 |                                    0.0001939 |                                       0.0002971 |                                   0.0046013 |                                    1.08e-05 |                              0.0000691 |                               0.0000499 |                                2.50e-06 |                                     3.41e-05 |                                    6.66e-05 |                                   4.2e-06 |                                   0.0092876 |                                0.0003287 |                                  7.74e-05 |                                                 0.0000458 |                                       3.91e-05 |                                 8.82e-05 |                                2.58e-05 |                            7.32e-05 |                                    9.57e-05 |
| ShotgunWGS-ControlPig5GutMicrobiome-Day7  | 5   | Control | Day 7       | Control Day 7         |   0.0012072 |      7.49e-05 |   0.0006482 |    2.08e-05 |     0.0002852 |    0.0002562 |     0.0000902 |       0.0096632 |    5.2e-06 |       2.93e-05 |     6.67e-05 |          7.91e-05 |      0.0001110 |    0.0000990 |  0.0003871 |        0.0001449 |     0.0002917 |      0.0005789 |   0.0002028 |      5.83e-05 |  0.0000729 |      1.53e-05 | 0.0005483 |   8.8e-06 | 1.01e-05 |       0.0001247 |     0.0001608 | 9.80e-06 |     6.2e-06 |   0.0003953 |   0.0001673 |   0.0001110 |    0.0000983 |       4.95e-05 |        0.0003809 |  0.0002106 | 0.0011405 |       0.0000879 |    0.0039948 |       6.67e-05 |         2.30e-06 |   0.0001087 |      0.0002058 |   3.42e-05 | 0.0002653 |      5.47e-05 | 0.0001459 |      4.53e-05 |    0.0010314 |    0.0004802 |        0.0003692 |    0.0012463 |     0.0028801 |  1.40e-05 |     0.0003718 | 0.0001426 |       0.0001703 |     0.0001426 |  0.0001498 |   0.0001354 |    0.0002393 |     1.6e-06 |    2.87e-05 |    2.83e-05 |      4.07e-05 | 0.0041645 |     3.22e-05 | 0.0001322 |     5.99e-05 |     7.16e-05 |    7.33e-05 | 2.0e-06 | 0.0072689 |   0.0733720 |   6.71e-05 | 0.0002191 |    0.0000957 |  2.18e-05 |     3.29e-05 |   2.08e-05 |    0.0000710 |       0.0215087 |       0.0000977 |        3.13e-05 | 0.0054932 |  0.0002702 | 0.0000853 |    7.80e-06 |           2.00e-06 |       0.0000661 |   0.0013678 |      0.0002178 |     0.0004516 |      0.0000697 |      3.42e-05 | 0.0000840 | 1.17e-05 | 6.58e-05 | 0.0003572 |    0.0008638 |    0.0094401 |       3.65e-05 |        0.0018766 |            0.0026456 |        0.0000794 |   9.80e-06 |     1.33e-05 |     0.0019229 | 1.17e-05 |                 0.0001003 |                7.52e-05 |                 0.0002631 |               1.86e-05 |               2.70e-05 |               0.0004040 |                1.86e-05 |               2.15e-05 |             0.0001853 |                 8.5e-06 |                2.34e-05 |               5.70e-05 |                  5.76e-05 |                    1.73e-05 |            8.80e-06 |           1.6e-06 |             0.0003002 |          1.01e-05 |      0.0006124 |        0.0011304 |        2.51e-05 |       6.41e-05 |       0.0012369 |      6.28e-05 |  5.80e-05 |   0.0001778 |    0.0000622 |        0.0009132 |  0.0001084 |    7.80e-06 |    5.2e-06 |     0.0000794 |    0.0005675 |  2.80e-05 |      3.48e-05 |      3.68e-05 |   1.3e-06 |     0.0002136 |  0.0005825 |    0.0002299 |      0.0001426 |     2.9e-06 |       0.0001296 |        0.0001185 |        0.0000967 |       7.49e-05 |    1.01e-05 |   0.0001530 |       9.40e-06 |   0.0001009 |   2.60e-06 |   0.0886641 |     1.01e-05 |   0.0033316 | 0.0001162 | 0.0000632 |     9.67e-05 |       4.88e-05 |    6.20e-06 |     0.0003871 |   0.0061049 |         0.0000973 |        0.0000983 |       0.0005945 | 4.95e-05 |    0.0001410 |     4.88e-05 |   0.0001107 |       0.0008725 |        1.37e-05 |   0.0003116 |         2.0e-06 |   5.2e-06 | 0.0001953 |    1.3e-06 |  0.0004184 |           1.33e-05 | 0.0004356 |      4.6e-06 |     0.0002680 |     0.0001211 |       0.0008107 |       0.0000736 |   0.0002686 | 0.0001071 |     0.0001791 |    1.14e-05 |     0.0000908 |         0.0001980 |          0.0043077 |        0.0001374 |     0.0001680 |       0.0000850 |        0.0001449 |            1.82e-05 |    0.0002191 |        0.0027140 |     0.0024063 |           8.89e-05 |        6.77e-05 |        7.80e-06 |      0.0001556 |    0.0001696 |         0.0003060 | 0.0042192 |      7.94e-05 | 0.0001563 |    0.0002920 |      5.01e-05 |        3.97e-05 | 0.0001195 | 0.0107353 |   0.0004021 |    0.0001354 |   0.0025304 |  2.02e-05 |  1.56e-05 |      0.0001465 |   1.47e-05 |   1.0e-06 |         1.3e-06 |   1.20e-05 |      1.63e-05 | 0.0000716 |    0.0001644 |    0.0029729 |       2.47e-05 |                2.0e-06 |    0.0001755 |    4.23e-05 |     1.27e-05 | 0.0000967 |      0.0000726 |     0.0000964 |   0.0005408 |      0.0053714 |   0.0589401 |       0.0005489 |        0.0373568 |  0.0000879 |    2.54e-05 |    3.06e-05 |        0.0002038 |   0.0013153 |  0.0002188 |       1.60e-05 |  0.0007654 |      0.0010939 |   0.0001726 | 0.0002898 |    1.17e-05 |     0.0021713 |    5.93e-05 |         7.0e-07 |   0.0011454 | 5.63e-05 | 3.32e-05 |    0.0000687 |   0.0017737 | 0.0014238 |         5.34e-05 | 1.17e-05 |   2.47e-05 |   0.0001211 |          4.33e-05 |      4.98e-05 | 3.74e-05 |    0.0e+00 | 0.0004001 |      4.69e-05 |      0.0001400 |    4.9e-06 |   0.0004304 | 0.0001504 |        1.14e-05 |     0.0002432 |  0.0000902 |   2.80e-05 |      2.67e-05 |  1.37e-05 |        1.79e-05 |      1.27e-05 |  4.30e-05 |      1.66e-05 |    2.60e-05 |       6.87e-05 |   1.30e-05 |      2.25e-05 |       0.0005639 |         3.61e-05 |    0.0006163 |      0.0013590 |       7.81e-05 |     0.0000980 |     0.0001325 | 3.45e-05 |   0.0001957 | 1.14e-05 |  0.0035511 |       2.21e-05 |        3.26e-05 |         4.75e-05 |     1.33e-05 |       2.08e-05 |   6.12e-05 |  0.0001546 |   1.37e-05 |     7.8e-06 |  0.0006264 |        3.13e-05 |   4.40e-05 |   4.40e-05 |         0.0001130 | 0.0000733 |   0.0000791 |  6.22e-05 |           2.41e-05 |   0.0001022 | 1.99e-05 |  0.0002061 |      7.50e-06 | 0.0000423 | 5.01e-05 | 0.0001149 |  5.47e-05 |      4.62e-05 |   4.10e-05 |         7.0e-07 |  2.34e-05 |  6.2e-06 |   5.5e-06 |     0.1788730 |   0.0006108 |           0.0000264 |   6.12e-05 | 0.0001481 |      0.0003174 |        0.0003262 |  0.0001784 | 0.0000563 |   2.05e-05 |     5.67e-05 |     7e-07 |  0.0001579 |        0.0e+00 |  0.0000967 |    0.0006316 |   0.0003917 |    3.09e-05 | 0.0012903 | 4.20e-06 |      7.8e-06 |   1.66e-05 | 3.16e-05 | 2.18e-05 |      0.0003080 |   0.0001374 |    1.53e-05 |     0.0001374 |        0.0001944 |   5.90e-06 |  0.0000866 |  0.0002360 |   4.13e-05 |    0.0002015 |   0.0001843 |      2.44e-05 |       3.16e-05 | 0.0001589 |   0.0052210 |   0.0001566 |  0.0000560 |     0.0001074 |       1.63e-05 |          0.0006486 |          0.0001735 |     6.87e-05 |        0.0001107 |     0.0003865 |          0.0003588 |      0.0001175 |        3.09e-05 |          4.46e-05 |     0.0000970 |     3.84e-05 |      8.37e-05 |    0.0000957 |      0.0005125 |      0.0001104 |         4.66e-05 |        0.0001208 |           0.0001293 |            8.10e-06 |       2.64e-05 |          4.23e-05 |   0.0001205 |        9.83e-05 |      2.38e-05 |        0.0002002 |     4.59e-05 |     0.0001299 |     1.95e-05 |     3.42e-05 |      6.38e-05 |     4.79e-05 |    4.2e-06 |    3.91e-05 |    2.21e-05 |    8.33e-05 |   1.86e-05 |       7.03e-05 |   0.0001032 |   0.0165768 |  0.0003607 |        2.3e-06 | 1.82e-05 | 0.0011503 |  3.48e-05 | 0.0000544 |        0.0001628 |     0.0004845 |  0.0003305 |  0.0001605 |         2.0e-06 |  1.47e-05 |    5.18e-05 |     1.99e-05 |      1.3e-06 |      0.0004503 |  1.79e-05 |     2.44e-05 | 4.92e-05 | 5.9e-06 | 0.0002989 |      1.07e-05 |    3.09e-05 |       2.93e-05 |   1.60e-05 |     0.0000964 |   0.0001029 |    3.78e-05 |     0.0001638 |    0.0001514 |       1.01e-05 |     9.60e-05 |   6.28e-05 | 0.0000674 |    0.0001345 |     5.27e-05 |  2.02e-05 | 0.0001996 |       0.0000944 |      6.5e-06 |     2.47e-05 |   3.61e-05 |      4.66e-05 |      0.0004350 |    0.0000713 |       1.95e-05 |   1.3e-06 |  0.0002256 |    2.96e-05 | 0.0029943 | 0.0002793 |    0.0015901 | 1.11e-05 |     2.38e-05 |       2.87e-05 |     3.65e-05 |   0.0002859 |       0.0000007 |         1.0e-06 |          7.0e-07 |     0.0023002 |    0.0012714 | 0.0001608 |       0.0071679 |      9.80e-06 |          2.0e-06 |   7.36e-05 |   1.27e-05 |     0.0002260 |     7.42e-05 |     3.26e-05 |   0.0001774 |    9.1e-06 |      0.0002246 |   0.0003835 |  0.0008931 |   1.11e-05 |  0.0006974 |   0.0002618 |     0.0011070 |    9.10e-06 |     0.0004975 |          0.0004252 |  1.17e-05 |      5.60e-05 | 0.0001966 |    5.20e-06 |      1.53e-05 |      9.80e-06 |         3.35e-05 |            1.0e-06 |      0.0003484 |    0.0001599 |     2.60e-05 | 3.9e-06 |    3.32e-05 |  7.98e-05 |    0.0001416 |  0.0000283 |     3.09e-05 |   3.9e-06 |    0.0002618 |   0.0002728 |        0.0001182 |  6.2e-06 |     0.0025476 | 2.0e-06 |  0.1651427 |       0.0002930 |         0.0002146 |         4.82e-05 | 0.0001400 |   0.0000824 |         0.0003090 |   0.0009214 |        0.0005196 |     1.73e-05 |     0.0001494 |      4.92e-05 |    0.0002377 |      0.0005493 |    2.90e-06 |    3.22e-05 |  0.0001390 | 0.0002292 |      6.8e-06 |      7.0e-07 | 6.87e-05 |      3.71e-05 | 0.0002116 |   0.0002188 |   0.0002009 |       2.67e-05 |    1.0e-06 |      0.0001927 |        0.0003321 |      0.0001905 |    0.0001885 |   6.97e-05 |      9.40e-06 |  0.0001322 |     0.0002061 | 0.0153119 |  1.11e-05 |   0.0004421 |   0.0001289 |   2.83e-05 |    5.18e-05 | 0.0000697 |   0.0001914 | 0.0001286 |    0.0422350 |          1.3e-06 |          1.27e-05 |             7.0e-07 |          3.65e-05 |      2.31e-05 |      0.0001872 |         0.0001143 |      2.6e-06 |  1.11e-05 |    0.0001234 |   0.0001127 |  0.0002725 |    0.0001127 | 0.0001856 |        1.50e-05 |       8.5e-06 |            1.92e-05 |     4.9e-06 |  0.0006779 |     1.66e-05 |   0.0068267 | 0.0001872 |  0.0010484 | 0.0000680 |      0.0009028 |     5.80e-05 |    2.21e-05 |     0.0001973 | 0.0030057 | 5.67e-05 | 0.0001107 |     0.0001895 |        0.0001761 |    5.18e-05 |     9.90e-05 |     7.55e-05 |   0.0007739 | 0.0005623 |       5.31e-05 |      0.0013756 |        1.76e-05 | 3.58e-05 |        0.0000866 |    4.53e-05 |       0.0002100 |     0.0133708 |    0.0005942 |         0.0000876 |       0.0141685 |      1.86e-05 |   4.07e-05 |      3.81e-05 |            0.0001328 |    0.0000964 |         8.47e-05 |   9.31e-05 |       0.0009048 |     0.0013508 |     0.0001280 |       0.0002507 |      0.0009702 |        0.0003187 |  0.0002953 |        1.37e-05 |         2.3e-06 |    6.80e-06 |       5.57e-05 |    5.99e-05 |         2.0e-06 |    7.80e-06 |      1.92e-05 |      6.2e-06 | 8.30e-05 |   4.6e-06 |       0.0002735 |         0.0001993 |   0.0007003 |          0.0017239 |             0.0008599 |     0.0001302 |    0.0001061 |      6.80e-05 |    0.0001748 |     2.80e-05 |           0.0000944 |    2.67e-05 |        8.14e-05 |        6.02e-05 |     7.07e-05 |       3.6e-06 |            0.0005691 |   0.0008680 |   0.0002946 |       1.6e-06 |           0.0001087 |  0.0006219 | 0.0001566 |        0.0001133 |    0.0001055 |      0.0000912 |  3.61e-05 | 0.0003002 |   7.20e-06 | 0.0011174 |     0.0001091 |   0.0001520 |   8.80e-06 |   1.92e-05 | 5.37e-05 |    1.50e-05 |     3.39e-05 | 4.6e-06 |    0.0002653 |      4.2e-06 |   5.21e-05 | 1.20e-05 |        5.90e-06 |  0.0001048 |   0.0042954 |         0.0001094 |         5.93e-05 |      1.3e-06 | 0.0009344 |   0.0003793 | 3.22e-05 |      7.2e-06 | 3.22e-05 |  3.87e-05 |       7.50e-06 |  4.72e-05 | 0.0000964 |     8.24e-05 |   0.0002540 |   0.0000628 |   0.0001169 | 0.0000964 | 1.69e-05 | 0.0003506 |    0.0003139 |           5.5e-06 |  6.25e-05 |            3.9e-06 |                                           4.33e-05 |                                       0.0003220 |                                        4.04e-05 |                                    3.68e-05 |                            0.0002735 |                                 0.0003340 |                                        9.8e-06 |                                   0.0001953 |                                      1.40e-05 |                                            7.20e-06 |                                 3.32e-05 |                                  8.50e-06 |                                                           0.0001201 |                                 0.0034580 |                                        6.19e-05 |                                  7.00e-05 |                                       0.0043846 |                                 0.0001260 |                                 0.0001605 |                                     0.0001384 |                                    0.0001136 |                                       0.0001976 |                                   0.0048931 |                                    1.04e-05 |                              0.0000563 |                               0.0000619 |                                2.00e-06 |                                     2.41e-05 |                                    3.09e-05 |                                   3.3e-06 |                                   0.0079575 |                                0.0002338 |                                  4.20e-05 |                                                 0.0000641 |                                       4.53e-05 |                                 3.97e-05 |                                1.33e-05 |                            6.51e-05 |                                    4.20e-06 |
| ShotgunWGS-TomatoPig18GutMicrobiome-Day7  | 18  | Tomato  | Day 7       | Tomato Day 7          |   0.0006501 |      8.19e-05 |   0.0003281 |    1.85e-05 |     0.0001486 |    0.0001896 |     0.0001094 |       0.0059482 |    3.4e-06 |       4.77e-05 |     7.80e-05 |          7.57e-05 |      0.0000959 |    0.0001464 |  0.0006512 |        0.0000813 |     0.0002064 |      0.0003556 |   0.0002956 |      9.65e-05 |  0.0000561 |      2.75e-05 | 0.0001997 |   6.7e-06 | 7.30e-06 |       0.0000645 |     0.0001498 | 9.50e-06 |     5.0e-06 |   0.0004022 |   0.0002872 |   0.0000583 |    0.0001049 |       7.74e-05 |        0.0002737 |  0.0001032 | 0.0009177 |       0.0000780 |    0.0024922 |       5.05e-05 |         5.16e-05 |   0.0000415 |      0.0001099 |   1.91e-05 | 0.0001262 |      8.30e-05 | 0.0001285 |      2.97e-05 |    0.0007864 |    0.0002524 |        0.0003326 |    0.0008061 |     0.0016144 |  1.51e-05 |     0.0002620 | 0.0001027 |       0.0001369 |     0.0000830 |  0.0001088 |   0.0001027 |    0.0003999 |     1.1e-06 |    2.08e-05 |    4.43e-05 |      5.83e-05 | 0.0085061 |     2.64e-05 | 0.0001105 |     6.17e-05 |     5.89e-05 |    4.71e-05 | 2.8e-06 | 0.0049581 |   0.0944738 |   3.76e-05 | 0.0001200 |    0.0000819 |  1.46e-05 |     4.09e-05 |   1.12e-05 |    0.0001139 |       0.0585603 |       0.0000987 |        2.52e-05 | 0.0040758 |  0.0002507 | 0.0000752 |    9.50e-06 |           7.30e-06 |       0.0000931 |   0.0011628 |      0.0002592 |     0.0002816 |      0.0001212 |      5.72e-05 | 0.0000752 | 2.75e-05 | 1.80e-05 | 0.0002120 |    0.0007567 |    0.0061681 |       3.25e-05 |        0.0012464 |            0.0018191 |        0.0000611 |   7.90e-06 |     6.20e-06 |     0.0006967 | 5.60e-06 |                 0.0000724 |                4.88e-05 |                 0.0002345 |               1.07e-05 |               1.40e-05 |               0.0002406 |                4.50e-06 |               1.40e-05 |             0.0001975 |                 3.9e-06 |                2.41e-05 |               2.69e-05 |                  3.87e-05 |                    1.51e-05 |            2.20e-06 |           0.0e+00 |             0.0004336 |          5.60e-06 |      0.0006372 |        0.0006221 |        1.12e-05 |       3.20e-05 |       0.0006283 |      9.98e-05 |  3.25e-05 |   0.0002098 |    0.0001127 |        0.0004213 |  0.0001251 |    6.20e-06 |    2.2e-06 |     0.0000954 |    0.0006995 |  2.19e-05 |      3.59e-05 |      2.86e-05 |   0.0e+00 |     0.0001868 |  0.0005901 |    0.0002087 |      0.0001402 |     9.0e-06 |       0.0000970 |        0.0000740 |        0.0001139 |       9.42e-05 |    7.30e-06 |   0.0002575 |       1.07e-05 |   0.0002227 |   8.40e-06 |   0.0567714 |     6.70e-06 |   0.0087798 | 0.0000836 | 0.0001027 |     8.41e-05 |       3.76e-05 |    6.70e-06 |     0.0002182 |   0.0031340 |         0.0000567 |        0.0000965 |       0.0010394 | 4.09e-05 |    0.0001144 |     3.93e-05 |   0.0001060 |       0.0009912 |        5.60e-06 |   0.0003332 |         6.0e-07 |   3.4e-06 | 0.0003949 |    2.8e-06 |  0.0003444 |           1.12e-05 | 0.0004426 |      6.2e-06 |     0.0001515 |     0.0000583 |       0.0005839 |       0.0001290 |   0.0002485 | 0.0001391 |     0.0001301 |    2.08e-05 |     0.0000595 |         0.0001172 |          0.0027587 |        0.0001094 |     0.0001335 |       0.0000611 |        0.0000976 |            1.85e-05 |    0.0001677 |        0.0016553 |     0.0010131 |           5.10e-05 |        6.45e-05 |        7.90e-06 |      0.0001217 |    0.0000578 |         0.0001761 | 0.0018118 |      5.67e-05 | 0.0001111 |    0.0002373 |      3.48e-05 |        4.49e-05 | 0.0001234 | 0.0048387 |   0.0004465 |    0.0000724 |   0.0024844 |  1.68e-05 |  1.12e-05 |      0.0000869 |   1.57e-05 |   0.0e+00 |         6.0e-07 |   3.40e-06 |      7.90e-06 | 0.0000174 |    0.0001975 |    0.0023369 |       1.68e-05 |                0.0e+00 |    0.0000954 |    2.30e-05 |     1.57e-05 | 0.0000898 |      0.0000359 |     0.0001200 |   0.0090541 |      0.0048308 |   0.0372300 |       0.0003798 |        0.0766393 |  0.0000337 |    1.96e-05 |    1.74e-05 |        0.0001374 |   0.0011617 |  0.0000948 |       2.24e-05 |  0.0004942 |      0.0016200 |   0.0001251 | 0.0004297 |    1.01e-05 |     0.0013552 |    3.70e-05 |         2.8e-06 |   0.0027166 | 2.47e-05 | 5.10e-05 |    0.0001167 |   0.0011976 | 0.0011370 |         8.98e-05 | 2.20e-06 |   3.93e-05 |   0.0001296 |          4.09e-05 |      5.95e-05 | 5.16e-05 |    6.0e-07 | 0.0003803 |      7.52e-05 |      0.0000836 |    3.4e-06 |   0.0002418 | 0.0000948 |        7.90e-06 |     0.0001257 |  0.0000808 |   3.03e-05 |      1.12e-05 |  7.90e-06 |        1.12e-05 |      1.18e-05 |  3.37e-05 |      1.46e-05 |    1.23e-05 |       8.02e-05 |   7.90e-06 |      2.19e-05 |       0.0003607 |         3.81e-05 |    0.0005200 |      0.0007629 |       9.70e-05 |     0.0000797 |     0.0001369 | 2.24e-05 |   0.0001027 | 8.40e-06 |  0.0027587 |       1.07e-05 |        1.46e-05 |         3.98e-05 |     7.90e-06 |       4.71e-05 |   5.55e-05 |  0.0000684 |   5.00e-06 |     3.4e-06 |  0.0003579 |        7.29e-05 |   9.37e-05 |   4.94e-05 |         0.0001105 | 0.0001335 |   0.0000578 |  3.03e-05 |           2.47e-05 |   0.0001829 | 5.00e-06 |  0.0003046 |      1.29e-05 | 0.0000639 | 4.82e-05 | 0.0000797 |  9.82e-05 |      3.48e-05 |   6.00e-05 |         1.1e-06 |  2.52e-05 |  2.8e-06 |   9.0e-06 |     0.1042975 |   0.0005228 |           0.0001430 |   4.15e-05 | 0.0001015 |      0.0003433 |        0.0004684 |  0.0001329 | 0.0001071 |   1.74e-05 |     2.86e-05 |     0e+00 |  0.0001318 |        1.1e-06 |  0.0001879 |    0.0004465 |   0.0003758 |    4.04e-05 | 0.0009199 | 1.51e-05 |      6.2e-06 |   2.36e-05 | 1.68e-05 | 1.63e-05 |      0.0002059 |   0.0000735 |    7.30e-06 |     0.0001127 |        0.0001509 |   5.60e-06 |  0.0000684 |  0.0002485 |   4.60e-05 |    0.0001257 |   0.0001167 |      1.29e-05 |       2.58e-05 | 0.0001357 |   0.0064643 |   0.0001402 |  0.0000370 |     0.0001481 |       1.12e-05 |          0.0003876 |          0.0000987 |     3.93e-05 |        0.0000892 |     0.0002395 |          0.0002721 |      0.0000959 |        2.19e-05 |          2.08e-05 |     0.0000645 |     3.14e-05 |      6.84e-05 |    0.0000533 |      0.0002850 |      0.0000959 |         3.87e-05 |        0.0000791 |           0.0000825 |            1.07e-05 |       1.01e-05 |          4.04e-05 |   0.0002311 |        7.24e-05 |      2.80e-05 |        0.0002474 |     7.24e-05 |     0.0001245 |     8.40e-06 |     7.01e-05 |      4.04e-05 |     3.98e-05 |    5.6e-06 |    4.99e-05 |    2.30e-05 |    8.13e-05 |   1.46e-05 |       8.47e-05 |   0.0000864 |   0.0049536 |  0.0005694 |        3.4e-06 | 1.35e-05 | 0.0006754 |  1.63e-05 | 0.0000191 |        0.0001784 |     0.0007500 |  0.0002221 |  0.0001834 |         0.0e+00 |  1.57e-05 |    8.86e-05 |     2.13e-05 |      6.0e-07 |      0.0002850 |  1.07e-05 |     7.90e-06 | 2.36e-05 | 5.0e-06 | 0.0001778 |      7.90e-06 |    3.20e-05 |       1.07e-05 |   9.50e-06 |     0.0000628 |   0.0001172 |    3.31e-05 |     0.0001285 |    0.0001464 |       8.40e-06 |     6.28e-05 |   4.54e-05 | 0.0000993 |    0.0002754 |     6.84e-05 |  1.29e-05 | 0.0001946 |       0.0001335 |      7.3e-06 |     3.03e-05 |   4.66e-05 |      4.43e-05 |      0.0002765 |    0.0000774 |       1.23e-05 |   0.0e+00 |  0.0001930 |    3.48e-05 | 0.0061911 | 0.0004028 |    0.0009514 | 7.30e-06 |     2.19e-05 |       2.02e-05 |     2.92e-05 |   0.0000864 |       0.0000673 |         1.1e-06 |          8.4e-06 |     0.0015790 |    0.0017703 | 0.0000987 |       0.0087316 |      8.40e-06 |          6.0e-07 |   6.23e-05 |   1.29e-05 |     0.0004970 |     8.02e-05 |     3.03e-05 |   0.0000639 |    9.5e-06 |      0.0001380 |   0.0004347 |  0.0009979 |   1.68e-05 |  0.0005082 |   0.0002822 |     0.0006815 |    7.30e-06 |     0.0002535 |          0.0002137 |  1.40e-05 |      4.21e-05 | 0.0001397 |    6.70e-06 |      2.02e-05 |      9.50e-06 |         3.98e-05 |            1.1e-06 |      0.0001616 |    0.0000869 |     1.40e-05 | 3.9e-06 |    2.19e-05 |  9.42e-05 |    0.0001700 |  0.0000409 |     3.93e-05 |   6.2e-06 |    0.0002889 |   0.0005312 |        0.0001862 |  3.4e-06 |     0.0026594 | 2.2e-06 |  0.2682445 |       0.0004056 |         0.0003629 |         3.37e-05 | 0.0000746 |   0.0000337 |         0.0001874 |   0.0008173 |        0.0003281 |     8.40e-06 |     0.0000993 |      5.22e-05 |    0.0000993 |      0.0003136 |    3.90e-06 |    2.36e-05 |  0.0001004 | 0.0002132 |      6.2e-06 |      1.1e-06 | 4.09e-05 |      6.23e-05 | 0.0002042 |   0.0002563 |   0.0002760 |       3.20e-05 |    0.0e+00 |      0.0002008 |        0.0003012 |      0.0001526 |    0.0001767 |   7.18e-05 |      2.20e-06 |  0.0001027 |     0.0002047 | 0.0088533 |  9.50e-06 |   0.0004117 |   0.0001498 |   3.76e-05 |    5.16e-05 | 0.0000864 |   0.0002025 | 0.0001212 |    0.0264808 |          6.0e-07 |          5.00e-06 |             3.9e-06 |          4.60e-05 |      6.51e-05 |      0.0002350 |         0.0002031 |      1.1e-06 |  1.07e-05 |    0.0001105 |   0.0001851 |  0.0006002 |    0.0001700 | 0.0003214 |        1.35e-05 |       2.2e-06 |            2.24e-05 |     4.5e-06 |  0.0004824 |     1.29e-05 |   0.0013463 | 0.0001386 |  0.0006204 | 0.0011297 |      0.0005278 |     4.71e-05 |    3.59e-05 |     0.0001986 | 0.0028008 | 2.58e-05 | 0.0001077 |     0.0002216 |        0.0001481 |    7.12e-05 |     1.29e-04 |     9.48e-05 |   0.0004717 | 0.0008341 |       5.67e-05 |      0.0016879 |        1.29e-05 | 4.43e-05 |        0.0000813 |    6.34e-05 |       0.0002008 |     0.0057126 |    0.0008779 |         0.0001212 |       0.0215249 |      2.52e-05 |   2.97e-05 |      2.41e-05 |            0.0001419 |    0.0000561 |         3.93e-05 |   4.21e-05 |       0.0006226 |     0.0019566 |     0.0001156 |       0.0001705 |      0.0007130 |        0.0001666 |  0.0001784 |        2.08e-05 |         2.8e-06 |    2.20e-06 |       4.94e-05 |    4.54e-05 |         2.2e-06 |    5.60e-06 |      1.57e-05 |      3.4e-06 | 6.45e-05 |   1.7e-06 |       0.0001924 |         0.0001251 |   0.0003489 |          0.0011656 |             0.0005121 |     0.0001111 |    0.0001694 |      8.98e-05 |    0.0000959 |     1.74e-05 |           0.0000909 |    2.30e-05 |        7.91e-05 |        8.47e-05 |     4.94e-05 |       2.2e-06 |            0.0003237 |   0.0003781 |   0.0001750 |       5.6e-06 |           0.0001060 |  0.0004005 | 0.0001228 |        0.0000909 |    0.0000847 |      0.0000471 |  4.99e-05 | 0.0000729 |   5.60e-06 | 0.0008106 |     0.0000875 |   0.0001043 |   1.63e-05 |   2.13e-05 | 6.23e-05 |    1.35e-05 |     4.66e-05 | 5.6e-06 |    0.0001178 |      1.7e-06 |   3.76e-05 | 2.13e-05 |        1.46e-05 |  0.0001733 |   0.0028142 |         0.0001761 |         8.47e-05 |      1.7e-06 | 0.0005396 |   0.0001329 | 4.04e-05 |      3.4e-06 | 2.69e-05 |  3.81e-05 |       6.20e-06 |  3.59e-05 | 0.0000701 |     8.64e-05 |   0.0002721 |   0.0000376 |   0.0001840 | 0.0000595 | 1.01e-05 | 0.0003068 |    0.0004129 |           6.7e-06 |  6.06e-05 |            2.8e-06 |                                           7.80e-05 |                                       0.0001761 |                                        2.97e-05 |                                    5.60e-06 |                            0.0006776 |                                 0.0002193 |                                        2.8e-06 |                                   0.0000668 |                                      8.40e-06 |                                            9.50e-06 |                                 2.30e-05 |                                  6.70e-06 |                                                           0.0000707 |                                 0.0024249 |                                        4.71e-05 |                                  4.21e-05 |                                       0.0022275 |                                 0.0000791 |                                 0.0001526 |                                     0.0001638 |                                    0.0000920 |                                       0.0001531 |                                   0.0029915 |                                    3.40e-06 |                              0.0000589 |                               0.0000561 |                                1.07e-05 |                                     1.80e-05 |                                    1.40e-05 |                                   4.5e-06 |                                   0.0069882 |                                0.0002317 |                                  2.69e-05 |                                                 0.0001021 |                                       4.54e-05 |                                 6.70e-06 |                                5.00e-06 |                            4.88e-05 |                                    3.31e-05 |

``` r
# move Sample_Name to rownames, remove metadata
RelAbund.Genus.Filt.zerofilt.alphadiv <- RelAbund.Genus.Filt.zerofilt

rownames(RelAbund.Genus.Filt.zerofilt.alphadiv) <- RelAbund.Genus.Filt.zerofilt.alphadiv$Sample_Name  

RelAbund.Genus.Filt.zerofilt.alphadiv[1:5,1:8]
```

    ##                                                                         Sample_Name
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14 ShotgunWGS-ControlPig6GutMicrobiome-Day14
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0   ShotgunWGS-ControlPig8GutMicrobiome-Day0
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14 ShotgunWGS-ControlPig3GutMicrobiome-Day14
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7   ShotgunWGS-TomatoPig14GutMicrobiome-Day7
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7   ShotgunWGS-ControlPig5GutMicrobiome-Day7
    ##                                           Pig    Diet Time_Point
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0    8 Control      Day 0
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7   14  Tomato      Day 7
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7    5 Control      Day 7
    ##                                           Diet_By_Time_Point Abiotrophia
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14     Control Day 14 0.001305713
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0       Control Day 0 0.001347804
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14     Control Day 14 0.001066255
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7        Tomato Day 7 0.001311580
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7       Control Day 7 0.001207244
    ##                                           Acaryochloris  Acetivibrio
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14  6.983388e-05 0.0005122869
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0   9.904370e-05 0.0007097339
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14  6.914992e-05 0.0005363651
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7   1.090209e-04 0.0008422076
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7   7.488298e-05 0.0006482261

``` r
# remove metadata
RelAbund.Genus.Filt.zerofilt.alphadiv <- RelAbund.Genus.Filt.zerofilt.alphadiv %>%
  select(Abiotrophia:ncol(.))

RelAbund.Genus.Filt.zerofilt.alphadiv[1:5,1:5]
```

    ##                                           Abiotrophia Acaryochloris
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14 0.001305713  6.983388e-05
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0  0.001347804  9.904370e-05
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14 0.001066255  6.914992e-05
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7  0.001311580  1.090209e-04
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7  0.001207244  7.488298e-05
    ##                                            Acetivibrio  Acetobacter
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14 0.0005122869 1.700751e-05
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0  0.0007097339 2.047538e-05
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14 0.0005363651 1.553931e-05
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7  0.0008422076 3.412106e-05
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7  0.0006482261 2.083700e-05
    ##                                           Acetohalobium
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14  0.0002669664
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0   0.0003268918
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14  0.0002628733
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7   0.0003320562
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7   0.0002852065

``` r
rownames(RelAbund.Genus.Filt.zerofilt.alphadiv)
```

    ##  [1] "ShotgunWGS-ControlPig6GutMicrobiome-Day14" 
    ##  [2] "ShotgunWGS-ControlPig8GutMicrobiome-Day0"  
    ##  [3] "ShotgunWGS-ControlPig3GutMicrobiome-Day14" 
    ##  [4] "ShotgunWGS-TomatoPig14GutMicrobiome-Day7"  
    ##  [5] "ShotgunWGS-ControlPig5GutMicrobiome-Day7"  
    ##  [6] "ShotgunWGS-TomatoPig18GutMicrobiome-Day7"  
    ##  [7] "ShotgunWGS-TomatoPig16GutMicrobiome-Day7"  
    ##  [8] "ShotgunWGS-ControlPig10GutMicrobiome-Day7" 
    ##  [9] "ShotgunWGS-ControlPig2GutMicrobiome-Day0"  
    ## [10] "ShotgunWGS-TomatoPig18GutMicrobiome-Day0"  
    ## [11] "ShotgunWGS-ControlPig10GutMicrobiome-Day0" 
    ## [12] "ShotgunWGS-ControlPig7GutMicrobiome-Day0"  
    ## [13] "ShotgunWGS-ControlPig8GutMicrobiome-Day14" 
    ## [14] "ShotgunWGS-TomatoPig11GutMicrobiome-Day0"  
    ## [15] "ShotgunWGS-TomatoPig19GutMicrobiome-Day0"  
    ## [16] "ShotgunWGS-TomatoPig17GutMicrobiome-Day14" 
    ## [17] "ShotgunWGS-ControlPig9GutMicrobiome-Day14" 
    ## [18] "ShotgunWGS-ControlPig10GutMicrobiome-Day14"
    ## [19] "ShotgunWGS-TomatoPig19GutMicrobiome-Day7"  
    ## [20] "ShotgunWGS-ControlPig5GutMicrobiome-Day14" 
    ## [21] "ShotgunWGS-ControlPig2GutMicrobiome-Day7"  
    ## [22] "ShotgunWGS-ControlPig6GutMicrobiome-Day7"  
    ## [23] "ShotgunWGS-TomatoPig12GutMicrobiome-Day0"  
    ## [24] "ShotgunWGS-TomatoPig14GutMicrobiome-Day0"  
    ## [25] "ShotgunWGS-ControlPig7GutMicrobiome-Day14" 
    ## [26] "ShotgunWGS-TomatoPig11GutMicrobiome-Day14" 
    ## [27] "ShotgunWGS-TomatoPig20GutMicrobiome-Day0"  
    ## [28] "ShotgunWGS-ControlPig9GutMicrobiome-Day0"  
    ## [29] "ShotgunWGS-TomatoPig11GutMicrobiome-Day7"  
    ## [30] "ShotgunWGS-TomatoPig13GutMicrobiome-Day7"  
    ## [31] "ShotgunWGS-TomatoPig17GutMicrobiome-Day0"  
    ## [32] "ShotgunWGS-TomatoPig19GutMicrobiome-Day14" 
    ## [33] "ShotgunWGS-TomatoPig13GutMicrobiome-Day0"  
    ## [34] "ShotgunWGS-ControlPig2GutMicrobiome-Day14" 
    ## [35] "ShotgunWGS-ControlPig1GutMicrobiome-Day7"  
    ## [36] "ShotgunWGS-TomatoPig15GutMicrobiome-Day7"  
    ## [37] "ShotgunWGS-TomatoPig15GutMicrobiome-Day0"  
    ## [38] "ShotgunWGS-TomatoPig12GutMicrobiome-Day7"  
    ## [39] "ShotgunWGS-TomatoPig14GutMicrobiome-Day14" 
    ## [40] "ShotgunWGS-TomatoPig20GutMicrobiome-Day14" 
    ## [41] "ShotgunWGS-ControlPig1GutMicrobiome-Day0"  
    ## [42] "ShotgunWGS-ControlPig4GutMicrobiome-Day14" 
    ## [43] "ShotgunWGS-ControlPig6GutMicrobiome-Day0"  
    ## [44] "ShotgunWGS-TomatoPig16GutMicrobiome-Day0"  
    ## [45] "ShotgunWGS-TomatoPig16GutMicrobiome-Day14" 
    ## [46] "ShotgunWGS-TomatoPig18GutMicrobiome-Day14" 
    ## [47] "ShotgunWGS-ControlPig7GutMicrobiome-Day7"  
    ## [48] "ShotgunWGS-ControlPig4GutMicrobiome-Day7"  
    ## [49] "ShotgunWGS-TomatoPig13GutMicrobiome-Day14" 
    ## [50] "ShotgunWGS-ControlPig8GutMicrobiome-Day7"  
    ## [51] "ShotgunWGS-TomatoPig15GutMicrobiome-Day14" 
    ## [52] "ShotgunWGS-TomatoPig12GutMicrobiome-Day14" 
    ## [53] "ShotgunWGS-TomatoPig20GutMicrobiome-Day7"  
    ## [54] "ShotgunWGS-ControlPig1GutMicrobiome-Day14" 
    ## [55] "ShotgunWGS-ControlPig3GutMicrobiome-Day0"  
    ## [56] "ShotgunWGS-ControlPig5GutMicrobiome-Day0"  
    ## [57] "ShotgunWGS-ControlPig4GutMicrobiome-Day0"  
    ## [58] "ShotgunWGS-ControlPig9GutMicrobiome-Day7"  
    ## [59] "ShotgunWGS-ControlPig3GutMicrobiome-Day7"  
    ## [60] "ShotgunWGS-TomatoPig17GutMicrobime-Day7"

### Calculate alpha diversity

``` r
# run alpha diversity on phyla
genera.filt.div <- diversity(RelAbund.Genus.Filt.zerofilt.alphadiv, index = "shannon")

# convert to df
genera.filt.div.df <- as.data.frame(genera.filt.div)

# make column name 'shannon.phyla.filt'
colnames(genera.filt.div.df) <- "shannon.genera.filt"

head(genera.filt.div.df)
```

    ##                                           shannon.genera.filt
    ## ShotgunWGS-ControlPig6GutMicrobiome-Day14            3.454094
    ## ShotgunWGS-ControlPig8GutMicrobiome-Day0             3.724200
    ## ShotgunWGS-ControlPig3GutMicrobiome-Day14            3.405044
    ## ShotgunWGS-TomatoPig14GutMicrobiome-Day7             3.639925
    ## ShotgunWGS-ControlPig5GutMicrobiome-Day7             3.526991
    ## ShotgunWGS-TomatoPig18GutMicrobiome-Day7             3.281356

Combine shannon alpha diversity results with metadata

``` r
# compile genera metadata
genera.metadata <- RelAbund.Genus.Filt.zerofilt[,1:5]

# combine with metadata
genera.filt.div.df.meta <- cbind(genera.metadata, genera.filt.div.df)

head(genera.filt.div.df.meta)
```

    ##                                 Sample_Name Pig    Diet Time_Point
    ## 1 ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## 2  ShotgunWGS-ControlPig8GutMicrobiome-Day0   8 Control      Day 0
    ## 3 ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ## 4  ShotgunWGS-TomatoPig14GutMicrobiome-Day7  14  Tomato      Day 7
    ## 5  ShotgunWGS-ControlPig5GutMicrobiome-Day7   5 Control      Day 7
    ## 6  ShotgunWGS-TomatoPig18GutMicrobiome-Day7  18  Tomato      Day 7
    ##   Diet_By_Time_Point shannon.genera.filt
    ## 1     Control Day 14            3.454094
    ## 2      Control Day 0            3.724200
    ## 3     Control Day 14            3.405044
    ## 4       Tomato Day 7            3.639925
    ## 5      Control Day 7            3.526991
    ## 6       Tomato Day 7            3.281356

### Plotting

X axis by diet

``` r
alpha.diversity.genera.bydiet <- genera.filt.div.df.meta %>%
  ggplot(aes(x = Diet, y = shannon.genera.filt, fill = Diet_By_Time_Point)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Diet_By_Time_Point), color = "black", alpha = 0.7, position=position_jitterdodge()) +
  scale_fill_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  scale_color_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  labs(x=NULL, 
       y="Shannon diversity index", 
       title = "Alpha Diversity",
       subtitle = "Shannon Index, Genera Level", 
       fill="Diet & Time Point")

alpha.diversity.genera.bydiet
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-75-1.png)<!-- -->

``` r
ggsave("Figures/AlphaDiversityGenera_ByDiet_Boxplot.png", 
       plot = alpha.diversity.genera.bydiet, 
       dpi = 800, 
       width = 10, 
       height = 6)
```

X-axis by Day

``` r
genera.filt.div.df.meta <- genera.filt.div.df.meta %>%
  mutate(Time_Point = fct_relevel(Time_Point, c("Day 0", "Day 7", "Day 14")))

alpha.diversity.genera.bytime <- genera.filt.div.df.meta %>%
  ggplot(aes(x = Time_Point, y = shannon.genera.filt, fill = Diet_By_Time_Point)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Diet_By_Time_Point), color = "black", alpha = 0.7, position=position_jitterdodge()) +
  scale_fill_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  scale_color_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  labs(x=NULL, 
       y="Shannon diversity index", 
       title = "Alpha Diversity",
       subtitle = "Shannon Index, Genera Level", 
       fill="Diet & Time Point")

alpha.diversity.genera.bytime
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

``` r
ggsave("Figures/AlphaDiversityGenera_ByTime_Boxplot.png", 
       plot = alpha.diversity.genera.bytime, 
       dpi = 800, 
       width = 7, 
       height = 5)
```

### Statistics

Repeated measures ANOVA on Shannon alpha diversity

``` r
# must remove columns that aren't used in anova
head(genera.filt.div.df.meta) 
```

    ##                                 Sample_Name Pig    Diet Time_Point
    ## 1 ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## 2  ShotgunWGS-ControlPig8GutMicrobiome-Day0   8 Control      Day 0
    ## 3 ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ## 4  ShotgunWGS-TomatoPig14GutMicrobiome-Day7  14  Tomato      Day 7
    ## 5  ShotgunWGS-ControlPig5GutMicrobiome-Day7   5 Control      Day 7
    ## 6  ShotgunWGS-TomatoPig18GutMicrobiome-Day7  18  Tomato      Day 7
    ##   Diet_By_Time_Point shannon.genera.filt
    ## 1     Control Day 14            3.454094
    ## 2      Control Day 0            3.724200
    ## 3     Control Day 14            3.405044
    ## 4       Tomato Day 7            3.639925
    ## 5      Control Day 7            3.526991
    ## 6       Tomato Day 7            3.281356

``` r
genera.filt.div.df.meta.foranova <- genera.filt.div.df.meta[,-c(1,5)]

head(genera.filt.div.df.meta.foranova)
```

    ##   Pig    Diet Time_Point shannon.genera.filt
    ## 1   6 Control     Day 14            3.454094
    ## 2   8 Control      Day 0            3.724200
    ## 3   3 Control     Day 14            3.405044
    ## 4  14  Tomato      Day 7            3.639925
    ## 5   5 Control      Day 7            3.526991
    ## 6  18  Tomato      Day 7            3.281356

``` r
genera.filt.alphadiv.anova <- 
  anova_test(data = genera.filt.div.df.meta.foranova,
             formula = shannon.genera.filt ~ Diet*Time_Point + Error(Pig/Time_Point),
             dv = shannon.genera.filt, 
             wid = Pig, 
             within = Time_Point, 
             between = Diet)

get_anova_table(genera.filt.alphadiv.anova)
```

    ## ANOVA Table (type II tests)
    ## 
    ##            Effect DFn DFd     F     p p<.05   ges
    ## 1            Diet   1  18 2.888 0.106       0.066
    ## 2      Time_Point   2  36 0.254 0.777       0.008
    ## 3 Diet:Time_Point   2  36 0.037 0.964       0.001

-   Non-significant effect of diet (p = 0.106)  
-   Non-significant effect of timepoint (0.777)  
-   Non-significant interaction of diet:time point (p = 0.964).

Check for normality

``` r
shapiro.test(genera.filt.div.df.meta.foranova$shannon.genera.filt)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  genera.filt.div.df.meta.foranova$shannon.genera.filt
    ## W = 0.97004, p-value = 0.1466

Normal.

No need for posthoc test since no model parameters are significant.

## ALDEx2

Quick introduction to anatomy of the aldex function

The aldex function does every step - data transformation and
statistics  
variable.name &lt;- aldex(reads.data, variables.vector, mc.samples=\#,
test=“t”/“kw”, effect=T/F)  
reads.data - your reads/count data, unchanged  
variables.vector - a vector of the variables corresponding to sample
groups, in SAME order as sample names (and therefore columns)  
mc.samples - here you tell the function how many Monte Carlo sampels to
use with an integer (128 is typical)  
test - which test do you want, t-test and wilcoxon, or anova-like and
kruskal wallace? (will always do the parametric and non-parametric) t =
t-test and wilcoxon kw = anova-like and kruskal wallace  
effect - do you want to incude effect results in output?

Key to aldex outputs - taken directly from vignette

-   we.ep - Expected P value of Welch’s t test
-   we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t
    test  
-   wi.ep - Expected P value of Wilcoxon rank test
-   wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon
    test
-   kw.ep - Expected P value of Kruskal-Wallace test
-   kw.eBH - Expected Benjamini-Hochberg corrected P value of
    Kruskal-Wallace test
-   glm.ep - Expected P value of glm test
-   glm.eBH - Expected Benjamini-Hochberg corrected P value of glm test
-   rab.all - median clr value for all samples in the feature
-   rab.win.NS - median clr value for the NS group of samples
-   rab.win.S - median clr value for the S group of samples
-   dif.btw - median difference in clr values between S and NS groups
-   dif.win - median of the largest difference in clr values within S
    and NS groups
-   effect - median effect size: diff.btw / max(diff.win) for all
    instances
-   overlap - proportion of effect size that overlaps 0 (i.e. no effect)

ALDEx2 takes counts, not relative abundance.

We are using Benjamini Hochberg corrected pvalues, or `we.eBH` for
t-tests (i.e., subsetting by time), and Benjamini-Hochberg corrected
pvalues of the glm test `glm.eBH` for ANOVA tests (i.e., subsetting by
diet)

Downloading ALDEx2

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")
```

### Wrangling

Since we use counts for ALDEx2, we need to filter our counts data to
include only the genera we ended up using in our final analysis

``` r
# this is the data set filtered to remove inplausible phyla, but still includes genera with a lot of missing values 
Genus.Counts.Filt[1:10,1:10]
```

    ## # A tibble: 10 × 10
    ##    domain    phylum   class order family genus `ShotgunWGS-Co…` `ShotgunWGS-Co…`
    ##    <chr>     <chr>    <chr> <chr> <chr>  <chr>            <dbl>            <dbl>
    ##  1 Viruses   unclass… uncl… Caud… Podov… AHJD…               29                0
    ##  2 Bacteria  Firmicu… Baci… Lact… Aeroc… Abio…             5067             5661
    ##  3 Eukaryota unclass… uncl… uncl… uncla… Acan…                0                0
    ##  4 Bacteria  Cyanoba… uncl… uncl… uncla… Acar…              271              416
    ##  5 Bacteria  Firmicu… Clos… Clos… Rumin… Acet…             1988             2981
    ##  6 Bacteria  Proteob… Alph… Rhod… Aceto… Acet…               66               86
    ##  7 Bacteria  Firmicu… Clos… Hala… Halob… Acet…             1036             1373
    ##  8 Bacteria  Teneric… Moll… Acho… Achol… Acho…              779             1269
    ##  9 Bacteria  Proteob… Beta… Burk… Alcal… Achr…              192              298
    ## 10 Bacteria  Firmicu… Nega… Sele… Acida… Acid…            50181            39909
    ## # … with 2 more variables: `ShotgunWGS-ControlPig3GutMicrobiome-Day14` <dbl>,
    ## #   `ShotgunWGS-TomatoPig14GutMicrobiome-Day7` <dbl>

``` r
dim(Genus.Counts.Filt)
```

    ## [1] 895  66

``` r
# final genera list (after filtering for zeros)
final_genera[1:10,]
```

    ##  [1] "Abiotrophia"     "Acaryochloris"   "Acetivibrio"     "Acetobacter"    
    ##  [5] "Acetohalobium"   "Acholeplasma"    "Achromobacter"   "Acidaminococcus"
    ##  [9] "Acidilobus"      "Acidimicrobium"

``` r
# how many final genera do we have?
dim(final_genera)
```

    ## [1] 755   1

``` r
# join to create a df with genera in rows, samples in columns
# filtered for genera used in this analysis
genera_counts_foraldex <- inner_join(final_genera, Genus.Counts.Filt,
                                     by = "genus")

dim(genera_counts_foraldex)
```

    ## [1] 755  66

``` r
# remove non-necessary metadata
genera_counts_foraldex <- genera_counts_foraldex[,-c(2:6)]

genera_counts_foraldex[1:10, 1:4]
```

    ##              genus ShotgunWGS-ControlPig6GutMicrobiome-Day14
    ## 1      Abiotrophia                                      5067
    ## 2    Acaryochloris                                       271
    ## 3      Acetivibrio                                      1988
    ## 4      Acetobacter                                        66
    ## 5    Acetohalobium                                      1036
    ## 6     Acholeplasma                                       779
    ## 7    Achromobacter                                       192
    ## 8  Acidaminococcus                                     50181
    ## 9       Acidilobus                                        10
    ## 10  Acidimicrobium                                        59
    ##    ShotgunWGS-ControlPig8GutMicrobiome-Day0
    ## 1                                      5661
    ## 2                                       416
    ## 3                                      2981
    ## 4                                        86
    ## 5                                      1373
    ## 6                                      1269
    ## 7                                       298
    ## 8                                     39909
    ## 9                                        20
    ## 10                                      126
    ##    ShotgunWGS-ControlPig3GutMicrobiome-Day14
    ## 1                                       4117
    ## 2                                        267
    ## 3                                       2071
    ## 4                                         60
    ## 5                                       1015
    ## 6                                        817
    ## 7                                        197
    ## 8                                      31994
    ## 9                                         25
    ## 10                                        81

``` r
# add genera as rownames
rownames(genera_counts_foraldex) <- genera_counts_foraldex$genus

# remove genera as column for cleaner data
genera_counts_foraldex <- genera_counts_foraldex %>%
  select(-genus)
```

### Subset by Time

#### Day 0

Look at the effect of diet at day 0.

``` r
# subset day 0 only
Day0.Counts.Genera.filt <- genera_counts_foraldex %>% 
  select(ends_with("Day0"))
```

ALDEx2 function needs a factor of variables

``` r
# order alphabetically so making the meta data vector is easier
Day0.Counts.Genera.filt <- Day0.Counts.Genera.filt[order(colnames(Day0.Counts.Genera.filt))]

Diets.Day0.Genera <- as.vector(c(rep("Control", times=10), rep("Tomato", times=10)))

# check and make sure it came out right
Diets.Day0.Genera
```

    ##  [1] "Control" "Control" "Control" "Control" "Control" "Control" "Control"
    ##  [8] "Control" "Control" "Control" "Tomato"  "Tomato"  "Tomato"  "Tomato" 
    ## [15] "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"

Run t-test

``` r
filt.Genera.Day0.ByDiet.aldex <- aldex(Day0.Counts.Genera.filt, 
                                       Diets.Day0.Genera, 
                                       mc.samples = 1000, 
                                       test = "t", 
                                       effect = TRUE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.ttest: doing t-test

    ## aldex.effect: calculating effect sizes

``` r
filt.Genera.Day0.ByDiet.aldex <- 
  filt.Genera.Day0.ByDiet.aldex[order(filt.Genera.Day0.ByDiet.aldex$we.eBH, 
                                      decreasing = FALSE),]

kable(head(filt.Genera.Day0.ByDiet.aldex))
```

|                        |    rab.all | rab.win.Control | rab.win.Tomato |   diff.btw |  diff.win |     effect |   overlap |     we.ep |    we.eBH |     wi.ep |    wi.eBH |
|:-----------------------|-----------:|----------------:|---------------:|-----------:|----------:|-----------:|----------:|----------:|----------:|----------:|----------:|
| Methanoplanus          |  0.3201847 |       0.4587166 |      0.1585091 | -0.2744457 | 0.2929463 | -0.8900506 | 0.1761648 | 0.0161507 | 0.6571223 | 0.0206049 | 0.6391049 |
| Herbaspirillum         | -0.5564932 |      -0.7855297 |     -0.3447862 |  0.4894062 | 0.5715672 |  0.8331084 | 0.1700000 | 0.0108595 | 0.6613759 | 0.0156222 | 0.6367983 |
| Caenorhabditis         | -0.8681507 |      -0.7231369 |     -1.0313465 | -0.3344422 | 0.4358249 | -0.7236616 | 0.1917616 | 0.0233795 | 0.6654581 | 0.0326584 | 0.6445105 |
| Gallionella            | -0.8191944 |      -0.9806255 |     -0.6565452 |  0.3378006 | 0.4275946 |  0.7312681 | 0.1807638 | 0.0233197 | 0.6656534 | 0.0259145 | 0.6331800 |
| Epsilon15-like viruses | -5.4155796 |      -6.3291178 |     -4.7770370 |  1.5282103 | 1.9392933 |  0.7162934 | 0.1777644 | 0.0342860 | 0.6681890 | 0.0276429 | 0.6289355 |
| Collinsella            |  6.4097209 |       5.9122154 |      6.8773284 |  0.9055678 | 0.8808264 |  0.9277446 | 0.1979604 | 0.0098821 | 0.6683595 | 0.0208505 | 0.6554576 |

Create a histogram of pvalues of `we.eBH`

``` r
hist(filt.Genera.Day0.ByDiet.aldex$we.eBH,
     breaks = 20,
     main = "Histogram of p-values on the effect of diet at day 0 on genera",
     xlab = "Benjamini Hochberg corrected p-value (we.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-87-1.png)<!-- -->

`we.eBH` is the Benjamini-Hochberg corrected p-value, no significantly
different genera at day 0.

#### Day 7

Look at the effect of diet at day 7.

``` r
# subset day 7 only
Day7.Counts.Genera.filt <- genera_counts_foraldex %>% 
  select(ends_with("Day7"))
```

ALDEx2 function needs a factor of variables

``` r
# order alphabetically so making the meta data vector is easier
Day7.Counts.Genera.filt <- Day7.Counts.Genera.filt[order(colnames(Day7.Counts.Genera.filt))]

Diets.Day7.Genera <- as.vector(c(rep("Control", times=10), rep("Tomato", times=10)))

# check and make sure it came out right
Diets.Day7.Genera
```

    ##  [1] "Control" "Control" "Control" "Control" "Control" "Control" "Control"
    ##  [8] "Control" "Control" "Control" "Tomato"  "Tomato"  "Tomato"  "Tomato" 
    ## [15] "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"

Run t-test

``` r
filt.Genera.Day7.ByDiet.aldex <- aldex(Day7.Counts.Genera.filt, 
                                       Diets.Day7.Genera, 
                                       mc.samples = 1000, 
                                       test = "t", 
                                       effect = TRUE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.ttest: doing t-test

    ## aldex.effect: calculating effect sizes

``` r
filt.Genera.Day7.ByDiet.aldex <- 
  filt.Genera.Day7.ByDiet.aldex[order(filt.Genera.Day7.ByDiet.aldex$we.eBH, 
                                      decreasing = FALSE),]

kable(head(filt.Genera.Day7.ByDiet.aldex))
```

|                                      |    rab.all | rab.win.Control | rab.win.Tomato |  diff.btw |  diff.win |    effect |   overlap |     we.ep |    we.eBH |     wi.ep |    wi.eBH |
|:-------------------------------------|-----------:|----------------:|---------------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
| unclassified (derived from Bacteria) |  2.5026985 |       1.9874193 |      3.0144669 | 0.8894342 | 0.5458126 | 1.6215552 | 0.0466000 | 0.0000347 | 0.0256786 | 0.0001885 | 0.0816149 |
| Staphylococcus                       |  4.3981480 |       4.2392110 |      4.6647900 | 0.3930480 | 0.2840870 | 1.3716187 | 0.0507898 | 0.0003334 | 0.1056045 | 0.0002354 | 0.0877397 |
| Alphatorquevirus                     | -2.7795012 |      -3.9782818 |     -2.1901574 | 1.8729947 | 1.5103759 | 1.2115892 | 0.0950000 | 0.0013080 | 0.2092870 | 0.0017119 | 0.2325936 |
| Lambda-like viruses                  | -0.5334489 |      -1.3489239 |      0.2540438 | 1.5967109 | 1.1783139 | 1.1878787 | 0.1027794 | 0.0013951 | 0.2330323 | 0.0017119 | 0.2490163 |
| Clavibacter                          |  0.7748387 |       0.5505235 |      1.0398215 | 0.6786631 | 0.6909331 | 0.9901394 | 0.1345731 | 0.0034204 | 0.3646270 | 0.0059299 | 0.4272165 |
| Kluyveromyces                        | -2.7486549 |      -3.1193205 |     -2.3047942 | 0.8797267 | 0.9542796 | 0.8722283 | 0.1534000 | 0.0137741 | 0.5092141 | 0.0189920 | 0.5001200 |

One genera was significantly different by diet at day 7 - unclassified
(derived from bacteria), padj = 0.025

``` r
hist(filt.Genera.Day7.ByDiet.aldex$we.eBH,
     breaks = 20,
     main = "Histogram of p-values on the effect of diet at day 7 on genera",
     xlab = "Benjamini Hochberg corrected p-value (we.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-92-1.png)<!-- -->

What is the directionality of the change?

``` r
filt.Genera.Day7.ByDiet.aldex %>%
  select(rab.win.Control, rab.win.Tomato, we.eBH) %>%
  filter(we.eBH <= 0.05)
```

    ##                                      rab.win.Control rab.win.Tomato     we.eBH
    ## unclassified (derived from Bacteria)        1.987419       3.014467 0.02567862

Unclassified (derived from Bacteria) is higher in Tomato than Control.

#### Day 14

Look at the effect of diet on day 14.

``` r
# subset day 14 only
Day14.Counts.Genera.filt <- genera_counts_foraldex %>% 
  select(ends_with("Day14"))
```

ALDEx2 function needs a factor of variables

``` r
# order alphabetically so making the meta data vector is easier
Day14.Counts.Genera.filt <- Day14.Counts.Genera.filt[order(colnames(Day14.Counts.Genera.filt))]

Diets.Day14.Genera <- as.vector(c(rep("Control", times=10), rep("Tomato", times=10)))

# check and make sure it came out right
Diets.Day14.Genera
```

    ##  [1] "Control" "Control" "Control" "Control" "Control" "Control" "Control"
    ##  [8] "Control" "Control" "Control" "Tomato"  "Tomato"  "Tomato"  "Tomato" 
    ## [15] "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"

Run t-test

``` r
filt.Genera.Day14.ByDiet.aldex <- aldex(Day14.Counts.Genera.filt, 
                                       Diets.Day14.Genera, 
                                       mc.samples = 1000, 
                                       test = "t", 
                                       effect = TRUE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.ttest: doing t-test

    ## aldex.effect: calculating effect sizes

``` r
filt.Genera.Day14.ByDiet.aldex <- 
  filt.Genera.Day14.ByDiet.aldex[order(filt.Genera.Day14.ByDiet.aldex$we.eBH, 
                                      decreasing = FALSE),]

kable(head(filt.Genera.Day14.ByDiet.aldex))
```

|                                      |    rab.all | rab.win.Control | rab.win.Tomato |  diff.btw |  diff.win |   effect |  overlap |    we.ep |    we.eBH |     wi.ep |    wi.eBH |
|:-------------------------------------|-----------:|----------------:|---------------:|----------:|----------:|---------:|---------:|---------:|----------:|----------:|----------:|
| Lambda-like viruses                  | -0.6765788 |       -2.014760 |      1.1151072 | 3.1367093 | 0.5031233 | 6.152552 | 0.000014 | 0.00e+00 | 0.0000000 | 0.0000108 | 0.0017403 |
| Staphylococcus                       |  4.4868458 |        4.203031 |      4.8287360 | 0.6626948 | 0.2017804 | 3.372094 | 0.000014 | 0.00e+00 | 0.0000020 | 0.0000108 | 0.0017403 |
| Alphatorquevirus                     | -2.7478328 |       -4.775638 |     -0.9952775 | 3.7616825 | 1.1479595 | 3.213701 | 0.000014 | 1.00e-07 | 0.0000166 | 0.0000108 | 0.0017403 |
| unclassified (derived from Bacteria) |  2.1148119 |        1.507804 |      2.9099037 | 1.3704930 | 0.5614568 | 2.482695 | 0.000014 | 3.00e-07 | 0.0000489 | 0.0000108 | 0.0017403 |
| Loa                                  | -3.2560289 |       -4.570576 |     -2.2214430 | 2.3520024 | 1.2624078 | 1.786604 | 0.030000 | 6.60e-05 | 0.0053620 | 0.0001009 | 0.0075989 |
| Plasmodium                           | -0.4456836 |       -1.018278 |     -0.0948831 | 0.9253498 | 0.4646515 | 1.835934 | 0.060188 | 8.46e-05 | 0.0068115 | 0.0004468 | 0.0195314 |

``` r
hist(filt.Genera.Day14.ByDiet.aldex$we.eBH,
     breaks = 20,
     main = "Histogram of p-values on the effect of diet at day 14 on genera",
     xlab = "Benjamini Hochberg corrected p-value (we.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-98-1.png)<!-- -->

How many significant genera are there?

``` r
filt.Day14.Genera.aldex.sig <- filt.Genera.Day14.ByDiet.aldex[which(filt.Genera.Day14.ByDiet.aldex$we.eBH<0.05),]

length(rownames(filt.Day14.Genera.aldex.sig))
```

    ## [1] 14

Which genera are they?

``` r
sig_day14_genera_aldex2 <- as.data.frame(cbind(rownames(filt.Day14.Genera.aldex.sig),
                                 filt.Day14.Genera.aldex.sig$we.eBH))

sig_day14_genera_aldex2 <- sig_day14_genera_aldex2 %>%
  rename(Genera = V1,
         we.eBH_pvalue = V2)

sig_day14_genera_aldex2
```

    ##                                  Genera        we.eBH_pvalue
    ## 1                   Lambda-like viruses 3.90702968822257e-08
    ## 2                        Staphylococcus 1.98833990718406e-06
    ## 3                      Alphatorquevirus 1.66493484376377e-05
    ## 4  unclassified (derived from Bacteria) 4.88751914488331e-05
    ## 5                                   Loa  0.00536199352280373
    ## 6                            Plasmodium  0.00681151568341803
    ## 7                     Propionibacterium  0.00938576734150738
    ## 8                         Saccharomyces   0.0162702517857512
    ## 9                      Stenotrophomonas   0.0215525889312409
    ## 10                           Malassezia   0.0222884879030356
    ## 11                          Roseiflexus   0.0224147224012945
    ## 12                               Brugia   0.0315547210562422
    ## 13                        Streptococcus   0.0316064329890332
    ## 14                      Vanderwaltozyma   0.0355364448144129

-   Lambda-like viruses
-   Staphylococcus
-   Alphatorquervirus
-   unclassified (derived from Bacteria)
-   Loa
-   Plasmodium
-   Propionibacterium
-   Saccharomyces, Stenotrophomonas
-   Malassezia
-   Roseiflexus
-   Brugia
-   Strepococcus
-   Vanderwaltozyma

What is the directionality of the change?

``` r
filt.Genera.Day14.ByDiet.aldex %>%
  select(rab.win.Control, rab.win.Tomato, we.eBH) %>%
  filter(we.eBH <= 0.05)
```

    ##                                      rab.win.Control rab.win.Tomato
    ## Lambda-like viruses                      -2.01475973      1.1151072
    ## Staphylococcus                            4.20303084      4.8287360
    ## Alphatorquevirus                         -4.77563841     -0.9952775
    ## unclassified (derived from Bacteria)      1.50780388      2.9099037
    ## Loa                                      -4.57057605     -2.2214430
    ## Plasmodium                               -1.01827823     -0.0948831
    ## Propionibacterium                         1.40929760      2.0214784
    ## Saccharomyces                            -1.92007049      0.3880179
    ## Stenotrophomonas                          0.01556489      0.4808232
    ## Malassezia                               -3.61348687     -1.7990132
    ## Roseiflexus                               2.61037473      2.3131767
    ## Brugia                                   -3.16504650     -1.6096471
    ## Streptococcus                            10.54449623      8.0218160
    ## Vanderwaltozyma                          -4.08112899     -1.1653644
    ##                                            we.eBH
    ## Lambda-like viruses                  3.907030e-08
    ## Staphylococcus                       1.988340e-06
    ## Alphatorquevirus                     1.664935e-05
    ## unclassified (derived from Bacteria) 4.887519e-05
    ## Loa                                  5.361994e-03
    ## Plasmodium                           6.811516e-03
    ## Propionibacterium                    9.385767e-03
    ## Saccharomyces                        1.627025e-02
    ## Stenotrophomonas                     2.155259e-02
    ## Malassezia                           2.228849e-02
    ## Roseiflexus                          2.241472e-02
    ## Brugia                               3.155472e-02
    ## Streptococcus                        3.160643e-02
    ## Vanderwaltozyma                      3.553644e-02

All significantly different genera are higher in tomato as compared to
control.

### Subset by diet

#### Control

``` r
# subset control only samples across all time points, should be n=30
Control.Counts.Genera.filt <- genera_counts_foraldex %>% 
  select(contains("Control"))

dim(Control.Counts.Genera.filt)
```

    ## [1] 755  30

ALDEx2 function needs a factor of variables

``` r
# results in pigs at different time points being grouped together
Control.Counts.Genera.filt <- Control.Counts.Genera.filt[order(colnames(Control.Counts.Genera.filt))]

# then time point by "alphabetical" where 14 comes before 7
# ex, first few are Pig 10 Day 0, Pig 10 Day 14, Pig 10 Day 7, Pig 1 Day 0, Pig 1 Day 14, etc
TimePoints.Control.Genera <- as.vector(rep(c("Day0", "Day14", "Day7"), times=10))

# check and make sure it looks right
TimePoints.Control.Genera
```

    ##  [1] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [10] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [19] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [28] "Day0"  "Day14" "Day7"

More than two conditions this time, use the ANOVA-like test, Kruskal
Wallis

``` r
filt.Genera.Control.ByTime.aldex <- aldex(Control.Counts.Genera.filt, 
                                          TimePoints.Control.Genera, 
                                          mc.samples = 1000, 
                                          test = "kw", 
                                          effect = FALSE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.glm: doing Kruskal-Wallace and glm test (ANOVA-like)

    ## operating in serial mode

We are looking at `glm.eBH` for the BH corrected ANOVA pvalue

``` r
filt.Genera.Control.ByTime.aldex <- 
  filt.Genera.Control.ByTime.aldex[order(filt.Genera.Control.ByTime.aldex$glm.eBH, 
                                         decreasing = FALSE),]

kable(head(filt.Genera.Control.ByTime.aldex))
```

|                     |     kw.ep |    kw.eBH |    glm.ep |   glm.eBH |
|:--------------------|----------:|----------:|----------:|----------:|
| Oribacterium        | 0.0005354 | 0.1104829 | 0.0000000 | 0.0000158 |
| Streptococcus       | 0.0000760 | 0.0558923 | 0.0000000 | 0.0000159 |
| Lactococcus         | 0.0004159 | 0.1045718 | 0.0000065 | 0.0015823 |
| Granulicatella      | 0.0007024 | 0.1276212 | 0.0000951 | 0.0149421 |
| T4-like viruses     | 0.0056658 | 0.3088761 | 0.0011659 | 0.0822780 |
| Schizosaccharomyces | 0.0094127 | 0.3725284 | 0.0026917 | 0.1301247 |

``` r
hist(filt.Genera.Control.ByTime.aldex$glm.eBH,
     breaks = 20,
     main = "Histogram of p-values on the effect of time within the control diet on genera",
     xlab = "Benjamini Hochberg corrected p-value (glm.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-106-1.png)<!-- -->

How many significantly different genera are there?

``` r
filt.Genera.Control.ByTime.aldex.sig <- 
  filt.Genera.Control.ByTime.aldex[which(filt.Genera.Control.ByTime.aldex$glm.eBH<0.05),]

length(rownames(filt.Genera.Control.ByTime.aldex.sig))
```

    ## [1] 4

4 sig genera

Which genera are they?

``` r
sig_control_genera_aldex2 <- as.data.frame(cbind(rownames(filt.Genera.Control.ByTime.aldex.sig),
                                 filt.Genera.Control.ByTime.aldex.sig$glm.eBH))

sig_control_genera_aldex2 <- sig_control_genera_aldex2 %>%
  rename(Genera = V1,
         glm.eBH_pval = V2)

sig_control_genera_aldex2
```

    ##           Genera         glm.eBH_pval
    ## 1   Oribacterium 1.58068029273393e-05
    ## 2  Streptococcus 1.58746766478819e-05
    ## 3    Lactococcus  0.00158226394942239
    ## 4 Granulicatella   0.0149420529335598

-   Oribacterium
-   Streptococcus
-   Lactococcus
-   Granulicatella

#### Tomato

``` r
# subset tomato only samples across all time points, should be n=30
Tomato.Counts.Genera.filt <- genera_counts_foraldex %>% 
  select(contains("Tomato"))
```

ALDEx2 function needs a factor of variables

``` r
# results in pigs at different time points being grouped together
Tomato.Counts.Genera.filt <- Tomato.Counts.Genera.filt[order(colnames(Tomato.Counts.Genera.filt))]

# then time point by "alphabetical" where 14 comes before 7
# ex, first few are Pig 10 Day 0, Pig 10 Day 14, Pig 10 Day 7, Pig 1 Day 0, Pig 1 Day 14, etc
TimePoints.Tomato.Genera <- as.vector(rep(c("Day0", "Day14", "Day7"), times=10))

# check and make sure it looks right
TimePoints.Tomato.Genera
```

    ##  [1] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [10] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [19] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [28] "Day0"  "Day14" "Day7"

More than two conditions this time, use the ANOVA-like test

``` r
filt.Genera.Tomato.ByTime.aldex <- aldex(Tomato.Counts.Genera.filt, 
                                          TimePoints.Tomato.Genera, 
                                          mc.samples = 1000, 
                                          test = "kw", 
                                          effect = FALSE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.glm: doing Kruskal-Wallace and glm test (ANOVA-like)

    ## operating in serial mode

We are looking at `glm.eBH` for the BH corrected ANOVA pvalue

``` r
filt.Genera.Tomato.ByTime.aldex <- 
  filt.Genera.Tomato.ByTime.aldex[order(filt.Genera.Tomato.ByTime.aldex$glm.eBH, 
                                         decreasing = FALSE),]

kable(head(filt.Genera.Tomato.ByTime.aldex))
```

|                                      |     kw.ep |    kw.eBH |    glm.ep |   glm.eBH |
|:-------------------------------------|----------:|----------:|----------:|----------:|
| Staphylococcus                       | 0.0000565 | 0.0342115 | 0.0000000 | 0.0000001 |
| Alphatorquevirus                     | 0.0001143 | 0.0390495 | 0.0000003 | 0.0000802 |
| Lambda-like viruses                  | 0.0002278 | 0.0582011 | 0.0000017 | 0.0004113 |
| unclassified (derived from Bacteria) | 0.0019396 | 0.2946385 | 0.0000066 | 0.0012508 |
| Streptococcus                        | 0.0025595 | 0.3275011 | 0.0014981 | 0.1707533 |
| Crocosphaera                         | 0.0117640 | 0.4694176 | 0.0104829 | 0.3025469 |

``` r
hist(filt.Genera.Tomato.ByTime.aldex$glm.eBH,
     breaks = 20,
     main = "Histogram of p-values on the effect of time within the tomato diet on genera",
     xlab = "Benjamini Hochberg corrected p-value (glm.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-113-1.png)<!-- -->

How many significantly different genera are there?

``` r
filt.Genera.Tomato.ByTime.aldex.sig <- 
  filt.Genera.Tomato.ByTime.aldex[which(filt.Genera.Tomato.ByTime.aldex$glm.eBH<0.05),]

length(rownames(filt.Genera.Tomato.ByTime.aldex.sig))
```

    ## [1] 4

4 sig genera

Which genera are they?

``` r
sig_tomato_genera_aldex2 <- as.data.frame(cbind(rownames(filt.Genera.Tomato.ByTime.aldex.sig),
                                 filt.Genera.Tomato.ByTime.aldex.sig$glm.eBH))

sig_tomato_genera_aldex2 <- sig_tomato_genera_aldex2 %>%
  rename(Genera = V1,
         glm.eBH_pval = V2)

sig_tomato_genera_aldex2
```

    ##                                 Genera         glm.eBH_pval
    ## 1                       Staphylococcus 1.09222660374843e-07
    ## 2                     Alphatorquevirus 8.01595357034658e-05
    ## 3                  Lambda-like viruses 0.000411303309560209
    ## 4 unclassified (derived from Bacteria)  0.00125079010124612

-   Staphylococcus
-   Alphatorquevirus
-   Lambda-like viruses
-   unclassified (derived from Bacteria)

Any overlap between sig differences at day 14 and by diet?

Control over time and day 14 overlap

``` r
intersect(sig_day14_genera_aldex2$Genera, sig_control_genera_aldex2$Genera)
```

    ## [1] "Streptococcus"

Streptococcus

Tomato over time and day 14 overlap

``` r
intersect(sig_day14_genera_aldex2$Genera, sig_tomato_genera_aldex2$Genera)
```

    ## [1] "Lambda-like viruses"                 
    ## [2] "Staphylococcus"                      
    ## [3] "Alphatorquevirus"                    
    ## [4] "unclassified (derived from Bacteria)"

-   Lambda-like viruses
-   Staphylococcus
-   Alphatorquevirus
-   unclassified (derived from Bacteria)

# Phyla-level annotation

Read in phyla level data, annotated from MG-RAST. In “Phyla” tab of
Supplementary Information.

``` r
Phyla.Counts <- read_excel("Goggans_etal_2021_tomato_pig_microbiome_WGS.xlsx",
                                       sheet = "TableS3.Phyla")

str(Phyla.Counts)
```

    ## tibble [60 × 62] (S3: tbl_df/tbl/data.frame)
    ##  $ domain                                    : chr [1:60] "Bacteria" "Bacteria" "Eukaryota" "Bacteria" ...
    ##  $ phylum                                    : chr [1:60] "Acidobacteria" "Actinobacteria" "Apicomplexa" "Aquificae" ...
    ##  $ ShotgunWGS-ControlPig6GutMicrobiome-Day14 : num [1:60] 2874 186789 231 1953 368 ...
    ##  $ ShotgunWGS-ControlPig8GutMicrobiome-Day0  : num [1:60] 3717 277130 384 2254 992 ...
    ##  $ ShotgunWGS-ControlPig3GutMicrobiome-Day14 : num [1:60] 2663 126155 190 1642 386 ...
    ##  $ ShotgunWGS-TomatoPig14GutMicrobiome-Day7  : num [1:60] 880 39557 168 647 211 ...
    ##  $ ShotgunWGS-ControlPig5GutMicrobiome-Day7  : num [1:60] 2016 142345 171 1418 400 ...
    ##  $ ShotgunWGS-TomatoPig18GutMicrobiome-Day7  : num [1:60] 1377 181295 101 658 201 ...
    ##  $ ShotgunWGS-TomatoPig16GutMicrobiome-Day7  : num [1:60] 1570 58263 224 892 340 ...
    ##  $ ShotgunWGS-ControlPig10GutMicrobiome-Day7 : num [1:60] 1298 109273 287 900 211 ...
    ##  $ ShotgunWGS-ControlPig2GutMicrobiome-Day0  : num [1:60] 3114 159425 529 2095 764 ...
    ##  $ ShotgunWGS-TomatoPig18GutMicrobiome-Day0  : num [1:60] 2604 168472 170 1422 393 ...
    ##  $ ShotgunWGS-ControlPig10GutMicrobiome-Day0 : num [1:60] 3118 163425 231 1192 426 ...
    ##  $ ShotgunWGS-ControlPig7GutMicrobiome-Day0  : num [1:60] 2796 70967 389 1137 605 ...
    ##  $ ShotgunWGS-ControlPig8GutMicrobiome-Day14 : num [1:60] 2222 91465 325 1561 364 ...
    ##  $ ShotgunWGS-TomatoPig11GutMicrobiome-Day0  : num [1:60] 2136 68481 402 1377 691 ...
    ##  $ ShotgunWGS-TomatoPig19GutMicrobiome-Day0  : num [1:60] 2017 207693 143 1265 361 ...
    ##  $ ShotgunWGS-TomatoPig17GutMicrobiome-Day14 : num [1:60] 836 26050 147 633 212 ...
    ##  $ ShotgunWGS-ControlPig9GutMicrobiome-Day14 : num [1:60] 2612 172091 181 1645 393 ...
    ##  $ ShotgunWGS-ControlPig10GutMicrobiome-Day14: num [1:60] 2136 122681 250 1536 445 ...
    ##  $ ShotgunWGS-TomatoPig19GutMicrobiome-Day7  : num [1:60] 1090 78218 168 774 304 ...
    ##  $ ShotgunWGS-ControlPig5GutMicrobiome-Day14 : num [1:60] 2693 263950 266 1713 577 ...
    ##  $ ShotgunWGS-ControlPig2GutMicrobiome-Day7  : num [1:60] 3420 101192 582 2369 766 ...
    ##  $ ShotgunWGS-ControlPig6GutMicrobiome-Day7  : num [1:60] 2216 159323 115 1383 303 ...
    ##  $ ShotgunWGS-TomatoPig12GutMicrobiome-Day0  : num [1:60] 2146 78205 221 1265 390 ...
    ##  $ ShotgunWGS-TomatoPig14GutMicrobiome-Day0  : num [1:60] 732 77377 292 585 223 ...
    ##  $ ShotgunWGS-ControlPig7GutMicrobiome-Day14 : num [1:60] 2079 142139 322 1335 392 ...
    ##  $ ShotgunWGS-TomatoPig11GutMicrobiome-Day14 : num [1:60] 570 25927 180 425 270 ...
    ##  $ ShotgunWGS-TomatoPig20GutMicrobiome-Day0  : num [1:60] 2472 82091 415 1534 647 ...
    ##  $ ShotgunWGS-ControlPig9GutMicrobiome-Day0  : num [1:60] 1607 88397 432 1085 423 ...
    ##  $ ShotgunWGS-TomatoPig11GutMicrobiome-Day7  : num [1:60] 278 17451 96 150 107 ...
    ##  $ ShotgunWGS-TomatoPig13GutMicrobiome-Day7  : num [1:60] 1100 56205 157 984 306 ...
    ##  $ ShotgunWGS-TomatoPig17GutMicrobiome-Day0  : num [1:60] 1562 74553 171 780 238 ...
    ##  $ ShotgunWGS-TomatoPig19GutMicrobiome-Day14 : num [1:60] 765 47957 182 551 237 ...
    ##  $ ShotgunWGS-TomatoPig13GutMicrobiome-Day0  : num [1:60] 2182 124473 280 1483 476 ...
    ##  $ ShotgunWGS-ControlPig2GutMicrobiome-Day14 : num [1:60] 3329 116448 325 2149 703 ...
    ##  $ ShotgunWGS-ControlPig1GutMicrobiome-Day7  : num [1:60] 1920 55849 156 1234 408 ...
    ##  $ ShotgunWGS-TomatoPig15GutMicrobiome-Day7  : num [1:60] 757 38904 204 583 317 ...
    ##  $ ShotgunWGS-TomatoPig15GutMicrobiome-Day0  : num [1:60] 2037 120272 320 1399 561 ...
    ##  $ ShotgunWGS-TomatoPig12GutMicrobiome-Day7  : num [1:60] 1279 87121 215 917 409 ...
    ##  $ ShotgunWGS-TomatoPig14GutMicrobiome-Day14 : num [1:60] 583 36948 69 444 102 ...
    ##  $ ShotgunWGS-TomatoPig20GutMicrobiome-Day14 : num [1:60] 496 29179 99 374 142 ...
    ##  $ ShotgunWGS-ControlPig1GutMicrobiome-Day0  : num [1:60] 2963 90535 278 1596 631 ...
    ##  $ ShotgunWGS-ControlPig4GutMicrobiome-Day14 : num [1:60] 2548 133556 181 1734 432 ...
    ##  $ ShotgunWGS-ControlPig6GutMicrobiome-Day0  : num [1:60] 2269 127508 314 1058 413 ...
    ##  $ ShotgunWGS-TomatoPig16GutMicrobiome-Day0  : num [1:60] 1935 110140 207 1018 425 ...
    ##  $ ShotgunWGS-TomatoPig16GutMicrobiome-Day14 : num [1:60] 817 33981 133 443 187 ...
    ##  $ ShotgunWGS-TomatoPig18GutMicrobiome-Day14 : num [1:60] 705 92977 122 507 148 ...
    ##  $ ShotgunWGS-ControlPig7GutMicrobiome-Day7  : num [1:60] 1131 69602 297 628 290 ...
    ##  $ ShotgunWGS-ControlPig4GutMicrobiome-Day7  : num [1:60] 1298 112714 203 983 325 ...
    ##  $ ShotgunWGS-TomatoPig13GutMicrobiome-Day14 : num [1:60] 566 36447 72 514 125 ...
    ##  $ ShotgunWGS-ControlPig8GutMicrobiome-Day7  : num [1:60] 2173 159187 378 1311 361 ...
    ##  $ ShotgunWGS-TomatoPig15GutMicrobiome-Day14 : num [1:60] 1186 49134 150 858 249 ...
    ##  $ ShotgunWGS-TomatoPig12GutMicrobiome-Day14 : num [1:60] 1122 64744 254 1030 427 ...
    ##  $ ShotgunWGS-TomatoPig20GutMicrobiome-Day7  : num [1:60] 1109 97728 149 670 211 ...
    ##  $ ShotgunWGS-ControlPig1GutMicrobiome-Day14 : num [1:60] 2350 83993 210 1719 446 ...
    ##  $ ShotgunWGS-ControlPig3GutMicrobiome-Day0  : num [1:60] 3314 428097 206 2366 519 ...
    ##  $ ShotgunWGS-ControlPig5GutMicrobiome-Day0  : num [1:60] 2998 242356 283 1895 758 ...
    ##  $ ShotgunWGS-ControlPig4GutMicrobiome-Day0  : num [1:60] 3042 223010 351 1777 685 ...
    ##  $ ShotgunWGS-ControlPig9GutMicrobiome-Day7  : num [1:60] 499 68424 784 329 171 ...
    ##  $ ShotgunWGS-ControlPig3GutMicrobiome-Day7  : num [1:60] 2620 340300 165 1993 484 ...
    ##  $ ShotgunWGS-TomatoPig17GutMicrobime-Day7   : num [1:60] 1340 71395 159 648 270 ...

## Data filtering

### Remove inplausible phyla

These phyla are not plausibly found in a rectal swab of a pig, and were
incorrectly annotated, so we are removing them.

``` r
Phyla.Counts.Filt <- Phyla.Counts %>%
  filter(phylum != "Chordata" , phylum != "Arthropoda" , phylum != "Cnidaria" , 
         phylum != "Porifera" , phylum != "Echinodermata", phylum != "Streptophyta",
         phylum != "Platyhelminthes")
```

Transpose.

``` r
Phyla.Counts.Filt.t <- as.tibble(t(Phyla.Counts.Filt))

# make phyla colnames
colnames(Phyla.Counts.Filt.t) <- Phyla.Counts.Filt.t[2,]

# remove domain, phylum rows
Phyla.Counts.Filt.t <- Phyla.Counts.Filt.t[3:62,]

# convert character to numeric
Phyla.Counts.Filt.t <- as.data.frame(apply((Phyla.Counts.Filt.t), 2, as.numeric))

str(Phyla.Counts.Filt.t[,1:5])
```

    ## 'data.frame':    60 obs. of  5 variables:
    ##  $ Acidobacteria : num  2874 3717 2663 880 2016 ...
    ##  $ Actinobacteria: num  186789 277130 126155 39557 142345 ...
    ##  $ Apicomplexa   : num  231 384 190 168 171 101 224 287 529 170 ...
    ##  $ Aquificae     : num  1953 2254 1642 647 1418 ...
    ##  $ Ascomycota    : num  1491 2196 1281 672 1178 ...

``` r
# add back sample names as column
Phyla.Counts.Filt.t <- Phyla.Counts.Filt.t %>%
  mutate(Sample_Name = AllSamples.Metadata$Sample_Name)

# move Sample_Name to first column
Phyla.Counts.Filt.t <- Phyla.Counts.Filt.t %>%
  relocate(Sample_Name)

kable(head(Phyla.Counts.Filt.t))
```

| Sample\_Name                              | Acidobacteria | Actinobacteria | Apicomplexa | Aquificae | Ascomycota | Bacillariophyta | Bacteroidetes | Basidiomycota | Blastocladiomycota | Candidatus Poribacteria | Chlamydiae | Chlorobi | Chloroflexi | Chlorophyta | Chromerida | Chrysiogenetes | Chytridiomycota | Crenarchaeota | Cyanobacteria | Deferribacteres | Deinococcus-Thermus | Dictyoglomi | Elusimicrobia | Euglenida | Euryarchaeota | Fibrobacteres | Firmicutes | Fusobacteria | Gemmatimonadetes | Glomeromycota | Hemichordata | Korarchaeota | Lentisphaerae | Microsporidia | Nanoarchaeota | Nematoda | Nitrospirae | Phaeophyceae | Placozoa | Planctomycetes | Proteobacteria | Spirochaetes | Synergistetes | Tenericutes | Thaumarchaeota | Thermotogae | Verrucomicrobia | Xanthophyceae | unclassified (derived from Bacteria) | unclassified (derived from Eukaryota) | unclassified (derived from Fungi) | unclassified (derived from Viruses) | unclassified (derived from other sequences) |
|:------------------------------------------|--------------:|---------------:|------------:|----------:|-----------:|----------------:|--------------:|--------------:|-------------------:|------------------------:|-----------:|---------:|------------:|------------:|-----------:|---------------:|----------------:|--------------:|--------------:|----------------:|--------------------:|------------:|--------------:|----------:|--------------:|--------------:|-----------:|-------------:|-----------------:|--------------:|-------------:|-------------:|--------------:|--------------:|--------------:|---------:|------------:|-------------:|---------:|---------------:|---------------:|-------------:|--------------:|------------:|---------------:|------------:|----------------:|--------------:|-------------------------------------:|--------------------------------------:|----------------------------------:|------------------------------------:|--------------------------------------------:|
| ShotgunWGS-ControlPig6GutMicrobiome-Day14 |          2874 |         186789 |         231 |      1953 |       1491 |             105 |       1424565 |           240 |                  0 |                      26 |        552 |     4889 |        7842 |         370 |          0 |            331 |               0 |           648 |          8838 |            1494 |                2481 |        1217 |           632 |         3 |         13175 |          4768 |    2059948 |        15350 |              211 |             0 |           26 |           75 |           765 |            49 |             4 |      178 |         551 |            0 |       68 |           1523 |         105309 |        11519 |          4453 |        2764 |             56 |        5014 |            3209 |             1 |                                 1197 |                                  1260 |                                 0 |                                1546 |                                          48 |
| ShotgunWGS-ControlPig8GutMicrobiome-Day0  |          3717 |         277130 |         384 |      2254 |       2196 |             184 |       1391417 |           405 |                  0 |                      49 |        829 |     6073 |        9612 |         571 |          0 |            416 |               0 |           902 |         13612 |            1994 |                3586 |        1554 |          1007 |         0 |         19176 |          6963 |    2223331 |        21242 |              340 |             0 |           18 |           97 |          1833 |            88 |             5 |      265 |         797 |            1 |       68 |           2632 |         154698 |        17463 |          7489 |        3731 |             80 |        7105 |            6282 |             1 |                                 1720 |                                  4189 |                                 1 |                                2626 |                                          33 |
| ShotgunWGS-ControlPig3GutMicrobiome-Day14 |          2663 |         126155 |         190 |      1642 |       1281 |             129 |       1260217 |           198 |                  0 |                      21 |        554 |     4469 |        7596 |         362 |          0 |            271 |               0 |           638 |         11276 |            1419 |                2321 |        1109 |           699 |         0 |         12790 |          4985 |    2266610 |        14356 |              241 |             0 |           16 |           52 |           774 |            56 |             4 |      146 |         563 |            2 |       31 |           1570 |         104879 |        10922 |          4148 |        2493 |             59 |        5052 |            3738 |             0 |                                 1148 |                                  1266 |                                 0 |                                2003 |                                          62 |
| ShotgunWGS-TomatoPig14GutMicrobiome-Day7  |           880 |          39557 |         168 |       647 |        672 |              65 |        415935 |           114 |                  0 |                      17 |        223 |     1565 |        2849 |         184 |          0 |            114 |               2 |           328 |          3426 |             543 |                 967 |         424 |           359 |         1 |          8299 |          1750 |     628580 |         5545 |               61 |             0 |            5 |           36 |           492 |            24 |             6 |      138 |         241 |            1 |       19 |            540 |          71783 |         5634 |          1870 |        1292 |             37 |        1998 |            1259 |             0 |                                 1161 |                                   662 |                                 6 |                                1010 |                                         115 |
| ShotgunWGS-ControlPig5GutMicrobiome-Day7  |          2016 |         142345 |         171 |      1418 |       1178 |             111 |        798569 |           182 |                  1 |                      22 |        505 |     3835 |        6249 |         388 |          0 |            273 |               0 |           652 |         10511 |            1166 |                2095 |         897 |           665 |         2 |         13289 |          4040 |    1919749 |        13260 |              211 |             0 |            8 |           66 |          1339 |            80 |             4 |      161 |         483 |            1 |       27 |           1674 |         114187 |        10757 |          4358 |        2309 |             55 |        4527 |            3442 |             1 |                                 1323 |                                  1369 |                                 0 |                                1475 |                                          13 |
| ShotgunWGS-TomatoPig18GutMicrobiome-Day7  |          1377 |         181295 |         101 |       658 |        799 |              64 |        690378 |           119 |                  0 |                      17 |        222 |     2198 |        3255 |         220 |          0 |             91 |               0 |           274 |          7460 |             445 |                1102 |         423 |           230 |         0 |          4970 |          2071 |     793943 |         5068 |              208 |             0 |            2 |           25 |           288 |            32 |             1 |      134 |         245 |            0 |       29 |           1096 |          69391 |         4728 |          1482 |         915 |             26 |        1710 |            2361 |             0 |                                 1431 |                                   518 |                                 0 |                                1267 |                                          59 |

Calculate relative abundance, and bind back to metadata.

``` r
Phyla.Counts.Filt.t.wtotal <- Phyla.Counts.Filt.t %>%
  mutate(Total.Counts = rowSums(Phyla.Counts.Filt.t[,2:ncol(Phyla.Counts.Filt.t)]))

dim(Phyla.Counts.Filt.t.wtotal)
```

    ## [1] 60 55

``` r
# create rel abund df
RelAbund.Phyla.Filt <- Phyla.Counts.Filt.t.wtotal[,2:54]/Phyla.Counts.Filt.t.wtotal$Total.Counts

# add back metadata
RelAbund.Phyla.Filt <- bind_cols(AllSamples.Metadata, RelAbund.Phyla.Filt)
```

### Counting missing data

``` r
# remove metadata
RelAbund.Phyla.Filt.nometadata <- RelAbund.Phyla.Filt %>%
  select_if(is.numeric) 

# create a list with the number of zeros for each genus
counting_zeros_phyla <- sapply(RelAbund.Phyla.Filt.nometadata, function(x){ (sum(x==0))})

# plot a histogram to look
counting_zeros_phyla_df <- as.data.frame(counting_zeros_phyla)

hist(counting_zeros_phyla_df$counting_zeros_phyla, 
     breaks = 61,
     main = "Histogram of Genera with Zero Relative Intensity",
     sub = "Starting at No Zeros",
     xlab = "Number of zero relative intensity values",
     ylab = "Frequency")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-122-1.png)<!-- -->

Big first bar is many phyla which have zero missing values.

``` r
# filter for any phyla with at least 1 missing value
counting_zeros_phyla_df_missingval <- counting_zeros_phyla_df %>%
  rownames_to_column(var = "rowname") %>%
  filter(counting_zeros_phyla > 0) %>%
  column_to_rownames(var = "rowname")

# how many genera have at least one missing value?
dim(counting_zeros_phyla_df_missingval)
```

    ## [1] 9 1

9 phyla have at least 1 missing value.

``` r
# histogram of number of zeros, starting at 1 zero
hist(counting_zeros_phyla_df_missingval$counting_zeros_phyla, 
     breaks = 60,
     main = "Histogram of Genera with Zero Relative Intensity",
     sub = "Starting at 1 Zero",
     xlab = "Number of zero relative intensity values",
     ylab = "Frequency")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-124-1.png)<!-- -->

``` r
# create table of number of phyla with more than 1 missing value
counting_zeros_phyla_df_missingval
```

    ##                                   counting_zeros_phyla
    ## Blastocladiomycota                                  55
    ## Chromerida                                          48
    ## Chytridiomycota                                     42
    ## Euglenida                                           31
    ## Glomeromycota                                       59
    ## Nanoarchaeota                                        8
    ## Phaeophyceae                                        35
    ## Xanthophyceae                                       29
    ## unclassified (derived from Fungi)                   46

### Filter for &lt;33% missingness

This would mean 33% missing values in our dataset.

``` r
# removing phyla that have 20 or more zeros
counting_zeros_phyla_df_missing20ormore <- counting_zeros_phyla_df %>%
  rownames_to_column(var = "rowname") %>%
  filter(counting_zeros_phyla >= 20) %>%
  column_to_rownames(var = "rowname")

# how many phyla have 20 or more missing value?
dim(counting_zeros_phyla_df_missing20ormore)
```

    ## [1] 8 1

8 phyla have more than 20 missing values.

``` r
# make a character vector from the rownames of previous data frame containing the phyla we want to get rid of
phyla.20zeros <- c(rownames(counting_zeros_phyla_df_missing20ormore))

# use select function to select all columns EXCEPT the ones in the character vector, we want to remove those
# and add in metadata
RelAbund.Phyla.Filt.zerofilt <- RelAbund.Phyla.Filt %>%
  select(everything(), -all_of(phyla.20zeros))

# check dimensions to make sure it filtered correctly
dim(RelAbund.Phyla.Filt.zerofilt)
```

    ## [1] 60 50

``` r
# removed 8, like we expected
```

Our final dataset has 45 phyla (because 5 columns are metadata).

Write final dataset genus rel abund to csv

``` r
write_csv(RelAbund.Phyla.Filt.zerofilt,
          file = "Phyla_RelAbund_Final_Filtered_WithMetadata.csv")
```

## Microbiome profile

See “Genera” section above for rarefaction curves and kronas plots

### Wrangling

Wrangling to enable collection of some summary statistics about our
microbiome profile.

Grab names of final phyla

``` r
# contains inplausible genera removed, but not removed for zeroes
dim(Phyla.Counts.Filt)
```

    ## [1] 53 62

``` r
Phyla.Counts.Filt[1:10, 1:5]
```

    ## # A tibble: 10 × 5
    ##    domain    phylum           `ShotgunWGS-Co…` `ShotgunWGS-Co…` `ShotgunWGS-Co…`
    ##    <chr>     <chr>                       <dbl>            <dbl>            <dbl>
    ##  1 Bacteria  Acidobacteria                2874             3717             2663
    ##  2 Bacteria  Actinobacteria             186789           277130           126155
    ##  3 Eukaryota Apicomplexa                   231              384              190
    ##  4 Bacteria  Aquificae                    1953             2254             1642
    ##  5 Eukaryota Ascomycota                   1491             2196             1281
    ##  6 Eukaryota Bacillariophyta               105              184              129
    ##  7 Bacteria  Bacteroidetes             1424565          1391417          1260217
    ##  8 Eukaryota Basidiomycota                 240              405              198
    ##  9 Eukaryota Blastocladiomyc…                0                0                0
    ## 10 Bacteria  Candidatus Pori…               26               49               21

``` r
# final filtered data
RelAbund.Phyla.Filt.zerofilt[1:10, 1:5]
```

    ## # A tibble: 10 × 5
    ##    Sample_Name                           Pig   Diet  Time_Point Diet_By_Time_Po…
    ##    <chr>                                 <fct> <fct> <fct>      <fct>           
    ##  1 ShotgunWGS-ControlPig6GutMicrobiome-… 6     Cont… Day 14     Control Day 14  
    ##  2 ShotgunWGS-ControlPig8GutMicrobiome-… 8     Cont… Day 0      Control Day 0   
    ##  3 ShotgunWGS-ControlPig3GutMicrobiome-… 3     Cont… Day 14     Control Day 14  
    ##  4 ShotgunWGS-TomatoPig14GutMicrobiome-… 14    Toma… Day 7      Tomato Day 7    
    ##  5 ShotgunWGS-ControlPig5GutMicrobiome-… 5     Cont… Day 7      Control Day 7   
    ##  6 ShotgunWGS-TomatoPig18GutMicrobiome-… 18    Toma… Day 7      Tomato Day 7    
    ##  7 ShotgunWGS-TomatoPig16GutMicrobiome-… 16    Toma… Day 7      Tomato Day 7    
    ##  8 ShotgunWGS-ControlPig10GutMicrobiome… 10    Cont… Day 7      Control Day 7   
    ##  9 ShotgunWGS-ControlPig2GutMicrobiome-… 2     Cont… Day 0      Control Day 0   
    ## 10 ShotgunWGS-TomatoPig18GutMicrobiome-… 18    Toma… Day 0      Tomato Day 0

``` r
dim(RelAbund.Phyla.Filt.zerofilt)
```

    ## [1] 60 50

``` r
# grab colnames which have all the final phyla
final_phyla <- colnames(RelAbund.Phyla.Filt.zerofilt)

final_phyla
```

    ##  [1] "Sample_Name"                                
    ##  [2] "Pig"                                        
    ##  [3] "Diet"                                       
    ##  [4] "Time_Point"                                 
    ##  [5] "Diet_By_Time_Point"                         
    ##  [6] "Acidobacteria"                              
    ##  [7] "Actinobacteria"                             
    ##  [8] "Apicomplexa"                                
    ##  [9] "Aquificae"                                  
    ## [10] "Ascomycota"                                 
    ## [11] "Bacillariophyta"                            
    ## [12] "Bacteroidetes"                              
    ## [13] "Basidiomycota"                              
    ## [14] "Candidatus Poribacteria"                    
    ## [15] "Chlamydiae"                                 
    ## [16] "Chlorobi"                                   
    ## [17] "Chloroflexi"                                
    ## [18] "Chlorophyta"                                
    ## [19] "Chrysiogenetes"                             
    ## [20] "Crenarchaeota"                              
    ## [21] "Cyanobacteria"                              
    ## [22] "Deferribacteres"                            
    ## [23] "Deinococcus-Thermus"                        
    ## [24] "Dictyoglomi"                                
    ## [25] "Elusimicrobia"                              
    ## [26] "Euryarchaeota"                              
    ## [27] "Fibrobacteres"                              
    ## [28] "Firmicutes"                                 
    ## [29] "Fusobacteria"                               
    ## [30] "Gemmatimonadetes"                           
    ## [31] "Hemichordata"                               
    ## [32] "Korarchaeota"                               
    ## [33] "Lentisphaerae"                              
    ## [34] "Microsporidia"                              
    ## [35] "Nanoarchaeota"                              
    ## [36] "Nematoda"                                   
    ## [37] "Nitrospirae"                                
    ## [38] "Placozoa"                                   
    ## [39] "Planctomycetes"                             
    ## [40] "Proteobacteria"                             
    ## [41] "Spirochaetes"                               
    ## [42] "Synergistetes"                              
    ## [43] "Tenericutes"                                
    ## [44] "Thaumarchaeota"                             
    ## [45] "Thermotogae"                                
    ## [46] "Verrucomicrobia"                            
    ## [47] "unclassified (derived from Bacteria)"       
    ## [48] "unclassified (derived from Eukaryota)"      
    ## [49] "unclassified (derived from Viruses)"        
    ## [50] "unclassified (derived from other sequences)"

``` r
# remove metadata colnames
final_phyla <- final_phyla[6:50]  

final_phyla <- as.data.frame(final_phyla)

final_phyla <- final_phyla %>%
  rename(phylum = final_phyla)
```

Get back domain and `inner_join` with `final_phyla` list

``` r
# pull from full dataset the domain and genus columns
Phyla.Counts.Filt.Domain.Phyla <- Phyla.Counts.Filt %>%
  select(domain, phylum)

head(Phyla.Counts.Filt.Domain.Phyla)
```

    ## # A tibble: 6 × 2
    ##   domain    phylum         
    ##   <chr>     <chr>          
    ## 1 Bacteria  Acidobacteria  
    ## 2 Bacteria  Actinobacteria 
    ## 3 Eukaryota Apicomplexa    
    ## 4 Bacteria  Aquificae      
    ## 5 Eukaryota Ascomycota     
    ## 6 Eukaryota Bacillariophyta

``` r
# want to join Genus.Counts.Filt.Domain.Genera with final_phyla
final_phyla_withdomain <- inner_join(final_phyla, Phyla.Counts.Filt.Domain.Phyla,
                                     by = "phylum")
```

### Count phyla

``` r
final_phyla_withdomain %>%
  count()
```

    ##    n
    ## 1 45

``` r
final_phyla_withdomain %>%
  group_by(domain) %>%
  count()
```

    ## # A tibble: 5 × 2
    ## # Groups:   domain [5]
    ##   domain              n
    ##   <chr>           <int>
    ## 1 Archaea             5
    ## 2 Bacteria           28
    ## 3 Eukaryota          10
    ## 4 other sequences     1
    ## 5 Viruses             1

### Most prevalent phyla

``` r
RelAbund.Phyla.Filt.zerofilt[1:5, 1:10]
```

    ## # A tibble: 5 × 10
    ##   Sample_Name              Pig   Diet  Time_Point Diet_By_Time_Po… Acidobacteria
    ##   <chr>                    <fct> <fct> <fct>      <fct>                    <dbl>
    ## 1 ShotgunWGS-ControlPig6G… 6     Cont… Day 14     Control Day 14        0.000741
    ## 2 ShotgunWGS-ControlPig8G… 8     Cont… Day 0      Control Day 0         0.000885
    ## 3 ShotgunWGS-ControlPig3G… 3     Cont… Day 14     Control Day 14        0.000690
    ## 4 ShotgunWGS-TomatoPig14G… 14    Toma… Day 7      Tomato Day 7          0.000732
    ## 5 ShotgunWGS-ControlPig5G… 5     Cont… Day 7      Control Day 7         0.000656
    ## # … with 4 more variables: Actinobacteria <dbl>, Apicomplexa <dbl>,
    ## #   Aquificae <dbl>, Ascomycota <dbl>

``` r
phyla_means <- RelAbund.Phyla.Filt.zerofilt %>%
  summarize_if(is.numeric, mean)

phyla_means_t <- t(phyla_means)
phyla_means_t <- as.data.frame(phyla_means_t)

phyla_means_t %>%
  rename(rel_abund_phyla = V1) %>%
  arrange(-rel_abund_phyla)
```

    ##                                             rel_abund_phyla
    ## Firmicutes                                     5.273546e-01
    ## Bacteroidetes                                  3.544950e-01
    ## Actinobacteria                                 4.660132e-02
    ## Proteobacteria                                 3.859541e-02
    ## Fusobacteria                                   4.300028e-03
    ## Euryarchaeota                                  4.280344e-03
    ## Spirochaetes                                   3.786741e-03
    ## Cyanobacteria                                  2.973060e-03
    ## Chloroflexi                                    2.143931e-03
    ## Fibrobacteres                                  1.647436e-03
    ## Thermotogae                                    1.475999e-03
    ## Synergistetes                                  1.429229e-03
    ## Chlorobi                                       1.371726e-03
    ## Verrucomicrobia                                1.113246e-03
    ## Tenericutes                                    8.325538e-04
    ## unclassified (derived from Viruses)            7.714176e-04
    ## Acidobacteria                                  7.491839e-04
    ## Deinococcus-Thermus                            7.206358e-04
    ## unclassified (derived from Eukaryota)          5.727411e-04
    ## unclassified (derived from Bacteria)           5.606341e-04
    ## Ascomycota                                     5.224676e-04
    ## Planctomycetes                                 5.021699e-04
    ## Aquificae                                      4.911229e-04
    ## Deferribacteres                                4.032208e-04
    ## Lentisphaerae                                  3.300053e-04
    ## Chlamydiae                                     3.262732e-04
    ## Dictyoglomi                                    3.130261e-04
    ## Elusimicrobia                                  2.375048e-04
    ## Crenarchaeota                                  2.006612e-04
    ## Nitrospirae                                    1.656603e-04
    ## Chlorophyta                                    1.318528e-04
    ## Apicomplexa                                    1.227757e-04
    ## Basidiomycota                                  8.845319e-05
    ## Nematoda                                       8.275259e-05
    ## Chrysiogenetes                                 8.110410e-05
    ## Gemmatimonadetes                               6.650590e-05
    ## Bacillariophyta                                3.943884e-05
    ## unclassified (derived from other sequences)    2.783963e-05
    ## Microsporidia                                  2.109075e-05
    ## Thaumarchaeota                                 1.895141e-05
    ## Korarchaeota                                   1.880226e-05
    ## Placozoa                                       1.426135e-05
    ## Candidatus Poribacteria                        1.034724e-05
    ## Hemichordata                                   4.906748e-06
    ## Nanoarchaeota                                  1.627075e-06

The most prevalent phyla are Firmicutes (52.7% average abundance),
Bacteroidetes (35.4%), Actinobacteria (4.7%), Proteobacteria (3.9%) and
Fusobaceria (0.43%).

What percent of the reads are from Bacteria?

``` r
final_phyla_bacteriaonly <- final_phyla_withdomain %>%
  filter(domain == "Bacteria")

final_phyla_bacteriaonly <- final_phyla_bacteriaonly$phylum

# select columns corresponding to bacteria
RelAbund.Phyla.Filt.zerofilt.baconly <- RelAbund.Phyla.Filt.zerofilt %>%
  select(contains(final_phyla_bacteriaonly)) 

# create rowsums
RelAbund.Phyla.Filt.zerofilt.baconly <- RelAbund.Phyla.Filt.zerofilt.baconly %>%
  mutate(rowsums = rowSums(RelAbund.Phyla.Filt.zerofilt.baconly[])) 

mean(RelAbund.Phyla.Filt.zerofilt.baconly$rowsums)
```

    ## [1] 0.9930777

``` r
sd(RelAbund.Phyla.Filt.zerofilt.baconly$rowsums)
```

    ## [1] 0.002045929

## PERMANOVA

### All samples, full model

Repeated measures, using Pig as a block and set permutations using
`how()` ORIGINAL BLOCK

``` r
set.seed(2021)
# create factors
factors_time_diet_pig_phyla <- RelAbund.Phyla.Filt.zerofilt %>% 
  select(Time_Point, Diet, Pig)

# create permutations
perm_time_diet_pig_phyla <- how(nperm = 9999)
setBlocks(perm_time_diet_pig_phyla) <- with(factors_time_diet_pig_phyla, Pig)

# run permanova
AllData.Phyla.Filt.permanova <- adonis2(RelAbund.Phyla.Filt.zerofilt[,-c(1:5)]~Diet*Time_Point,
                                        data = factors_time_diet_pig_phyla,
                                        permutations = perm_time_diet_pig_phyla,
                                        method = "bray")

AllData.Phyla.Filt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  with(factors_time_diet_pig_phyla, Pig) 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = RelAbund.Phyla.Filt.zerofilt[, -c(1:5)] ~ Diet * Time_Point, data = factors_time_diet_pig_phyla, permutations = perm_time_diet_pig_phyla, method = "bray")
    ##                 Df SumOfSqs      R2      F Pr(>F)   
    ## Diet             1 0.007675 0.02656 1.6746 0.0136 * 
    ## Time_Point       2 0.028496 0.09860 3.1087 0.0046 **
    ## Diet:Time_Point  2 0.005338 0.01847 0.5824 0.4915   
    ## Residual        54 0.247497 0.85637                 
    ## Total           59 0.289007 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

-   Diet: p = 0.0150, significant  
-   Time\_Point: p = 0.0054, significant  
-   Diet\*Time\_Point: p = 0.4870, non-significant

Interaction:

``` r
# create factors
Pig <- as.factor(RelAbund.Phyla.Filt.zerofilt$Pig)
Diet <- as.factor(RelAbund.Phyla.Filt.zerofilt$Diet)


# create permutations
perm_time_diet_pig_phyla <- how(within = Within(type="series", constant=TRUE),
                                plots = Plots(strata=Pig, type="free",))
# run permanova
AllData.Phyla.Filt.permanova <- adonis2(RelAbund.Phyla.Filt.zerofilt[,-c(1:5)]~Diet*Time_Point,
                                        data = factors_time_diet_pig_phyla,
                                        permutations = perm_time_diet_pig_phyla,
                                        method = "bray",
                                        by = "margin")

AllData.Phyla.Filt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Plots: Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = RelAbund.Phyla.Filt.zerofilt[, -c(1:5)] ~ Diet * Time_Point, data = factors_time_diet_pig_phyla, permutations = perm_time_diet_pig_phyla, method = "bray", by = "margin")
    ##                 Df SumOfSqs      R2      F Pr(>F)
    ## Diet:Time_Point  2 0.005338 0.01847 0.5824  0.515
    ## Residual        54 0.247497 0.85637              
    ## Total           59 0.289007 1.00000

Interaction not significant (p=.51), so remove from model

``` r
# create factors
Pig <- as.factor(RelAbund.Phyla.Filt.zerofilt$Pig)
Diet <- as.factor(RelAbund.Phyla.Filt.zerofilt$Diet)


# create permutations
perm_time_diet_pig_phyla <- how(within = Within(type="series", constant=TRUE),
                                plots = Plots(strata=Pig, type = "free"))
# run permanova
AllData.Phyla.Filt.permanova <- adonis2(RelAbund.Phyla.Filt.zerofilt[,-c(1:5)]~Diet + Time_Point,
                                        data = factors_time_diet_pig_phyla,
                                        permutations = perm_time_diet_pig_phyla,
                                        method = "bray",
                                        by = "margin")

AllData.Phyla.Filt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Plots: Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = RelAbund.Phyla.Filt.zerofilt[, -c(1:5)] ~ Diet + Time_Point, data = factors_time_diet_pig_phyla, permutations = perm_time_diet_pig_phyla, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)  
    ## Diet        1 0.007675 0.02656 1.7000  0.295  
    ## Time_Point  2 0.028496 0.09860 3.1558  0.015 *
    ## Residual   56 0.252835 0.87484                
    ## Total      59 0.289007 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Test for homogeneity of multivariate dispersions

``` r
dis <- vegdist(RelAbund.Phyla.Filt.zerofilt[,-c(1:5)], method = "bray")
mod <- betadisper(dis, Diet)
permutest(mod)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq    Mean Sq     F N.Perm Pr(>F)
    ## Groups     1 0.000153 0.00015277 0.121    999  0.726
    ## Residuals 58 0.073213 0.00126230

Non significant! good for our PERMANOVA test validity

### Post Hoc PERMANOVA within Time

#### Within Control Diet Only

Effect of control diet over time.

``` r
# filter data set for only control samples
control.RelAbund.Phyla.zerofilt <- subset(RelAbund.Phyla.Filt.zerofilt, Diet == "Control")

# create factors
factors_control_pig_phyla <- droplevels(control.RelAbund.Phyla.zerofilt %>% 
  select(Time_Point, Pig))

# create permutations
perm_control_pig_phyla <- how(within = Within(type="series", constant=TRUE),
                                plots = Plots(strata=factors_control_pig_phyla$Pig, type = "free"))

# run PERMANOVA
Control.ByTime.Phyla.zerofilt.permanova <- adonis2(control.RelAbund.Phyla.zerofilt[,-c(1:5)]~Time_Point,
        data = factors_control_pig_phyla,
        permutations = perm_control_pig_phyla, 
        method = "bray",
        by = "margin")

Control.ByTime.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_pig_phyla$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = control.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_pig_phyla, permutations = perm_control_pig_phyla, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## Time_Point  2 0.025943 0.17486 2.8609   0.01 **
    ## Residual   27 0.122422 0.82514                 
    ## Total      29 0.148365 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Significant effect of time (p = 0.005) within control samples. Beta
diversity changing with time. Now the question is where is the
difference coming from (ie. between which time points?)

##### Control T1 vs Control T2

``` r
# filter data set for only samples at T1 and T2
control.T1T2.RelAbund.Phyla.zerofilt <- subset(control.RelAbund.Phyla.zerofilt, Time_Point != "Day 14")

# create factors
factors_control_T1T2_pig_phyla <- droplevels(control.T1T2.RelAbund.Phyla.zerofilt %>% 
  select(Time_Point, Pig))

# create permutations
perm_control_T1T2_pig_phyla <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_control_T1T2_pig_phyla$Pig,
                                                 type = "free"))

# run PERMANOVA
Control.T1T2.Phyla.zerofilt.permanova <- adonis2(control.T1T2.RelAbund.Phyla.zerofilt[,-c(1:5)]~Time_Point,
        data = factors_control_T1T2_pig_phyla,
        permutations = perm_control_T1T2_pig_phyla, 
        method = "bray",
        by = "margin")

Control.T1T2.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_T1T2_pig_phyla$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = control.T1T2.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_T1T2_pig_phyla, permutations = perm_control_T1T2_pig_phyla, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)  
    ## Time_Point  1  0.01404 0.11986 2.4513   0.03 *
    ## Residual   18  0.10309 0.88014                
    ## Total      19  0.11713 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

p=.085 so not significant between T1 and T2

##### Control T1 vs Control T3

``` r
# filter data set for only samples at T1 and T3
control.T1T3.RelAbund.Phyla.zerofilt <- subset(control.RelAbund.Phyla.zerofilt, Time_Point != "Day 7")

# create factors
factors_control_T1T3_pig_phyla <- droplevels(control.T1T3.RelAbund.Phyla.zerofilt %>% 
  select(Time_Point, Pig))

# create permutations
perm_control_T1T3_pig_phyla <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_control_T1T3_pig_phyla$Pig,
                                                 type = "free"))

# run PERMANOVA
Control.T1T3.Phyla.zerofilt.permanova <- adonis2(control.T1T3.RelAbund.Phyla.zerofilt[,-c(1:5)]~Time_Point,
        data = factors_control_T1T3_pig_phyla,
        permutations = perm_control_T1T3_pig_phyla, 
        method = "bray",
        by = "margin")

Control.T1T3.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_T1T3_pig_phyla$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = control.T1T3.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_T1T3_pig_phyla, permutations = perm_control_T1T3_pig_phyla, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)  
    ## Time_Point  1 0.022462 0.28421 7.1469   0.02 *
    ## Residual   18 0.056572 0.71579                
    ## Total      19 0.079034 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

P = .02 so significant. There is a significant difference between T1 and
T3 in the control diet pigs

##### Control T2 vs Control T3

``` r
# filter data set for only samples at T2 and T3
control.T2T3.RelAbund.Phyla.zerofilt <- subset(control.RelAbund.Phyla.zerofilt, Time_Point != "Day 0")

# create factors
factors_control_T2T3_pig_phyla <- droplevels(control.T2T3.RelAbund.Phyla.zerofilt %>% 
  select(Time_Point, Pig))

# create permutations
perm_control_T2T3_pig_phyla <- how(within = Within(type="series", constant=TRUE),
                                   plots = Plots(strata=factors_control_T2T3_pig_phyla$Pig,
                                                 type = "free"))

# run PERMANOVA
Control.T2T3.Phyla.zerofilt.permanova <- adonis2(control.T2T3.RelAbund.Phyla.zerofilt[,-c(1:5)]~Time_Point,
        data = factors_control_T2T3_pig_phyla,
        permutations = perm_control_T2T3_pig_phyla, 
        method = "bray",
        by = "margin")

Control.T2T3.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_control_T2T3_pig_phyla$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = control.T2T3.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_control_T2T3_pig_phyla, permutations = perm_control_T2T3_pig_phyla, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2    F Pr(>F)
    ## Time_Point  1 0.002413 0.02755 0.51   0.33
    ## Residual   18 0.085178 0.97245            
    ## Total      19 0.087591 1.00000

P = .315 so not significant

#### Within Tomato Diet Only

Effect of tomato diet over time.

``` r
# filter data for only tomato samples
tomato.RelAbund.Phyla.zerofilt <- subset(RelAbund.Phyla.Filt.zerofilt, Diet == "Tomato")

# create factors
factors_tomato_pig_phyla <- tomato.RelAbund.Phyla.zerofilt %>% 
  select(Time_Point, Pig)

# create permutations
perm_tomato_pig_phyla <- how(within = Within(type="series", constant=TRUE),
                             plots = Plots(strata=factors_tomato_pig_phyla$Pig, type = "free"))

# run PERMANOVA
tomato.ByTime.Phyla.zerofilt.permanova <- adonis2(tomato.RelAbund.Phyla.zerofilt[,-c(1:5)]~Time_Point,
                                                  data = factors_tomato_pig_phyla,
                                                  permutations = perm_tomato_pig_phyla, 
                                                  method = "bray",
                                                  by = "margin")

tomato.ByTime.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under NA model
    ## Marginal effects of terms
    ## Plots: factors_tomato_pig_phyla$Pig, plot permutation: free
    ## Permutation: series constant permutation within each Plot
    ## Number of permutations: 199
    ## 
    ## adonis2(formula = tomato.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Time_Point, data = factors_tomato_pig_phyla, permutations = perm_tomato_pig_phyla, method = "bray", by = "margin")
    ##            Df SumOfSqs      R2      F Pr(>F)
    ## Time_Point  2 0.007891 0.05935 0.8517   0.34
    ## Residual   27 0.125075 0.94065              
    ## Total      29 0.132966 1.00000

Non-significant effect of time (p = 0.325) within tomato samples. So no
post hoc tests necessary.

### Subset by time

#### Day 0 Only

Effect of diet at day 0.

``` r
# filter data set for only day 0 samples
d0.RelAbund.Phyla.zerofilt <- subset(RelAbund.Phyla.Filt.zerofilt, Time_Point == "Day 0")

# create factors
# don't need to include pig, since no repeated measures here 
# only testing Diet within a time point
factors_day0_phyla <- d0.RelAbund.Phyla.zerofilt %>% 
  select(Diet)

# create permutations
perm_day0_phyla <- how(nperm = 9999)

# run PERMANOVA
d0.Phyla.zerofilt.permanova <- adonis2(d0.RelAbund.Phyla.zerofilt[,-c(1:5)]~Diet,
                                       data = factors_day0_phyla,
                                       permutations = perm_day0_phyla, 
                                       method = "bray")

d0.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = d0.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Diet, data = factors_day0_phyla, permutations = perm_day0_phyla, method = "bray")
    ##          Df SumOfSqs      R2     F Pr(>F)
    ## Diet      1 0.003515 0.04898 0.927 0.3692
    ## Residual 18 0.068249 0.95102             
    ## Total    19 0.071764 1.00000

Non-significant effect of diet (p=0.376) at day 0.

#### Day 7 Only

Effect of diet at day 7.

``` r
# filter data set for only day 7 samples
d7.RelAbund.Phyla.zerofilt <- subset(RelAbund.Phyla.Filt.zerofilt, Time_Point == "Day 7")

# create factors
# don't need to include pig, since no repeated measures here 
# only testing Diet within a time point
factors_day7_phyla <- d7.RelAbund.Phyla.zerofilt %>% 
  select(Diet)

# create permutations
perm_day7_phyla <- how(nperm = 9999)

# run PERMANOVA
d7.Phyla.zerofilt.permanova <- adonis2(d7.RelAbund.Phyla.zerofilt[,-c(1:5)]~Diet,
                                       data = factors_day7_phyla,
                                       permutations = perm_day7_phyla, 
                                       method = "bray")

d7.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = d7.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Diet, data = factors_day7_phyla, permutations = perm_day7_phyla, method = "bray")
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Diet      1 0.005267 0.04205 0.7901 0.4009
    ## Residual 18 0.119990 0.95795              
    ## Total    19 0.125257 1.00000

Non-significant effect of diet (p=0.4097) at day 7.

#### Day 14 Only

Effect of diet at day 14.

``` r
# filter data set for only day 14 samples
d14.RelAbund.Phyla.zerofilt <- subset(RelAbund.Phyla.Filt.zerofilt, Time_Point == "Day 14")

# create factors
# don't need to include pig, since no repeated measures here 
# only testing Diet within a time point
factors_day14_phyla <- d14.RelAbund.Phyla.zerofilt %>% 
  select(Diet)

# create permutations
perm_day14_phyla <- how(nperm = 9999)

# run PERMANOVA
d14.Phyla.zerofilt.permanova <- adonis2(d14.RelAbund.Phyla.zerofilt[,-c(1:5)]~Diet,
                                       data = factors_day14_phyla,
                                       permutations = perm_day14_phyla, 
                                       method = "bray")

d14.Phyla.zerofilt.permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = d14.RelAbund.Phyla.zerofilt[, -c(1:5)] ~ Diet, data = factors_day14_phyla, permutations = perm_day14_phyla, method = "bray")
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Diet      1 0.004232 0.06665 1.2854 0.2691
    ## Residual 18 0.059258 0.93335              
    ## Total    19 0.063490 1.00000

Non-significant effect of diet (p=0.256) at day 14.

## PCoA Beta Diversity

### All samples

``` r
# calculate distances
phyla.filt.dist.zeros <- vegdist(RelAbund.Phyla.Filt.zerofilt[6:ncol(RelAbund.Phyla.Filt.zerofilt)], 
                                 method = "bray")

# do multi-dimensional scaling (the PCoA calculations) on those distances
scale.phyla.filt.zerofilt <- cmdscale(phyla.filt.dist.zeros, k=2)

# make into data frame and bind metadata
scale.phyla.filt.zerofilt.df <- as.data.frame(cbind(scale.phyla.filt.zerofilt, AllSamples.Metadata))

# do PCoA again, but get eigen values
scale.phyla.filt.zerofilt.eig <- cmdscale(phyla.filt.dist.zeros, k=2, eig = TRUE)

# convert eigenvalues to percentages and assign to a variable
eigs.phyla.filt.zerofilt <- (100*((scale.phyla.filt.zerofilt.eig$eig)/(sum(scale.phyla.filt.zerofilt.eig$eig))))

# round the converted eigenvalues
round.eigs.phyla.filt.zerofilt <- round(eigs.phyla.filt.zerofilt, 3)
```

Plot

``` r
PCoA_phyla_20zeros_allsamples <- scale.phyla.filt.zerofilt.df %>%
ggplot(aes(x=`1`, y=`2`, fill = Diet_By_Time_Point)) +
  geom_point(size = 3, color = "black", shape = 21, alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4", "sienna1","firebrick3","tomato4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))+
  labs(x=paste("PC1: ", round.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Diet & Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level") 

PCoA_phyla_20zeros_allsamples
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-147-1.png)<!-- -->

``` r
ggsave("Figures/BetaDiversity_PCoA_Phyla_allsamples.png", 
       plot = PCoA_phyla_20zeros_allsamples, 
       dpi = 800, 
       width = 10, 
       height = 8)
```

#### Facet by time point

Re-level factors

``` r
scale.phyla.filt.zerofilt.df <- scale.phyla.filt.zerofilt.df %>% 
  mutate(Time_Point = fct_relevel(Time_Point, c("Day 0", "Day 7", "Day 14")))
```

``` r
PCoA_phyla_20zeros_facetbytime <- scale.phyla.filt.zerofilt.df %>%
ggplot(aes(x=`1`, y=`2`, fill = Diet_By_Time_Point)) +
  geom_hline(yintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_point(size = 3, color = "black", shape = 21, alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4", "sienna1","firebrick3","tomato4")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x=paste("PC1: ", round.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Diet & Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level, Subset by Time Point") +
  facet_wrap(~Time_Point)

PCoA_phyla_20zeros_facetbytime
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-150-1.png)<!-- -->

``` r
ggsave("Figures/BetaDiversity_PCoA_Phyla_FacetByTimePoint.png", 
       plot = PCoA_phyla_20zeros_facetbytime, 
       dpi = 800, 
       width = 10, 
       height = 6)
```

#### Facet by diet

``` r
PCoA_phyla_20zeros_facetbydiet <- scale.phyla.filt.zerofilt.df %>%
ggplot(aes(x=`1`, y=`2`, fill = Diet_By_Time_Point)) +
  geom_hline(yintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = 0, color = "light grey", linetype = "dashed", size = 0.3) +
  geom_point(size = 3, color = "black", shape = 21, alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4", "sienna1","firebrick3","tomato4")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x=paste("PC1: ", round.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Diet & Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level, Subset by Diet") +
  facet_wrap(~Diet)

PCoA_phyla_20zeros_facetbydiet
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-152-1.png)<!-- -->

``` r
ggsave("Figures/BetaDiversity_PCoA_Phyla_FacetByDiet.png", 
       plot = PCoA_phyla_20zeros_facetbydiet, 
       dpi = 800, 
       width = 10, 
       height = 8)
```

### Subset

Ended up not using this as part of the paper. Since the input is
different here (i.e., the PCoA only has the subset data as an input) the
output looks slightly different.

#### Control only

``` r
# calculate distances
control.phyla.filt.dist.zeros <- vegdist(control.RelAbund.Phyla.zerofilt[,-c(1:5)], method = "bray")

# do PCoA calculations
control.scale.phyla.filt.zerofilt <- cmdscale(control.phyla.filt.dist.zeros, k=2)

# filter meta data
meta.control <- subset(AllSamples.Metadata, Diet == "Control")

# make pcoa table into data frame and bind metadata to it
control.scale.phyla.filt.zerofilt.df <- as.data.frame(cbind(meta.control, control.scale.phyla.filt.zerofilt))

# do PCoA again, but get eigenvalues
control.scale.phyla.filt.zerofilt.eig <- cmdscale(control.phyla.filt.dist.zeros, k=2, eig = TRUE)

# convert eigenvalues to percentages and assign to a variable
control.eigs.phyla.filt.zerofilt <- (100*((control.scale.phyla.filt.zerofilt.eig$eig)/sum(control.scale.phyla.filt.zerofilt.eig$eig)))

# round the eigenvalues
round.control.eigs.phyla.filt.zerofilt <- round(control.eigs.phyla.filt.zerofilt, 3)
```

Re-level factors

``` r
control.scale.phyla.filt.zerofilt.df$Time_Point <- factor(control.scale.phyla.filt.zerofilt.df$Time_Point, 
                                                          levels = c("Day 0", "Day 7", "Day 14"))
```

Plot

``` r
control.scale.phyla.filt.zerofilt.df %>%
ggplot(aes(x=`1`, y=`2`, fill = Time_Point)) +
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=paste("PC1: ", round.control.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.control.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level, Control Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-156-1.png)<!-- -->

#### Tomato only

``` r
# calculate distances
tomato.phyla.filt.dist.zeros <- vegdist(tomato.RelAbund.Phyla.zerofilt[,-c(1:5)], method = "bray")

# do PCoA calculations
tomato.scale.phyla.filt.zerofilt <- cmdscale(tomato.phyla.filt.dist.zeros, k=2)

# filter meta data
meta.tomato <- subset(AllSamples.Metadata, Diet == "Tomato")

# make pcoa table into data frame and bind metadata to it
tomato.scale.phyla.filt.zerofilt.df <- as.data.frame(cbind(meta.tomato, tomato.scale.phyla.filt.zerofilt))

# do PCoA again, but get eigenvalues
tomato.scale.phyla.filt.zerofilt.eig <- cmdscale(tomato.phyla.filt.dist.zeros, k=2, eig = TRUE)

# convert eigenvalues to percentages and assign to a variable
tomato.eigs.phyla.filt.zerofilt <- (100*((tomato.scale.phyla.filt.zerofilt.eig$eig)/sum(tomato.scale.phyla.filt.zerofilt.eig$eig)))

# round the eigenvalues
round.tomato.eigs.phyla.filt.zerofilt <- round(tomato.eigs.phyla.filt.zerofilt, 3)
```

Re-level factors

``` r
tomato.scale.phyla.filt.zerofilt.df$Time_Point <- factor(tomato.scale.phyla.filt.zerofilt.df$Time_Point, 
                                                         levels = c("Day 0", "Day 7", "Day 14"))
```

Plot

``` r
tomato.scale.phyla.filt.zerofilt.df %>%
ggplot(aes(x=`1`, y=`2`, fill = Time_Point)) +
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values=c("sienna1","firebrick3","tomato4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=paste("PC1: ", round.tomato.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.tomato.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level, Tomato Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-159-1.png)<!-- -->

#### Day 0 Only

``` r
# calculate distances
d0.phyla.filt.dist.zeros <- vegdist(d0.RelAbund.Phyla.zerofilt[,-c(1:5)], method = "bray")

# do PCoA calculations
d0.scale.phyla.filt.zerofilt <- cmdscale(d0.phyla.filt.dist.zeros, k=2)

# filter meta data
meta.d0 <- subset(AllSamples.Metadata, Time_Point == "Day 0")

# make pcoa table into data frame and bind metadata to it
d0.scale.phyla.filt.zerofilt.df <- as.data.frame(cbind(meta.d0, d0.scale.phyla.filt.zerofilt))

# do PCoA again, but get eigenvalues
d0.scale.phyla.filt.zerofilt.eig <- cmdscale(d0.phyla.filt.dist.zeros, k=2, eig = TRUE)

# convert eigenvalues to percentages and assign to a variable
d0.eigs.phyla.filt.zerofilt <- (100*((d0.scale.phyla.filt.zerofilt.eig$eig)/sum(d0.scale.phyla.filt.zerofilt.eig$eig)))

# round the eigenvalues
round.d0.eigs.phyla.filt.zerofilt <- round(d0.eigs.phyla.filt.zerofilt, 3)
```

Plot

``` r
d0.scale.phyla.filt.zerofilt.df %>%
  ggplot(aes(x = `1`, y = `2`, fill = Diet)) +
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values=c("steelblue2", "tomato2")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=paste("PC1: ", round.d0.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.d0.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level, Day 0 Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-161-1.png)<!-- -->

#### Day 7 Only

``` r
# calculate distances
d7.phyla.filt.dist.zeros <- vegdist(d7.RelAbund.Phyla.zerofilt[,-c(1:5)], method = "bray")

# do PCoA calculations
d7.scale.phyla.filt.zerofilt <- cmdscale(d7.phyla.filt.dist.zeros, k=2)

# filter meta data
meta.d7 <- subset(AllSamples.Metadata, Time_Point == "Day 7")

# make pcoa table into data frame and bind metadata to it
d7.scale.phyla.filt.zerofilt.df <- as.data.frame(cbind(meta.d7, d7.scale.phyla.filt.zerofilt))

# do PCoA again, but get eigenvalues
d7.scale.phyla.filt.zerofilt.eig <- cmdscale(d7.phyla.filt.dist.zeros, k=2, eig = TRUE)

# convert eigenvalues to percentages and assign to a variable
d7.eigs.phyla.filt.zerofilt <- (100*((d7.scale.phyla.filt.zerofilt.eig$eig)/sum(d7.scale.phyla.filt.zerofilt.eig$eig)))

# round the eigenvalues
round.d7.eigs.phyla.filt.zerofilt <- round(d7.eigs.phyla.filt.zerofilt, 3)
```

Plot

``` r
d7.scale.phyla.filt.zerofilt.df %>%
  ggplot(aes(x = `1`, y = `2`, fill = Diet)) +
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values=c("steelblue2", "tomato2")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=paste("PC1: ", round.d7.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.d7.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level, Day 7 Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-163-1.png)<!-- -->

#### Day 14 Only

``` r
# calculate distances
d14.phyla.filt.dist.zeros <- vegdist(d14.RelAbund.Phyla.zerofilt[,-c(1:5)], method = "bray")
# do PCoA calculations
d14.scale.phyla.filt.zerofilt <- cmdscale(d14.phyla.filt.dist.zeros, k=2)

# filter meta data
meta.d14 <- subset(AllSamples.Metadata, Time_Point == "Day 14")

# make pcoa table into data frame and bind metadata to it
d14.scale.phyla.filt.zerofilt.df <- as.data.frame(cbind(meta.d14, d14.scale.phyla.filt.zerofilt))

# do PCoA again, but get eigenvalues
d14.scale.phyla.filt.zerofilt.eig <- cmdscale(d14.phyla.filt.dist.zeros, k=2, eig = TRUE)

# convert eigenvalues to percentages and assign to a variable
d14.eigs.phyla.filt.zerofilt <- (100*((d14.scale.phyla.filt.zerofilt.eig$eig)/sum(d14.scale.phyla.filt.zerofilt.eig$eig)))

# round the eigenvalues
round.d14.eigs.phyla.filt.zerofilt <- round(d14.eigs.phyla.filt.zerofilt, 3)
```

Plot

``` r
d14.scale.phyla.filt.zerofilt.df %>%
  ggplot(aes(x = `1`, y = `2`, fill = Diet)) +
  geom_point(size=3, shape = 21, color = "black", alpha = 0.9) +
  scale_fill_manual(values=c("steelblue2", "tomato2")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  labs(x=paste("PC1: ", round.d14.eigs.phyla.filt.zerofilt[1], "%"), 
       y=paste("PC2: ", round.d14.eigs.phyla.filt.zerofilt[2], "%"), 
       fill="Time Point",
       title = "Beta Diversity",
       subtitle = "Phyla Level, Day 14 Only")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-165-1.png)<!-- -->

## Bacteriodetes, Firmicutes, and B to F

Given a priori interest in the phyla Bacteriodetes and Firmicutes, we
are conducted repeated measures ANOVA analysis for their changes in our
samples. The ratio of Bacteriodetes to Firmicutes (B to F) is a commonly
used metric for assessing the health of the microbiome, with a higher B
to F being more beneficial.

### Wrangling

``` r
dim(RelAbund.Phyla.Filt.zerofilt)
```

    ## [1] 60 50

60 samples, and 45 phyla (5 columns are metadata).

Re-level Time\_Point

``` r
RelAbund.Phyla.Filt.zerofilt <- RelAbund.Phyla.Filt.zerofilt %>%
  mutate(Time_Point = fct_relevel(Time_Point, c("Day 0", "Day 7", "Day 14")))

levels(RelAbund.Phyla.Filt.zerofilt$Time_Point)
```

    ## [1] "Day 0"  "Day 7"  "Day 14"

Add column `Other_phyla` with the sum of all phyla that are not
Bacteroidetes or Firmicutes

``` r
RelAbund.Phyla.Filt.zerofilt.withother <- RelAbund.Phyla.Filt.zerofilt %>%
  mutate(Other_phyla = rowSums(select(.[6:ncol(.)], !contains(c("Bacteroidetes", "Firmicutes")))))

kable(head(RelAbund.Phyla.Filt.zerofilt.withother))
```

| Sample\_Name                              | Pig | Diet    | Time\_Point | Diet\_By\_Time\_Point | Acidobacteria | Actinobacteria | Apicomplexa | Aquificae | Ascomycota | Bacillariophyta | Bacteroidetes | Basidiomycota | Candidatus Poribacteria | Chlamydiae |  Chlorobi | Chloroflexi | Chlorophyta | Chrysiogenetes | Crenarchaeota | Cyanobacteria | Deferribacteres | Deinococcus-Thermus | Dictyoglomi | Elusimicrobia | Euryarchaeota | Fibrobacteres | Firmicutes | Fusobacteria | Gemmatimonadetes | Hemichordata | Korarchaeota | Lentisphaerae | Microsporidia | Nanoarchaeota |  Nematoda | Nitrospirae | Placozoa | Planctomycetes | Proteobacteria | Spirochaetes | Synergistetes | Tenericutes | Thaumarchaeota | Thermotogae | Verrucomicrobia | unclassified (derived from Bacteria) | unclassified (derived from Eukaryota) | unclassified (derived from Viruses) | unclassified (derived from other sequences) | Other\_phyla |
|:------------------------------------------|:----|:--------|:------------|:----------------------|--------------:|---------------:|------------:|----------:|-----------:|----------------:|--------------:|--------------:|------------------------:|-----------:|----------:|------------:|------------:|---------------:|--------------:|--------------:|----------------:|--------------------:|------------:|--------------:|--------------:|--------------:|-----------:|-------------:|-----------------:|-------------:|-------------:|--------------:|--------------:|--------------:|----------:|------------:|---------:|---------------:|---------------:|-------------:|--------------:|------------:|---------------:|------------:|----------------:|-------------------------------------:|--------------------------------------:|------------------------------------:|--------------------------------------------:|-------------:|
| ShotgunWGS-ControlPig6GutMicrobiome-Day14 | 6   | Control | Day 14      | Control Day 14        |     0.0007406 |      0.0481336 |   0.0000595 | 0.0005033 |  0.0003842 |        2.71e-05 |     0.3670956 |      6.18e-05 |                6.70e-06 |  0.0001422 | 0.0012598 |   0.0020208 |   0.0000953 |       8.53e-05 |     0.0001670 |     0.0022775 |       0.0003850 |           0.0006393 |   0.0003136 |     0.0001629 |     0.0033951 |     0.0012287 |  0.5308271 |    0.0039555 |        0.0000544 |      6.7e-06 |     1.93e-05 |     0.0001971 |      1.26e-05 |       1.0e-06 | 0.0000459 |   0.0001420 | 1.75e-05 |      0.0003925 |      0.0271370 |    0.0029683 |     0.0011475 |   0.0007123 |       1.44e-05 |   0.0012921 |       0.0008269 |                            0.0003085 |                             0.0003247 |                           0.0003984 |                                    1.24e-05 |    0.1020763 |
| ShotgunWGS-ControlPig8GutMicrobiome-Day0  | 8   | Control | Day 0       | Control Day 0         |     0.0008850 |      0.0659807 |   0.0000914 | 0.0005366 |  0.0005228 |        4.38e-05 |     0.3312767 |      9.64e-05 |                1.17e-05 |  0.0001974 | 0.0014459 |   0.0022885 |   0.0001359 |       9.90e-05 |     0.0002148 |     0.0032408 |       0.0004747 |           0.0008538 |   0.0003700 |     0.0002398 |     0.0045655 |     0.0016578 |  0.5293436 |    0.0050574 |        0.0000809 |      4.3e-06 |     2.31e-05 |     0.0004364 |      2.10e-05 |       1.2e-06 | 0.0000631 |   0.0001898 | 1.62e-05 |      0.0006266 |      0.0368314 |    0.0041577 |     0.0017830 |   0.0008883 |       1.90e-05 |   0.0016916 |       0.0014957 |                            0.0004095 |                             0.0009973 |                           0.0006252 |                                    7.90e-06 |    0.1393790 |
| ShotgunWGS-ControlPig3GutMicrobiome-Day14 | 3   | Control | Day 14      | Control Day 14        |     0.0006897 |      0.0326727 |   0.0000492 | 0.0004253 |  0.0003318 |        3.34e-05 |     0.3263817 |      5.13e-05 |                5.40e-06 |  0.0001435 | 0.0011574 |   0.0019673 |   0.0000938 |       7.02e-05 |     0.0001652 |     0.0029204 |       0.0003675 |           0.0006011 |   0.0002872 |     0.0001810 |     0.0033125 |     0.0012911 |  0.5870258 |    0.0037180 |        0.0000624 |      4.1e-06 |     1.35e-05 |     0.0002005 |      1.45e-05 |       1.0e-06 | 0.0000378 |   0.0001458 | 8.00e-06 |      0.0004066 |      0.0271624 |    0.0028287 |     0.0010743 |   0.0006457 |       1.53e-05 |   0.0013084 |       0.0009681 |                            0.0002973 |                             0.0003279 |                           0.0005188 |                                    1.61e-05 |    0.0865920 |
| ShotgunWGS-TomatoPig14GutMicrobiome-Day7  | 14  | Tomato  | Day 7       | Tomato Day 7          |     0.0007324 |      0.0329202 |   0.0001398 | 0.0005384 |  0.0005593 |        5.41e-05 |     0.3461498 |      9.49e-05 |                1.41e-05 |  0.0001856 | 0.0013024 |   0.0023710 |   0.0001531 |       9.49e-05 |     0.0002730 |     0.0028512 |       0.0004519 |           0.0008048 |   0.0003529 |     0.0002988 |     0.0069066 |     0.0014564 |  0.5231174 |    0.0046147 |        0.0000508 |      4.2e-06 |     3.00e-05 |     0.0004095 |      2.00e-05 |       5.0e-06 | 0.0001148 |   0.0002006 | 1.58e-05 |      0.0004494 |      0.0597393 |    0.0046887 |     0.0015563 |   0.0010752 |       3.08e-05 |   0.0016628 |       0.0010478 |                            0.0009662 |                             0.0005509 |                           0.0008405 |                                    9.57e-05 |    0.1307244 |
| ShotgunWGS-ControlPig5GutMicrobiome-Day7  | 5   | Control | Day 7       | Control Day 7         |     0.0006564 |      0.0463444 |   0.0000557 | 0.0004617 |  0.0003835 |        3.61e-05 |     0.2599966 |      5.93e-05 |                7.20e-06 |  0.0001644 | 0.0012486 |   0.0020345 |   0.0001263 |       8.89e-05 |     0.0002123 |     0.0034222 |       0.0003796 |           0.0006821 |   0.0002920 |     0.0002165 |     0.0043266 |     0.0013153 |  0.6250284 |    0.0043172 |        0.0000687 |      2.6e-06 |     2.15e-05 |     0.0004359 |      2.60e-05 |       1.3e-06 | 0.0000524 |   0.0001573 | 8.80e-06 |      0.0005450 |      0.0371768 |    0.0035022 |     0.0014189 |   0.0007518 |       1.79e-05 |   0.0014739 |       0.0011206 |                            0.0004307 |                             0.0004457 |                           0.0004802 |                                    4.20e-06 |    0.1149734 |
| ShotgunWGS-TomatoPig18GutMicrobiome-Day7  | 18  | Tomato  | Day 7       | Tomato Day 7          |     0.0007724 |      0.1016953 |   0.0000567 | 0.0003691 |  0.0004482 |        3.59e-05 |     0.3872593 |      6.68e-05 |                9.50e-06 |  0.0001245 | 0.0012329 |   0.0018259 |   0.0001234 |       5.10e-05 |     0.0001537 |     0.0041846 |       0.0002496 |           0.0006182 |   0.0002373 |     0.0001290 |     0.0027879 |     0.0011617 |  0.4453529 |    0.0028428 |        0.0001167 |      1.1e-06 |     1.40e-05 |     0.0001616 |      1.80e-05 |       6.0e-07 | 0.0000752 |   0.0001374 | 1.63e-05 |      0.0006148 |      0.0389241 |    0.0026521 |     0.0008313 |   0.0005133 |       1.46e-05 |   0.0009592 |       0.0013244 |                            0.0008027 |                             0.0002906 |                           0.0007107 |                                    3.31e-05 |    0.1673878 |

Add column B to F

``` r
RelAbund.Phyla.Filt.zerofilt.withother.BtoF <- RelAbund.Phyla.Filt.zerofilt.withother %>%
  mutate(BtoF = Bacteroidetes/Firmicutes)

kable(head(RelAbund.Phyla.Filt.zerofilt.withother.BtoF))
```

| Sample\_Name                              | Pig | Diet    | Time\_Point | Diet\_By\_Time\_Point | Acidobacteria | Actinobacteria | Apicomplexa | Aquificae | Ascomycota | Bacillariophyta | Bacteroidetes | Basidiomycota | Candidatus Poribacteria | Chlamydiae |  Chlorobi | Chloroflexi | Chlorophyta | Chrysiogenetes | Crenarchaeota | Cyanobacteria | Deferribacteres | Deinococcus-Thermus | Dictyoglomi | Elusimicrobia | Euryarchaeota | Fibrobacteres | Firmicutes | Fusobacteria | Gemmatimonadetes | Hemichordata | Korarchaeota | Lentisphaerae | Microsporidia | Nanoarchaeota |  Nematoda | Nitrospirae | Placozoa | Planctomycetes | Proteobacteria | Spirochaetes | Synergistetes | Tenericutes | Thaumarchaeota | Thermotogae | Verrucomicrobia | unclassified (derived from Bacteria) | unclassified (derived from Eukaryota) | unclassified (derived from Viruses) | unclassified (derived from other sequences) | Other\_phyla |      BtoF |
|:------------------------------------------|:----|:--------|:------------|:----------------------|--------------:|---------------:|------------:|----------:|-----------:|----------------:|--------------:|--------------:|------------------------:|-----------:|----------:|------------:|------------:|---------------:|--------------:|--------------:|----------------:|--------------------:|------------:|--------------:|--------------:|--------------:|-----------:|-------------:|-----------------:|-------------:|-------------:|--------------:|--------------:|--------------:|----------:|------------:|---------:|---------------:|---------------:|-------------:|--------------:|------------:|---------------:|------------:|----------------:|-------------------------------------:|--------------------------------------:|------------------------------------:|--------------------------------------------:|-------------:|----------:|
| ShotgunWGS-ControlPig6GutMicrobiome-Day14 | 6   | Control | Day 14      | Control Day 14        |     0.0007406 |      0.0481336 |   0.0000595 | 0.0005033 |  0.0003842 |        2.71e-05 |     0.3670956 |      6.18e-05 |                6.70e-06 |  0.0001422 | 0.0012598 |   0.0020208 |   0.0000953 |       8.53e-05 |     0.0001670 |     0.0022775 |       0.0003850 |           0.0006393 |   0.0003136 |     0.0001629 |     0.0033951 |     0.0012287 |  0.5308271 |    0.0039555 |        0.0000544 |      6.7e-06 |     1.93e-05 |     0.0001971 |      1.26e-05 |       1.0e-06 | 0.0000459 |   0.0001420 | 1.75e-05 |      0.0003925 |      0.0271370 |    0.0029683 |     0.0011475 |   0.0007123 |       1.44e-05 |   0.0012921 |       0.0008269 |                            0.0003085 |                             0.0003247 |                           0.0003984 |                                    1.24e-05 |    0.1020763 | 0.6915539 |
| ShotgunWGS-ControlPig8GutMicrobiome-Day0  | 8   | Control | Day 0       | Control Day 0         |     0.0008850 |      0.0659807 |   0.0000914 | 0.0005366 |  0.0005228 |        4.38e-05 |     0.3312767 |      9.64e-05 |                1.17e-05 |  0.0001974 | 0.0014459 |   0.0022885 |   0.0001359 |       9.90e-05 |     0.0002148 |     0.0032408 |       0.0004747 |           0.0008538 |   0.0003700 |     0.0002398 |     0.0045655 |     0.0016578 |  0.5293436 |    0.0050574 |        0.0000809 |      4.3e-06 |     2.31e-05 |     0.0004364 |      2.10e-05 |       1.2e-06 | 0.0000631 |   0.0001898 | 1.62e-05 |      0.0006266 |      0.0368314 |    0.0041577 |     0.0017830 |   0.0008883 |       1.90e-05 |   0.0016916 |       0.0014957 |                            0.0004095 |                             0.0009973 |                           0.0006252 |                                    7.90e-06 |    0.1393790 | 0.6258254 |
| ShotgunWGS-ControlPig3GutMicrobiome-Day14 | 3   | Control | Day 14      | Control Day 14        |     0.0006897 |      0.0326727 |   0.0000492 | 0.0004253 |  0.0003318 |        3.34e-05 |     0.3263817 |      5.13e-05 |                5.40e-06 |  0.0001435 | 0.0011574 |   0.0019673 |   0.0000938 |       7.02e-05 |     0.0001652 |     0.0029204 |       0.0003675 |           0.0006011 |   0.0002872 |     0.0001810 |     0.0033125 |     0.0012911 |  0.5870258 |    0.0037180 |        0.0000624 |      4.1e-06 |     1.35e-05 |     0.0002005 |      1.45e-05 |       1.0e-06 | 0.0000378 |   0.0001458 | 8.00e-06 |      0.0004066 |      0.0271624 |    0.0028287 |     0.0010743 |   0.0006457 |       1.53e-05 |   0.0013084 |       0.0009681 |                            0.0002973 |                             0.0003279 |                           0.0005188 |                                    1.61e-05 |    0.0865920 | 0.5559920 |
| ShotgunWGS-TomatoPig14GutMicrobiome-Day7  | 14  | Tomato  | Day 7       | Tomato Day 7          |     0.0007324 |      0.0329202 |   0.0001398 | 0.0005384 |  0.0005593 |        5.41e-05 |     0.3461498 |      9.49e-05 |                1.41e-05 |  0.0001856 | 0.0013024 |   0.0023710 |   0.0001531 |       9.49e-05 |     0.0002730 |     0.0028512 |       0.0004519 |           0.0008048 |   0.0003529 |     0.0002988 |     0.0069066 |     0.0014564 |  0.5231174 |    0.0046147 |        0.0000508 |      4.2e-06 |     3.00e-05 |     0.0004095 |      2.00e-05 |       5.0e-06 | 0.0001148 |   0.0002006 | 1.58e-05 |      0.0004494 |      0.0597393 |    0.0046887 |     0.0015563 |   0.0010752 |       3.08e-05 |   0.0016628 |       0.0010478 |                            0.0009662 |                             0.0005509 |                           0.0008405 |                                    9.57e-05 |    0.1307244 | 0.6617057 |
| ShotgunWGS-ControlPig5GutMicrobiome-Day7  | 5   | Control | Day 7       | Control Day 7         |     0.0006564 |      0.0463444 |   0.0000557 | 0.0004617 |  0.0003835 |        3.61e-05 |     0.2599966 |      5.93e-05 |                7.20e-06 |  0.0001644 | 0.0012486 |   0.0020345 |   0.0001263 |       8.89e-05 |     0.0002123 |     0.0034222 |       0.0003796 |           0.0006821 |   0.0002920 |     0.0002165 |     0.0043266 |     0.0013153 |  0.6250284 |    0.0043172 |        0.0000687 |      2.6e-06 |     2.15e-05 |     0.0004359 |      2.60e-05 |       1.3e-06 | 0.0000524 |   0.0001573 | 8.80e-06 |      0.0005450 |      0.0371768 |    0.0035022 |     0.0014189 |   0.0007518 |       1.79e-05 |   0.0014739 |       0.0011206 |                            0.0004307 |                             0.0004457 |                           0.0004802 |                                    4.20e-06 |    0.1149734 | 0.4159757 |
| ShotgunWGS-TomatoPig18GutMicrobiome-Day7  | 18  | Tomato  | Day 7       | Tomato Day 7          |     0.0007724 |      0.1016953 |   0.0000567 | 0.0003691 |  0.0004482 |        3.59e-05 |     0.3872593 |      6.68e-05 |                9.50e-06 |  0.0001245 | 0.0012329 |   0.0018259 |   0.0001234 |       5.10e-05 |     0.0001537 |     0.0041846 |       0.0002496 |           0.0006182 |   0.0002373 |     0.0001290 |     0.0027879 |     0.0011617 |  0.4453529 |    0.0028428 |        0.0001167 |      1.1e-06 |     1.40e-05 |     0.0001616 |      1.80e-05 |       6.0e-07 | 0.0000752 |   0.0001374 | 1.63e-05 |      0.0006148 |      0.0389241 |    0.0026521 |     0.0008313 |   0.0005133 |       1.46e-05 |   0.0009592 |       0.0013244 |                            0.0008027 |                             0.0002906 |                           0.0007107 |                                    3.31e-05 |    0.1673878 | 0.8695561 |

Convert data from wide to long (i.e. make data
[tidy](https://r4ds.had.co.nz/tidy-data.html))

``` r
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF %>%
  pivot_longer(cols = 6:ncol(.),
               names_to = "phylum",
               values_to = "rel_abund")

RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long[1:10,]
```

    ## # A tibble: 10 × 7
    ##    Sample_Name          Pig   Diet  Time_Point Diet_By_Time_Po… phylum rel_abund
    ##    <chr>                <fct> <fct> <fct>      <fct>            <chr>      <dbl>
    ##  1 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Acido…   7.41e-4
    ##  2 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Actin…   4.81e-2
    ##  3 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Apico…   5.95e-5
    ##  4 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Aquif…   5.03e-4
    ##  5 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Ascom…   3.84e-4
    ##  6 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Bacil…   2.71e-5
    ##  7 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Bacte…   3.67e-1
    ##  8 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Basid…   6.18e-5
    ##  9 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Candi…   6.70e-6
    ## 10 ShotgunWGS-ControlP… 6     Cont… Day 14     Control Day 14   Chlam…   1.42e-4

### Plotting

Stacked bar chart of B, F, and all the other phyla

``` r
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.BFandOther <-
  RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long %>%
    filter(phylum %in% c("Bacteroidetes", "Firmicutes", "Other_phyla"))

head(RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.BFandOther)
```

    ## # A tibble: 6 × 7
    ##   Sample_Name           Pig   Diet  Time_Point Diet_By_Time_Po… phylum rel_abund
    ##   <chr>                 <fct> <fct> <fct>      <fct>            <chr>      <dbl>
    ## 1 ShotgunWGS-ControlPi… 6     Cont… Day 14     Control Day 14   Bacte…     0.367
    ## 2 ShotgunWGS-ControlPi… 6     Cont… Day 14     Control Day 14   Firmi…     0.531
    ## 3 ShotgunWGS-ControlPi… 6     Cont… Day 14     Control Day 14   Other…     0.102
    ## 4 ShotgunWGS-ControlPi… 8     Cont… Day 0      Control Day 0    Bacte…     0.331
    ## 5 ShotgunWGS-ControlPi… 8     Cont… Day 0      Control Day 0    Firmi…     0.529
    ## 6 ShotgunWGS-ControlPi… 8     Cont… Day 0      Control Day 0    Other…     0.139

``` r
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.BFandOther %>%
  ggplot(aes(x=as.numeric(Pig), y=rel_abund, fill=phylum))+
  geom_col()+
  scale_fill_brewer(palette = "GnBu") +
  facet_grid(~Time_Point)+
  theme_classic()+
  labs(y="Relative Abundance", 
       fill="Phylum",
       x = "Pig") +
  theme(panel.grid = element_blank(), axis.text = element_text(color="black"),
        strip.text = element_text(color = "black", size = 14), 
        strip.background = element_blank())
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-171-1.png)<!-- -->

### B to F Ratio

#### Plotting

B to F boxplot with jitter

``` r
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.BtoF <- 
  RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long %>%
  filter(phylum == "BtoF")

head(RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.BtoF)
```

    ## # A tibble: 6 × 7
    ##   Sample_Name           Pig   Diet  Time_Point Diet_By_Time_Po… phylum rel_abund
    ##   <chr>                 <fct> <fct> <fct>      <fct>            <chr>      <dbl>
    ## 1 ShotgunWGS-ControlPi… 6     Cont… Day 14     Control Day 14   BtoF       0.692
    ## 2 ShotgunWGS-ControlPi… 8     Cont… Day 0      Control Day 0    BtoF       0.626
    ## 3 ShotgunWGS-ControlPi… 3     Cont… Day 14     Control Day 14   BtoF       0.556
    ## 4 ShotgunWGS-TomatoPig… 14    Toma… Day 7      Tomato Day 7     BtoF       0.662
    ## 5 ShotgunWGS-ControlPi… 5     Cont… Day 7      Control Day 7    BtoF       0.416
    ## 6 ShotgunWGS-TomatoPig… 18    Toma… Day 7      Tomato Day 7     BtoF       0.870

``` r
BtoF_Boxplot <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.BtoF %>%
  ggplot(aes(x=Diet, y=rel_abund, fill=Diet_By_Time_Point))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(fill = Diet_By_Time_Point), color = "black", alpha = 0.7, position=position_jitterdodge()) +
  scale_fill_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  scale_color_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(color = "black"), 
        panel.grid.minor = element_blank()) +
  labs(x=NULL, 
       y= "Bacteroidetes to Firmicutes", 
       fill="Diet & Time Point",
       title = "Ratio of Bacteroidetes to Firmicutes") 

BtoF_Boxplot
```

    ## Warning: Removed 5 rows containing non-finite values (stat_boxplot).

    ## Warning: Removed 5 rows containing missing values (geom_point).

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-172-1.png)<!-- -->

Saving

``` r
ggsave("Figures/BtoFRatio_Boxplot.png", 
       plot = BtoF_Boxplot, 
       dpi = 800, 
       width = 7, 
       height = 5)
```

#### Statistics

``` r
head(RelAbund.Phyla.Filt.zerofilt.withother.BtoF)
```

    ## # A tibble: 6 × 52
    ##   Sample_Name              Pig   Diet  Time_Point Diet_By_Time_Po… Acidobacteria
    ##   <chr>                    <fct> <fct> <fct>      <fct>                    <dbl>
    ## 1 ShotgunWGS-ControlPig6G… 6     Cont… Day 14     Control Day 14        0.000741
    ## 2 ShotgunWGS-ControlPig8G… 8     Cont… Day 0      Control Day 0         0.000885
    ## 3 ShotgunWGS-ControlPig3G… 3     Cont… Day 14     Control Day 14        0.000690
    ## 4 ShotgunWGS-TomatoPig14G… 14    Toma… Day 7      Tomato Day 7          0.000732
    ## 5 ShotgunWGS-ControlPig5G… 5     Cont… Day 7      Control Day 7         0.000656
    ## 6 ShotgunWGS-TomatoPig18G… 18    Toma… Day 7      Tomato Day 7          0.000772
    ## # … with 46 more variables: Actinobacteria <dbl>, Apicomplexa <dbl>,
    ## #   Aquificae <dbl>, Ascomycota <dbl>, Bacillariophyta <dbl>,
    ## #   Bacteroidetes <dbl>, Basidiomycota <dbl>, `Candidatus Poribacteria` <dbl>,
    ## #   Chlamydiae <dbl>, Chlorobi <dbl>, Chloroflexi <dbl>, Chlorophyta <dbl>,
    ## #   Chrysiogenetes <dbl>, Crenarchaeota <dbl>, Cyanobacteria <dbl>,
    ## #   Deferribacteres <dbl>, `Deinococcus-Thermus` <dbl>, Dictyoglomi <dbl>,
    ## #   Elusimicrobia <dbl>, Euryarchaeota <dbl>, Fibrobacteres <dbl>, …

``` r
# select only columns used for ANOVA
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.ForANOVA <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF %>%
  select(Pig, Diet, Time_Point, BtoF)
```

Repeated measures ANOVA of B to F ratio

``` r
BtoF.Ratio.ANOVA <- anova_test(data = RelAbund.Phyla.Filt.zerofilt.withother.BtoF.ForANOVA, 
                          formula = BtoF ~ Diet*Time_Point + Error(Pig/(Time_Point)),
                          dv = BtoF, wid = Pig, between = Diet, within = Time_Point)

get_anova_table(BtoF.Ratio.ANOVA)
```

    ## ANOVA Table (type II tests)
    ## 
    ##            Effect DFn DFd     F     p p<.05   ges
    ## 1            Diet   1  18 0.125 0.728       0.004
    ## 2      Time_Point   2  36 5.437 0.009     * 0.113
    ## 3 Diet:Time_Point   2  36 0.850 0.436       0.020

-   Significant effect of time point (p = 0.009)  
-   Nonsignificant effect of diet (p = 0.728)  
-   Nonsignificant effect of diet:timepoint (p = 0.436)

Use posthoc to see where is significant using a fdr p-value adjustment
for multiple testing, grouping by time (both diets)

``` r
BtoF.Ratio.ANOVA.posthoc <- pairwise_t_test(BtoF ~ Time_Point, 
                                 data = RelAbund.Phyla.Filt.zerofilt.withother.BtoF, 
                                 paired = TRUE, 
                                 p.adjust.method = "fdr") 

BtoF.Ratio.ANOVA.posthoc
```

    ## # A tibble: 3 × 10
    ##   .y.   group1 group2    n1    n2 statistic    df     p p.adj p.adj.signif
    ## * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl> <chr>       
    ## 1 BtoF  Day 0  Day 7     20    20      1.42    19 0.172 0.235 ns          
    ## 2 BtoF  Day 0  Day 14    20    20      3.15    19 0.005 0.016 *           
    ## 3 BtoF  Day 7  Day 14    20    20      1.23    19 0.235 0.235 ns

Significant difference is between day 0 and day 14 (padj = 0.016).

Use posthoc to see where is significant using a fdr p-value adjustment
for multiple testing, separating by diet

``` r
BtoF.Ratio.ANOVA.posthoc.bytime <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF.ForANOVA %>%
  group_by(Diet) %>%
  pairwise_t_test(BtoF ~ Time_Point, 
                  paired = TRUE, 
                  p.adjust.method = "fdr")

BtoF.Ratio.ANOVA.posthoc.bytime
```

    ## # A tibble: 6 × 11
    ##   Diet  .y.   group1 group2    n1    n2 statistic    df     p p.adj p.adj.signif
    ## * <fct> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl> <chr>       
    ## 1 Cont… BtoF  Day 0  Day 7     10    10     1.34      9 0.213 0.32  ns          
    ## 2 Cont… BtoF  Day 0  Day 14    10    10     3.19      9 0.011 0.033 *           
    ## 3 Cont… BtoF  Day 7  Day 14    10    10     0.589     9 0.571 0.571 ns          
    ## 4 Toma… BtoF  Day 0  Day 7     10    10     0.403     9 0.696 0.696 ns          
    ## 5 Toma… BtoF  Day 0  Day 14    10    10     1.35      9 0.211 0.633 ns          
    ## 6 Toma… BtoF  Day 7  Day 14    10    10     0.812     9 0.438 0.657 ns

Significant in control between day 0 and 14, padj = 0.033. All else
non-significant.

### Bacteroidetes and Firmicutes

#### Plotting

Boxplotting

``` r
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.OnlyBandF <- 
  RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long %>%
    filter(phylum == "Bacteroidetes" | phylum == "Firmicutes")

BandF_Boxplot <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF.long.OnlyBandF %>%
  ggplot(aes(x=Diet, y=rel_abund, fill=Diet_By_Time_Point))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(fill = Diet_By_Time_Point), color = "black", position=position_jitterdodge(), alpha = 0.7) +
  scale_fill_manual(values=c("skyblue1", "dodgerblue", "royalblue4", "sienna1","firebrick3","tomato4"))+
  ylim(0, 0.75) +
  theme_bw() +
  facet_wrap("phylum")+
  labs(x=NULL, y= "Relative Abundance", fill="Diet & Time Point") +
  theme(axis.text.x = element_text(size = 11, color = "black"), 
        axis.text.y = element_text(color = "black"), 
        panel.grid.minor = element_blank(), 
        strip.text.x = element_text(color = "black", size = 14),
        strip.background = element_rect(fill = "white"))

BandF_Boxplot
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-178-1.png)<!-- -->

Saving

``` r
ggsave("Figures/BacteroidetesAndFirmicutes_Boxplot.png", 
       plot = BandF_Boxplot, 
       dpi = 800, 
       width = 7, 
       height = 5)
```

#### Statistics

##### Bacteroidetes

``` r
# select only columns used for ANOVA
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.Bonly <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF %>%
  select(Pig, Diet, Time_Point, Bacteroidetes)
```

Bacteroidetes repeated measures ANOVA

``` r
Bacteroidetes.ANOVA <- anova_test(data = RelAbund.Phyla.Filt.zerofilt.withother.BtoF.Bonly, 
                          formula = Bacteroidetes ~ Diet*Time_Point + Error(Pig/(Time_Point)),
                          dv = Bacteroidetes, wid = Pig, between = Diet, within = Time_Point)

get_anova_table(Bacteroidetes.ANOVA)
```

    ## ANOVA Table (type II tests)
    ## 
    ##            Effect DFn DFd     F     p p<.05      ges
    ## 1            Diet   1  18 0.009 0.928       0.000272
    ## 2      Time_Point   2  36 4.131 0.024     * 0.089000
    ## 3 Diet:Time_Point   2  36 0.700 0.503       0.016000

-   Significant effect of time point (p = 0.024)  
-   Nonsignificant effect of diet (p = 0.928)  
-   onsignificant effect of diet:timepoint (p = 0.503)

Use posthoc to see where is significant using a fdr p-value adjustment
for multiple testing, grouping by time (both diets)

``` r
Bacteroidetes.ANOVA.posthoc <- pairwise_t_test(Bacteroidetes ~ Time_Point, 
                                 data = RelAbund.Phyla.Filt.zerofilt.withother.BtoF, 
                                 paired = TRUE, 
                                 p.adjust.method = "fdr") 

Bacteroidetes.ANOVA.posthoc
```

    ## # A tibble: 3 × 10
    ##   .y.         group1 group2    n1    n2 statistic    df     p p.adj p.adj.signif
    ## * <chr>       <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl> <chr>       
    ## 1 Bacteroide… Day 0  Day 7     20    20     1.62     19 0.122 0.183 ns          
    ## 2 Bacteroide… Day 0  Day 14    20    20     2.56     19 0.019 0.058 ns          
    ## 3 Bacteroide… Day 7  Day 14    20    20     0.572    19 0.574 0.574 ns

Borderline significant difference is between day 0 and day 14 (padj =
0.058).

Use posthoc to see where is significant using a fdr p-value adjustment
for multiple testing, separating by diet

``` r
Bacteroidetes.ANOVA.posthoc.bytime <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF.Bonly %>%
  group_by(Diet) %>%
  pairwise_t_test(Bacteroidetes ~ Time_Point, 
                  paired = TRUE, 
                  p.adjust.method = "fdr")

Bacteroidetes.ANOVA.posthoc.bytime
```

    ## # A tibble: 6 × 11
    ##   Diet  .y.   group1 group2    n1    n2 statistic    df     p p.adj p.adj.signif
    ## * <fct> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl> <chr>       
    ## 1 Cont… Bact… Day 0  Day 7     10    10     1.56      9 0.154 0.231 ns          
    ## 2 Cont… Bact… Day 0  Day 14    10    10     3.02      9 0.015 0.044 *           
    ## 3 Cont… Bact… Day 7  Day 14    10    10     0.152     9 0.883 0.883 ns          
    ## 4 Toma… Bact… Day 0  Day 7     10    10     0.486     9 0.638 0.638 ns          
    ## 5 Toma… Bact… Day 0  Day 14    10    10     1.08      9 0.308 0.638 ns          
    ## 6 Toma… Bact… Day 7  Day 14    10    10     0.495     9 0.633 0.638 ns

Significant difference in control between day 0 and 14, padj = 0.044.
All else non-significant.

##### Firmicutes

``` r
# select only columns used for ANOVA
RelAbund.Phyla.Filt.zerofilt.withother.BtoF.Fonly <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF %>%
  select(Pig, Diet, Time_Point, Firmicutes)
```

Firmicutes repeated measures ANOVA

``` r
Firmicutes.ANOVA <- anova_test(data = RelAbund.Phyla.Filt.zerofilt.withother.BtoF.Fonly, 
                          formula = Firmicutes ~ Diet*Time_Point + Error(Pig/(Time_Point)),
                          dv = Firmicutes, wid = Pig, between = Diet, within = Time_Point)

get_anova_table(Firmicutes.ANOVA)
```

    ## ANOVA Table (type II tests)
    ## 
    ##            Effect DFn DFd     F     p p<.05   ges
    ## 1            Diet   1  18 1.079 0.313       0.033
    ## 2      Time_Point   2  36 8.102 0.001     * 0.161
    ## 3 Diet:Time_Point   2  36 0.993 0.380       0.023

-   Significant effect of time point (p = 0.001)  
-   Nonsignificant effect of diet (p = 0.313)  
-   Nonsignificant effect of diet:timepoint (p = 0.380)

Use posthoc to see where is significant using a fdr p-value adjustment
for multiple testing, grouping by time (both diets)

``` r
Firmicutes.ANOVA.posthoc <- pairwise_t_test(Firmicutes ~ Time_Point, 
                                 data = RelAbund.Phyla.Filt.zerofilt.withother.BtoF, 
                                 paired = TRUE, 
                                 p.adjust.method = "fdr") 

Firmicutes.ANOVA.posthoc
```

    ## # A tibble: 3 × 10
    ##   .y.       group1 group2    n1    n2 statistic    df       p p.adj p.adj.signif
    ## * <chr>     <chr>  <chr>  <int> <int>     <dbl> <dbl>   <dbl> <dbl> <chr>       
    ## 1 Firmicut… Day 0  Day 7     20    20     -1.60    19 1.25e-1 0.125 ns          
    ## 2 Firmicut… Day 0  Day 14    20    20     -3.93    19 9.07e-4 0.003 **          
    ## 3 Firmicut… Day 7  Day 14    20    20     -1.67    19 1.11e-1 0.125 ns

Significant difference is between day 0 and day 14 (padj = 0.003).

Use posthoc to see where is significant using a fdr p-value adjustment
for multiple testing, separating by diet

``` r
Firmicutes.ANOVA.posthoc.bytime <- RelAbund.Phyla.Filt.zerofilt.withother.BtoF.Fonly %>%
  group_by(Diet) %>%
  pairwise_t_test(Firmicutes ~ Time_Point, 
                  paired = TRUE, 
                  p.adjust.method = "fdr")

Firmicutes.ANOVA.posthoc.bytime
```

    ## # A tibble: 6 × 11
    ##   Diet  .y.   group1 group2    n1    n2 statistic    df     p p.adj p.adj.signif
    ## * <fct> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl> <chr>       
    ## 1 Cont… Firm… Day 0  Day 7     10    10    -1.48      9 0.173 0.259 ns          
    ## 2 Cont… Firm… Day 0  Day 14    10    10    -3.25      9 0.01  0.03  *           
    ## 3 Cont… Firm… Day 7  Day 14    10    10    -0.843     9 0.421 0.421 ns          
    ## 4 Toma… Firm… Day 0  Day 7     10    10    -0.550     9 0.596 0.596 ns          
    ## 5 Toma… Firm… Day 0  Day 14    10    10    -1.96      9 0.082 0.246 ns          
    ## 6 Toma… Firm… Day 7  Day 14    10    10    -1.05      9 0.322 0.483 ns

Significant difference is in control between day 0 and 14, padj = 0.030.
All else non-significant.

## Alpha Diversity

Calculated alpha diversity of phyla based on relative abundance,
including all the filtering for implausible phyla and removing samples
with more than 33.33% missing samples

### Wrangling

``` r
dim(RelAbund.Phyla.Filt.zerofilt)
```

    ## [1] 60 50

``` r
RelAbund.Phyla.Filt.zerofilt[1:5,1:10]
```

    ## # A tibble: 5 × 10
    ##   Sample_Name              Pig   Diet  Time_Point Diet_By_Time_Po… Acidobacteria
    ##   <chr>                    <fct> <fct> <fct>      <fct>                    <dbl>
    ## 1 ShotgunWGS-ControlPig6G… 6     Cont… Day 14     Control Day 14        0.000741
    ## 2 ShotgunWGS-ControlPig8G… 8     Cont… Day 0      Control Day 0         0.000885
    ## 3 ShotgunWGS-ControlPig3G… 3     Cont… Day 14     Control Day 14        0.000690
    ## 4 ShotgunWGS-TomatoPig14G… 14    Toma… Day 7      Tomato Day 7          0.000732
    ## 5 ShotgunWGS-ControlPig5G… 5     Cont… Day 7      Control Day 7         0.000656
    ## # … with 4 more variables: Actinobacteria <dbl>, Apicomplexa <dbl>,
    ## #   Aquificae <dbl>, Ascomycota <dbl>

Wrangle

``` r
# move Sample_Name to rownames, remove metadata
RelAbund.Phyla.Filt.zerofilt.alphadiv <- RelAbund.Phyla.Filt.zerofilt

rownames(RelAbund.Phyla.Filt.zerofilt.alphadiv) <- RelAbund.Phyla.Filt.zerofilt.alphadiv$Sample_Name  
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
# remove metadata
RelAbund.Phyla.Filt.zerofilt.alphadiv <- RelAbund.Phyla.Filt.zerofilt.alphadiv %>%
  select(Acidobacteria:ncol(.))

RelAbund.Phyla.Filt.zerofilt.alphadiv[1:5,1:5]
```

    ## # A tibble: 5 × 5
    ##   Acidobacteria Actinobacteria Apicomplexa Aquificae Ascomycota
    ##           <dbl>          <dbl>       <dbl>     <dbl>      <dbl>
    ## 1      0.000741         0.0481   0.0000595  0.000503   0.000384
    ## 2      0.000885         0.0660   0.0000914  0.000537   0.000523
    ## 3      0.000690         0.0327   0.0000492  0.000425   0.000332
    ## 4      0.000732         0.0329   0.000140   0.000538   0.000559
    ## 5      0.000656         0.0463   0.0000557  0.000462   0.000384

### Calculate alpha diversity

``` r
# run alpha diversity on phyla
phyla.filt.div <- diversity(RelAbund.Phyla.Filt.zerofilt.alphadiv, index = "shannon")

# convert to df
phyla.filt.div.df <- as.data.frame(phyla.filt.div)

# make column name 'shannon.phyla.filt'
colnames(phyla.filt.div.df) <- "shannon.phyla.filt"

head(phyla.filt.div.df)
```

    ##   shannon.phyla.filt
    ## 1           1.122832
    ## 2           1.231877
    ## 3           1.062241
    ## 4           1.221226
    ## 5           1.108308
    ## 6           1.261677

Combine with metadata

``` r
# combine with metadata
phyla.filt.div.df.meta <- cbind(RelAbund.Phyla.Filt.zerofilt[,1:5], phyla.filt.div.df)

head(phyla.filt.div.df.meta)
```

    ##                                 Sample_Name Pig    Diet Time_Point
    ## 1 ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## 2  ShotgunWGS-ControlPig8GutMicrobiome-Day0   8 Control      Day 0
    ## 3 ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ## 4  ShotgunWGS-TomatoPig14GutMicrobiome-Day7  14  Tomato      Day 7
    ## 5  ShotgunWGS-ControlPig5GutMicrobiome-Day7   5 Control      Day 7
    ## 6  ShotgunWGS-TomatoPig18GutMicrobiome-Day7  18  Tomato      Day 7
    ##   Diet_By_Time_Point shannon.phyla.filt
    ## 1     Control Day 14           1.122832
    ## 2      Control Day 0           1.231877
    ## 3     Control Day 14           1.062241
    ## 4       Tomato Day 7           1.221226
    ## 5      Control Day 7           1.108308
    ## 6       Tomato Day 7           1.261677

### Plotting

X-axis by diet

``` r
alpha.diversity.phyla.bydiet <- phyla.filt.div.df.meta %>%
  ggplot(aes(x = Diet, y = shannon.phyla.filt, fill = Diet_By_Time_Point)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(fill = Diet_By_Time_Point), 
             color = "black", 
             alpha = 0.7, 
             position=position_jitterdodge()) +
  scale_fill_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  labs(x=NULL, 
       y="Shannon diversity index", 
       title = "Alpha Diversity",
       subtitle = "Shannon Index, Phyla Level", 
       fill="Diet & Time Point")

alpha.diversity.phyla.bydiet
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-192-1.png)<!-- -->

``` r
ggsave("Figures/AlphaDiversityPhyla_ByDiet_Boxplot.png", 
       plot = alpha.diversity.phyla.bydiet, 
       dpi = 800, 
       width = 10, 
       height = 6)
```

X-axis by time point

``` r
alpha.diversity.phyla.bytime <- phyla.filt.div.df.meta %>%
  ggplot(aes(x = Time_Point, y = shannon.phyla.filt, fill = Diet_By_Time_Point)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(color = "black", alpha = 0.7, position=position_jitterdodge()) +
  scale_fill_manual(values = c("skyblue1", "dodgerblue", "royalblue4", 
                               "sienna1","firebrick3","tomato4")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  labs(x=NULL, 
       y="Shannon diversity index", 
       title = "Alpha Diversity",
       subtitle = "Shannon Index, Phyla Level", 
       fill="Diet & Time Point")

alpha.diversity.phyla.bytime
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-194-1.png)<!-- -->

``` r
ggsave("Figures/AlphaDiversityPhyla_ByTime_Boxplot.png", 
       plot = alpha.diversity.phyla.bytime, 
       dpi = 800, 
       width = 7, 
       height = 5)
```

### Statistics

Repeated measures ANOVA

``` r
# must remove columns that aren't used in anova
head(phyla.filt.div.df.meta) 
```

    ##                                 Sample_Name Pig    Diet Time_Point
    ## 1 ShotgunWGS-ControlPig6GutMicrobiome-Day14   6 Control     Day 14
    ## 2  ShotgunWGS-ControlPig8GutMicrobiome-Day0   8 Control      Day 0
    ## 3 ShotgunWGS-ControlPig3GutMicrobiome-Day14   3 Control     Day 14
    ## 4  ShotgunWGS-TomatoPig14GutMicrobiome-Day7  14  Tomato      Day 7
    ## 5  ShotgunWGS-ControlPig5GutMicrobiome-Day7   5 Control      Day 7
    ## 6  ShotgunWGS-TomatoPig18GutMicrobiome-Day7  18  Tomato      Day 7
    ##   Diet_By_Time_Point shannon.phyla.filt
    ## 1     Control Day 14           1.122832
    ## 2      Control Day 0           1.231877
    ## 3     Control Day 14           1.062241
    ## 4       Tomato Day 7           1.221226
    ## 5      Control Day 7           1.108308
    ## 6       Tomato Day 7           1.261677

``` r
phyla.filt.div.foranova <- phyla.filt.div.df.meta[,-c(1,5)]

head(phyla.filt.div.foranova)
```

    ##   Pig    Diet Time_Point shannon.phyla.filt
    ## 1   6 Control     Day 14           1.122832
    ## 2   8 Control      Day 0           1.231877
    ## 3   3 Control     Day 14           1.062241
    ## 4  14  Tomato      Day 7           1.221226
    ## 5   5 Control      Day 7           1.108308
    ## 6  18  Tomato      Day 7           1.261677

``` r
phyla.filt.alphadiv.anova <- 
  anova_test(data = phyla.filt.div.foranova,
             formula = shannon.phyla.filt ~ Diet*Time_Point + Error(Pig/Time_Point),
             dv = shannon.phyla.filt, 
             wid = Pig, 
             within = Time_Point, 
             between = Diet)

get_anova_table(phyla.filt.alphadiv.anova)
```

    ## ANOVA Table (type II tests)
    ## 
    ##            Effect DFn DFd      F     p p<.05   ges
    ## 1            Diet   1  18 10.767 0.004     * 0.158
    ## 2      Time_Point   2  36  2.628 0.086       0.091
    ## 3 Diet:Time_Point   2  36  0.236 0.791       0.009

-   Significant effect of diet (p = 0.004)  
-   Non-significant effect of timepoint (0.086)  
-   Non-significant interaction of diet:time point (p = 0.791).

Check for normality

``` r
shapiro.test(phyla.filt.div.df.meta$shannon.phyla.filt)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  phyla.filt.div.df.meta$shannon.phyla.filt
    ## W = 0.99022, p-value = 0.9135

Normal.

Post-hoc tests

``` r
posthoc.morevariables <- phyla.filt.div.foranova %>%
  group_by(Time_Point) %>%
  anova_test(dv = shannon.phyla.filt, wid = Pig, between = Diet) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "fdr")
```

    ## Coefficient covariances computed by hccm()
    ## Coefficient covariances computed by hccm()
    ## Coefficient covariances computed by hccm()

``` r
posthoc.morevariables
```

    ## # A tibble: 3 × 9
    ##   Time_Point Effect   DFn   DFd     F     p `p<.05`   ges p.adj
    ##   <fct>      <chr>  <dbl> <dbl> <dbl> <dbl> <chr>   <dbl> <dbl>
    ## 1 Day 0      Diet       1    18  1.22 0.284 ""      0.063 0.284
    ## 2 Day 7      Diet       1    18  3.65 0.072 ""      0.169 0.108
    ## 3 Day 14     Diet       1    18  8.74 0.008 "*"     0.327 0.024

Significant effect of diet at day 14 (padj = 0.024)

``` r
posthoc.evenmorespecific <- phyla.filt.div.foranova %>%
  group_by(Time_Point) %>%
  pairwise_t_test(shannon.phyla.filt ~ Diet,
                  paired = TRUE,
                  p.adjust.method = "fdr")

posthoc.evenmorespecific
```

    ## # A tibble: 3 × 11
    ##   Time_Point .y.           group1 group2    n1    n2 statistic    df     p p.adj
    ## * <fct>      <chr>         <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl>
    ## 1 Day 0      shannon.phyl… Contr… Tomato    10    10     -1.03     9 0.328 0.328
    ## 2 Day 7      shannon.phyl… Contr… Tomato    10    10     -1.93     9 0.085 0.085
    ## 3 Day 14     shannon.phyl… Contr… Tomato    10    10     -3.21     9 0.011 0.011
    ## # … with 1 more variable: p.adj.signif <chr>

-   Significant effect between control and tomato at day 14 (padj =
    0.01)  
-   Nonsignificant at day 0 (padj = 0.328)  
-   Nonsignificant but getting close at day 7 (padj = 0.085)

## ALDEx2

Quick introduction to anatomy of the aldex function

The aldex function does every step - data transformation and
statistics  
variable.name &lt;- aldex(reads.data, variables.vector, mc.samples=\#,
test=“t”/“kw”, effect=T/F)  
reads.data - your reads/count data, un changed  
variables.vector - a vector of the variables corresponding to sample
groups, in SAME order as sample names (and therefore columns)  
mc.samples - here you tell the function how many Monte Carlo sampels to
use with an integer (128 is typical)  
test - which test do you want, t-test and wilcoxon, or anova-like and
kruskal wallace? (will always do the parametric and non-parametric) t =
t-test and wilcoxon kw = anova-like and kruskal wallace  
effect - do you want it to incude effect results in output?

Key to aldex outputs - taken directly from vignette

-   we.ep - Expected P value of Welch’s t test
-   we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t
    test  
-   wi.ep - Expected P value of Wilcoxon rank test
-   wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon
    test
-   kw.ep - Expected P value of Kruskal-Wallace test
-   kw.eBH - Expected Benjamini-Hochberg corrected P value of
    Kruskal-Wallace test
-   glm.ep - Expected P value of glm test
-   glm.eBH - Expected Benjamini-Hochberg corrected P value of glm test
-   rab.all - median clr value for all samples in the feature
-   rab.win.NS - median clr value for the NS group of samples
-   rab.win.S - median clr value for the S group of samples
-   dif.btw - median difference in clr values between S and NS groups
-   dif.win - median of the largest difference in clr values within S
    and NS groups
-   effect - median effect size: diff.btw / max(diff.win) for all
    instances
-   overlap - proportion of effect size that overlaps 0 (i.e. no effect)

ALDEx2 takes counts, not relative abundance.

We are using Benjamini Hochberg corrected pvalues, or `we.eBH` for
t-tests (i.e., subsetting by time), and Benjamini-Hochberg corrected
pvalues of the glm test `glm.eBH` for ANOVA tests (i.e., subsetting by
diet)

Downloading ALDEx2

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")
```

### Wrangling

Since we use counts for ALDEx2, we need to filter our counts data to
include only the phyla we ended up using in our final analysis

``` r
# this data set filtered to remove inplausible phyla, but still includes phyla with a lot of missing values 
Phyla.Counts.Filt[1:10,1:10]
```

    ## # A tibble: 10 × 10
    ##    domain    phylum           `ShotgunWGS-Co…` `ShotgunWGS-Co…` `ShotgunWGS-Co…`
    ##    <chr>     <chr>                       <dbl>            <dbl>            <dbl>
    ##  1 Bacteria  Acidobacteria                2874             3717             2663
    ##  2 Bacteria  Actinobacteria             186789           277130           126155
    ##  3 Eukaryota Apicomplexa                   231              384              190
    ##  4 Bacteria  Aquificae                    1953             2254             1642
    ##  5 Eukaryota Ascomycota                   1491             2196             1281
    ##  6 Eukaryota Bacillariophyta               105              184              129
    ##  7 Bacteria  Bacteroidetes             1424565          1391417          1260217
    ##  8 Eukaryota Basidiomycota                 240              405              198
    ##  9 Eukaryota Blastocladiomyc…                0                0                0
    ## 10 Bacteria  Candidatus Pori…               26               49               21
    ## # … with 5 more variables: `ShotgunWGS-TomatoPig14GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-ControlPig5GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-TomatoPig18GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-TomatoPig16GutMicrobiome-Day7` <dbl>,
    ## #   `ShotgunWGS-ControlPig10GutMicrobiome-Day7` <dbl>

``` r
dim(Phyla.Counts.Filt)
```

    ## [1] 53 62

``` r
# final phyla list, after filtering for number of zeros
final_phyla[1:10,]
```

    ##  [1] "Acidobacteria"           "Actinobacteria"         
    ##  [3] "Apicomplexa"             "Aquificae"              
    ##  [5] "Ascomycota"              "Bacillariophyta"        
    ##  [7] "Bacteroidetes"           "Basidiomycota"          
    ##  [9] "Candidatus Poribacteria" "Chlamydiae"

``` r
# how many final phyla do we have?
dim(final_phyla)
```

    ## [1] 45  1

``` r
# join to create a df with phyla in rows, samples in columns
# filtered for genera used in this analysis
phyla_counts_foraldex <- inner_join(final_phyla, Phyla.Counts.Filt,
                                    by = "phylum")

dim(phyla_counts_foraldex)
```

    ## [1] 45 62

``` r
phyla_counts_foraldex[1:10, 1:4]
```

    ##                     phylum    domain ShotgunWGS-ControlPig6GutMicrobiome-Day14
    ## 1            Acidobacteria  Bacteria                                      2874
    ## 2           Actinobacteria  Bacteria                                    186789
    ## 3              Apicomplexa Eukaryota                                       231
    ## 4                Aquificae  Bacteria                                      1953
    ## 5               Ascomycota Eukaryota                                      1491
    ## 6          Bacillariophyta Eukaryota                                       105
    ## 7            Bacteroidetes  Bacteria                                   1424565
    ## 8            Basidiomycota Eukaryota                                       240
    ## 9  Candidatus Poribacteria  Bacteria                                        26
    ## 10              Chlamydiae  Bacteria                                       552
    ##    ShotgunWGS-ControlPig8GutMicrobiome-Day0
    ## 1                                      3717
    ## 2                                    277130
    ## 3                                       384
    ## 4                                      2254
    ## 5                                      2196
    ## 6                                       184
    ## 7                                   1391417
    ## 8                                       405
    ## 9                                        49
    ## 10                                      829

``` r
# add phyla as rownames
rownames(phyla_counts_foraldex) <- phyla_counts_foraldex$phylum

# remove phylum, domain as columns for cleaner data
phyla_counts_foraldex <- phyla_counts_foraldex %>%
  select(-phylum, -domain)
```

### Subset by Time

#### Day 0

``` r
# subset day 0 only
Day0.Counts.Phyla.filt <- phyla_counts_foraldex %>% 
  select(ends_with("Day0"))
```

ALDEx2 function needs a factor of variables

``` r
# order alphabetically so making the meta data vector is easier
Day0.Counts.Phyla.filt <- Day0.Counts.Phyla.filt[order(colnames(Day0.Counts.Phyla.filt))]

Diets.Day0.Phyla <- as.vector(c(rep("Control", times=10), rep("Tomato", times=10)))

# check and make sure it came out right
Diets.Day0.Phyla
```

    ##  [1] "Control" "Control" "Control" "Control" "Control" "Control" "Control"
    ##  [8] "Control" "Control" "Control" "Tomato"  "Tomato"  "Tomato"  "Tomato" 
    ## [15] "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"

Run t-tests

``` r
# runs very slowly
# set cache = TRUE to save results
filt.Phyla.Day0.ByDiet.aldex <- aldex(Day0.Counts.Phyla.filt, 
                                       Diets.Day0.Phyla, 
                                       mc.samples = 1000, 
                                       test = "t", 
                                       effect = TRUE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.ttest: doing t-test

    ## aldex.effect: calculating effect sizes

``` r
filt.Phyla.Day0.ByDiet.aldex <- 
  filt.Phyla.Day0.ByDiet.aldex[order(filt.Phyla.Day0.ByDiet.aldex$we.eBH, 
                                      decreasing = FALSE),]

kable(head(filt.Phyla.Day0.ByDiet.aldex))
```

|                                       |    rab.all | rab.win.Control | rab.win.Tomato |   diff.btw |  diff.win |     effect |   overlap |     we.ep |    we.eBH |     wi.ep |    wi.eBH |
|:--------------------------------------|-----------:|----------------:|---------------:|-----------:|----------:|-----------:|----------:|----------:|----------:|----------:|----------:|
| Proteobacteria                        |  6.5857669 |       6.3145651 |      6.7528732 |  0.3746102 | 0.4414428 |  0.8688289 | 0.1646000 | 0.0083533 | 0.3221071 | 0.0110988 | 0.4064564 |
| unclassified (derived from Eukaryota) |  0.4397911 |       0.9033182 |      0.3721023 | -0.4816460 | 0.3775456 | -0.8855126 | 0.2372000 | 0.0259926 | 0.4121663 | 0.0520324 | 0.5629625 |
| Chlorophyta                           | -1.5919906 |      -1.6628591 |     -1.5087876 |  0.1732418 | 0.3127644 |  0.5175202 | 0.2666000 | 0.0670792 | 0.5280984 | 0.1077618 | 0.6142581 |
| Tenericutes                           |  0.9391923 |       1.0432714 |      0.8885737 | -0.1583287 | 0.2037393 | -0.6764046 | 0.2609478 | 0.0687257 | 0.5476915 | 0.0862333 | 0.6132091 |
| Basidiomycota                         | -2.2166312 |      -2.1306207 |     -2.3137243 | -0.1761152 | 0.2950344 | -0.5693332 | 0.2640000 | 0.0993361 | 0.5700791 | 0.1116548 | 0.6113181 |
| Elusimicrobia                         | -0.7859985 |      -0.8511297 |     -0.7338424 |  0.1259770 | 0.2378991 |  0.4818211 | 0.2789442 | 0.0942131 | 0.5839636 | 0.1284662 | 0.6345670 |

No significantly different phyla

``` r
hist(filt.Phyla.Day0.ByDiet.aldex$we.eBH,
     breaks = 45,
     main = "Histogram of p-values on the effect of diet at day 0 on phyla",
     xlab = "Benjamini Hochberg corrected p-value (we.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-206-1.png)<!-- -->

`we.eBH` is the Benjamini-Hochberg corrected p-value, and nothing is
&lt; 0.05

#### Day 7

``` r
# subset day 7 only
Day7.Counts.Phyla.filt <- phyla_counts_foraldex %>% 
  select(ends_with("Day7"))
```

ALDEx2 function needs a factor of variables

``` r
# order alphabetically so making the meta data vector is easier
Day7.Counts.Phyla.filt <- Day7.Counts.Phyla.filt[order(colnames(Day7.Counts.Phyla.filt))]

Diets.Day7.Phyla <- as.vector(c(rep("Control", times=10), rep("Tomato", times=10)))

# check and make sure it came out right
Diets.Day7.Phyla
```

    ##  [1] "Control" "Control" "Control" "Control" "Control" "Control" "Control"
    ##  [8] "Control" "Control" "Control" "Tomato"  "Tomato"  "Tomato"  "Tomato" 
    ## [15] "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"

Run t-tests

``` r
# runs very slowly
# set cache = TRUE to save results
filt.Phyla.Day7.ByDiet.aldex <- aldex(Day7.Counts.Phyla.filt, 
                                       Diets.Day7.Phyla, 
                                       mc.samples = 1000, 
                                       test = "t", 
                                       effect = TRUE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.ttest: doing t-test

    ## aldex.effect: calculating effect sizes

``` r
filt.Phyla.Day7.ByDiet.aldex <- 
  filt.Phyla.Day7.ByDiet.aldex[order(filt.Phyla.Day7.ByDiet.aldex$we.eBH, 
                                      decreasing = FALSE),]

kable(head(filt.Phyla.Day7.ByDiet.aldex))
```

|                                      |    rab.all | rab.win.Control | rab.win.Tomato |   diff.btw |  diff.win |     effect |   overlap |     we.ep |    we.eBH |     wi.ep |    wi.eBH |
|:-------------------------------------|-----------:|----------------:|---------------:|-----------:|----------:|-----------:|----------:|----------:|----------:|----------:|----------:|
| unclassified (derived from Bacteria) |  0.6375097 |       0.3089784 |      1.0304263 |  0.6599291 | 0.3822769 |  1.5816905 | 0.0495901 | 0.0000614 | 0.0027560 | 0.0002156 | 0.0095993 |
| Chrysiogenetes                       | -2.2763611 |      -2.1569750 |     -2.4322655 | -0.2743199 | 0.3633857 | -0.7063056 | 0.2112000 | 0.0499035 | 0.3550747 | 0.0541899 | 0.3126401 |
| Ascomycota                           |  0.2293522 |       0.1351592 |      0.3842469 |  0.2343789 | 0.4023438 |  0.5322208 | 0.2358000 | 0.0434242 | 0.4152844 | 0.0628626 | 0.3649693 |
| Firmicutes                           | 10.4275493 |      10.6466196 |     10.2738925 | -0.3318916 | 0.4093850 | -0.6862980 | 0.2097580 | 0.0410232 | 0.4212510 | 0.0286925 | 0.2980693 |
| Basidiomycota                        | -2.1959373 |      -2.3467395 |     -2.0850216 |  0.2574169 | 0.4261325 |  0.5450211 | 0.2453509 | 0.0808754 | 0.4670333 | 0.0853502 | 0.3993034 |
| Bacillariophyta                      | -3.3771542 |      -3.5139314 |     -3.2717009 |  0.2555352 | 0.4569006 |  0.5304439 | 0.2777445 | 0.1266018 | 0.5031681 | 0.1621232 | 0.4897163 |

One phyla was significantly different by diet at day 7 - unclassified
(derived from bacteria), padj = 0.002

``` r
hist(filt.Phyla.Day7.ByDiet.aldex$we.eBH,
     breaks = 45,
     main = "Histogram of p-values on the effect of diet at day 7 on phyla",
     xlab = "Benjamini Hochberg corrected p-value (we.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-211-1.png)<!-- -->

#### Day 14

``` r
# subset day 14 only
Day14.Counts.Phyla.filt <- phyla_counts_foraldex %>% 
  select(ends_with("Day14"))
```

ALDEx2 function needs a factor of variables

``` r
# order alphabetically so making the meta data vector is easier
Day14.Counts.Phyla.filt <- Day14.Counts.Phyla.filt[order(colnames(Day14.Counts.Phyla.filt))]

Diets.Day14.Phyla <- as.vector(c(rep("Control", times=10), rep("Tomato", times=10)))

# check and make sure it came out right
Diets.Day14.Phyla
```

    ##  [1] "Control" "Control" "Control" "Control" "Control" "Control" "Control"
    ##  [8] "Control" "Control" "Control" "Tomato"  "Tomato"  "Tomato"  "Tomato" 
    ## [15] "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"  "Tomato"

Run t-tests

``` r
filt.Phyla.Day14.ByDiet.aldex <- aldex(Day14.Counts.Phyla.filt, 
                                       Diets.Day14.Phyla, 
                                       mc.samples = 1000, 
                                       test = "t", 
                                       effect = TRUE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.ttest: doing t-test

    ## aldex.effect: calculating effect sizes

``` r
filt.Phyla.Day14.ByDiet.aldex <- 
  filt.Phyla.Day14.ByDiet.aldex[order(filt.Phyla.Day14.ByDiet.aldex$we.eBH, 
                                      decreasing = FALSE),]

kable(head(filt.Phyla.Day14.ByDiet.aldex))
```

|                                      |    rab.all | rab.win.Control | rab.win.Tomato |   diff.btw |  diff.win |     effect |   overlap |     we.ep |    we.eBH |     wi.ep |    wi.eBH |
|:-------------------------------------|-----------:|----------------:|---------------:|-----------:|----------:|-----------:|----------:|----------:|----------:|----------:|----------:|
| unclassified (derived from Bacteria) |  0.3878007 |       0.0656786 |      0.9924115 |  0.9148040 | 0.3457392 |  2.5428325 | 0.0000140 | 0.0000029 | 0.0001302 | 0.0000108 | 0.0004467 |
| Nematoda                             | -2.6749909 |      -2.9590197 |     -2.2600299 |  0.7407293 | 0.5273953 |  1.3785241 | 0.0388000 | 0.0006183 | 0.0083108 | 0.0003181 | 0.0043752 |
| Apicomplexa                          | -2.0607631 |      -2.4278818 |     -1.5566779 |  0.8009668 | 0.6160238 |  1.2318439 | 0.0901820 | 0.0011385 | 0.0131840 | 0.0015406 | 0.0147930 |
| Deinococcus-Thermus                  |  0.8934331 |       0.9663785 |      0.7934825 | -0.1756483 | 0.1458435 | -1.1470314 | 0.0980000 | 0.0062779 | 0.0330378 | 0.0047251 | 0.0286456 |
| Proteobacteria                       |  6.4735169 |       6.4218100 |      6.6020117 |  0.1992688 | 0.2071126 |  0.9684222 | 0.1323735 | 0.0064963 | 0.0402455 | 0.0063370 | 0.0366220 |
| Firmicutes                           | 10.6381375 |      10.6954982 |     10.4067027 | -0.3137346 | 0.3339467 | -0.8710526 | 0.1691662 | 0.0102429 | 0.0586362 | 0.0130723 | 0.0637137 |

``` r
hist(filt.Phyla.Day14.ByDiet.aldex$we.eBH,
     breaks = 45,
     main = "Histogram of p-values on the effect of diet at day 14 on phyla",
     xlab = "Benjamini Hochberg corrected p-value (we.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-216-1.png)<!-- -->

How many significant phyla are there?

``` r
filt.Phyla.Day14.ByDiet.aldex.sig <- 
  filt.Phyla.Day14.ByDiet.aldex[which(filt.Phyla.Day14.ByDiet.aldex$we.eBH<0.05),]

length(rownames(filt.Phyla.Day14.ByDiet.aldex.sig))
```

    ## [1] 5

5 sig phyla

Which phyla are they?

``` r
sig_day14_phyla_aldex2 <- as.data.frame(cbind(rownames(filt.Phyla.Day14.ByDiet.aldex.sig),
                                 filt.Phyla.Day14.ByDiet.aldex.sig$we.eBH))

sig_day14_phyla_aldex2
```

    ##                                     V1                   V2
    ## 1 unclassified (derived from Bacteria) 0.000130184199239088
    ## 2                             Nematoda  0.00831079482853596
    ## 3                          Apicomplexa   0.0131839873014319
    ## 4                  Deinococcus-Thermus   0.0330378200022153
    ## 5                       Proteobacteria   0.0402455424599716

-   unclassified (derived from Bacteria)
-   Nematoda
-   Apicomplexa
-   Deinococcus-Thermus
-   Proteobacteria

What is the directionality of the change?

``` r
filt.Phyla.Day14.ByDiet.aldex %>%
  select(rab.win.Control, rab.win.Tomato, we.eBH) %>%
  filter(we.eBH <= 0.05)
```

    ##                                      rab.win.Control rab.win.Tomato
    ## unclassified (derived from Bacteria)      0.06567863      0.9924115
    ## Nematoda                                 -2.95901972     -2.2600299
    ## Apicomplexa                              -2.42788185     -1.5566779
    ## Deinococcus-Thermus                       0.96637845      0.7934825
    ## Proteobacteria                            6.42180997      6.6020117
    ##                                            we.eBH
    ## unclassified (derived from Bacteria) 0.0001301842
    ## Nematoda                             0.0083107948
    ## Apicomplexa                          0.0131839873
    ## Deinococcus-Thermus                  0.0330378200
    ## Proteobacteria                       0.0402455425

Higher in control:  
\* Deinococcus-Thermus

Higher in tomato:  
\* unclassified (derived from Bacteria)  
\* Nematoda  
\* Apicomplexa  
\* Proteobacteria

### Subset by diet

#### Control

``` r
# subset control only samples across all time points, should be n=30
Control.Counts.Phyla.filt <- phyla_counts_foraldex %>% 
  select(contains("Control"))
```

ALDEx2 function needs a factor of variables

``` r
# results in pigs at different time points being grouped together
Control.Counts.Phyla.filt <- Control.Counts.Phyla.filt[order(colnames(Control.Counts.Phyla.filt))]

# then time point by "alphabetical" where 14 comes before 7
# ex, first few are Pig 10 Day 0, Pig 10 Day 14, Pig 10 Day 7, Pig 1 Day 0, Pig 1 Day 14, etc

TimePoints.Control.Phyla <- as.vector(rep(c("Day0", "Day14", "Day7"), times=10))

# check and make sure it looks right
TimePoints.Control.Phyla
```

    ##  [1] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [10] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [19] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [28] "Day0"  "Day14" "Day7"

More than two conditions this time, use the ANOVA-like test

``` r
filt.Phyla.Control.ByTime.aldex <- aldex(Control.Counts.Phyla.filt, 
                                          TimePoints.Control.Phyla, 
                                          mc.samples = 1000, 
                                          test = "kw", 
                                          effect = FALSE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.glm: doing Kruskal-Wallace and glm test (ANOVA-like)

    ## operating in serial mode

We are looking at `glm.eBH` for the BH corrected ANOVA pval

``` r
filt.Phyla.Control.ByTime.aldex <- 
  filt.Phyla.Control.ByTime.aldex[order(filt.Phyla.Control.ByTime.aldex$glm.eBH, 
                                         decreasing = FALSE),]

kable(head(filt.Phyla.Control.ByTime.aldex))
```

|                                       |     kw.ep |    kw.eBH |    glm.ep |   glm.eBH |
|:--------------------------------------|----------:|----------:|----------:|----------:|
| unclassified (derived from Bacteria)  | 0.0106850 | 0.1837993 | 0.0026968 | 0.0769864 |
| unclassified (derived from Eukaryota) | 0.0505786 | 0.2504227 | 0.0079943 | 0.1208772 |
| Dictyoglomi                           | 0.0266105 | 0.2116903 | 0.0224764 | 0.1449505 |
| Verrucomicrobia                       | 0.0385119 | 0.2325773 | 0.0257789 | 0.1537795 |
| Firmicutes                            | 0.0340934 | 0.2268098 | 0.0239987 | 0.1538562 |
| Fibrobacteres                         | 0.2588618 | 0.5236352 | 0.0316510 | 0.1592722 |

``` r
hist(filt.Phyla.Control.ByTime.aldex$glm.eBH,
     breaks = 45,
     main = "Histogram of p-values on the effect of time on control pigs on phyla",
     xlab = "Benjamini Hochberg corrected p-value (glm.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-224-1.png)<!-- -->

How many significant phyla are there?

``` r
filt.Phyla.Control.ByTime.aldex.sig <- 
  filt.Phyla.Control.ByTime.aldex[which(filt.Phyla.Control.ByTime.aldex$glm.eBH<0.05),]

length(rownames(filt.Phyla.Control.ByTime.aldex.sig))
```

    ## [1] 0

0 sig phyla

#### Tomato

``` r
# subset tomato only samples across all time points, should be n=30
Tomato.Counts.Phyla.filt <- phyla_counts_foraldex %>% 
  select(contains("Tomato"))
```

ALDEx2 function needs a factor of variables

``` r
# results in pigs at different time points being grouped together
Tomato.Counts.Phyla.filt <- Tomato.Counts.Phyla.filt[order(colnames(Tomato.Counts.Phyla.filt))]

# then time point by "alphabetical" where 14 comes before 7
# ex, first few are Pig 10 Day 0, Pig 10 Day 14, Pig 10 Day 7, Pig 1 Day 0, Pig 1 Day 14, etc

TimePoints.Tomato.Phyla <- as.vector(rep(c("Day0", "Day14", "Day7"), times=10))

# check and make sure it looks right
TimePoints.Tomato.Phyla
```

    ##  [1] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [10] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [19] "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7"  "Day0"  "Day14" "Day7" 
    ## [28] "Day0"  "Day14" "Day7"

More than two conditions this time, use the ANOVA-like test

``` r
filt.Phyla.Tomato.ByTime.aldex <- aldex(Tomato.Counts.Phyla.filt, 
                                          TimePoints.Tomato.Phyla, 
                                          mc.samples = 1000, 
                                          test = "kw", 
                                          effect = FALSE)
```

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## operating in serial mode

    ## computing center with all features

    ## aldex.glm: doing Kruskal-Wallace and glm test (ANOVA-like)

    ## operating in serial mode

We are looking at glm.eBH for the BH corrected ANOVA pval

``` r
filt.Phyla.Tomato.ByTime.aldex <- 
  filt.Phyla.Tomato.ByTime.aldex[order(filt.Phyla.Tomato.ByTime.aldex$glm.eBH, 
                                         decreasing = FALSE),]

kable(head(filt.Phyla.Tomato.ByTime.aldex))
```

|                                      |     kw.ep |    kw.eBH |    glm.ep |   glm.eBH |
|:-------------------------------------|----------:|----------:|----------:|----------:|
| unclassified (derived from Bacteria) | 0.0015220 | 0.0593092 | 0.0000096 | 0.0004317 |
| Tenericutes                          | 0.0051471 | 0.1075460 | 0.0061510 | 0.1147076 |
| Aquificae                            | 0.0368299 | 0.3521759 | 0.0264341 | 0.2773093 |
| Hemichordata                         | 0.1795195 | 0.5739416 | 0.1122691 | 0.4800990 |
| Basidiomycota                        | 0.1334207 | 0.5773943 | 0.1081339 | 0.5509906 |
| Ascomycota                           | 0.1932876 | 0.6842660 | 0.0918677 | 0.5952485 |

``` r
hist(filt.Phyla.Tomato.ByTime.aldex$glm.eBH,
     breaks = 45,
     main = "Histogram of p-values on the effect of time on tomato pigs on phyla",
     xlab = "Benjamini Hochberg corrected p-value (glm.eBH)")
```

![](Goggans_TomatoPigMicrobiomeAnalysis_GithubDoc_021921_EmmaEdits_files/figure-gfm/unnamed-chunk-230-1.png)<!-- -->

How many significant phyla are there?

``` r
filt.Phyla.Tomato.ByTime.aldex.sig <- 
  filt.Phyla.Tomato.ByTime.aldex[which(filt.Phyla.Tomato.ByTime.aldex$glm.eBH<0.05),]

length(rownames(filt.Phyla.Tomato.ByTime.aldex.sig))
```

    ## [1] 1

1 sig phyla, unclassified (derived from Bacteria)
