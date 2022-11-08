# tomato-pig-microbiome
### Whole genome sequencing analysis of gut microbiome of pigs consuming control or tomato containing diets.

Title: Short term tomato consumption alters the pig gut microbiome towards a more favorable profile

Running title: Tomato consumption alters the pig gut microbiome

Mallory L. Goggans1, Emma A. Bilbrey2, Cristian D. Quiroz-Moreno2, David M. Francis3, Sheila K. Jacobi4, Jasna Kovac5,6, Jessica L. Cooperstone1,2,*

1Food Science and Technology, The Ohio State University, Columbus, OH 43210.
2Horticulture and Crop Science, The Ohio State University, Columbus, OH 43210.
3Horticulture and Crop Science, The Ohio State University, Wooster, OH 44691.
4Animal Sciences, The Ohio State University, Columbus, OH 43210.
5Food Science, The Pennsylvania State University, University Park, PA, 16802.
6 Microbiome Center, Huck Institutes of the Life Sciences, The Pennsylvania State University, University Park, PA 16802.
*to whom correspondence should be addressed, [cooperstone.1@osu.edu](cooperstone.1@osu.edu) 

ABSTRACT
Diets rich in fruits and vegetables have been shown to exert positive effects on the gut microbiome. However, little is known about the specific effect of individual fruits or vegetables on gut microbe profiles. This study aims to elucidate the effects of tomato consumption on the gut microbiome, as tomatoes account for 22% of vegetable consumption in Western diets, and their consumption has been associated with positive health outcomes. Using piglets as a physiologically relevant model of human metabolism, 20 animals were assigned to either a control or a tomato powder-supplemented diet (both macronutrient matched and isocaloric) for 14 days. The microbiome was sampled rectally at three time points: day 0 (baseline), day 7 (midpoint), and day 14 (end of study). DNA was sequenced using shotgun metagenomics, and reads were annotated using MG-RAST. There were no differences in body weight or feed intake between our two treatment groups. There was a microbial shift which included a higher ratio of Bacteroidota to Bacillota (formerly known as Bacteroidetes and Firmicutes, respectively) and higher alpha-diversity in tomato-fed animals, indicating a shift to a more desirable phenotype. Analyses at both the phylum and genus levels showed global microbiome profile changes (permutational multivariate analysis of variance [PERMANOVA], P ≤ 0.05) over time but not with tomato consumption. These data suggest that short-term tomato consumption can beneficially influence the gut microbial profile, warranting further investigation in humans.

Full work published in Goggans et al., Microbiology Spectrum, 2022, [https://doi.org/10.1128/spectrum.02506-22](https://doi.org/10.1128/spectrum.02506-22) 

Raw data files annotated using MG-RAST.

MG-RAST project link: https://www.mg-rast.org/linkin.cgi?project=mgp93233

NCBI SRA fastq files: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA601162/

Data analysis completed using both phylum and genus levels of taxonomic identification.

R Markdown and Github doc containing all code are included in this repository. Input data is in the [supplement](https://journals.asm.org/doi/suppl/10.1128/spectrum.02506-22/suppl_file/spectrum.02506-22-s0001.xlsx) to our [paper](https://doi.org/10.1128/spectrum.02506-22)/.

GENERA:

- Data filtering
- Microbiome profile
- Rarefaction Curves
- Krona plots
- PERMANOVA
- PCoA Beta Diversity
- Alpha diversity
- Compositional analysis using ALDEx2


PHYLA:

- Data filtering
- Microbiome Profile
- PERMANOVA
- PCoA Beta Diversity
- Analysis of Bacteroidetes, Firmicutes, and ratio of B to F
- Alpha diversity
- Compositional analysis using ALDEx2



