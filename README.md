# Repository for the paper Insights into Diversity, Host-Range, and Temporal Stability of _Bacteroides_ and  _Phocaeicola_ Prophages

This repository contains scripts used to discover prophages from WGS data and R scripts used to analyze them and visualization.


## Prophage discovery part
Script are within the folder __ProphageDiscovery__.

The WGS data were retrieved from the BioProjects: PRJNA843113, PRJNA636979.


We have used [PHASTER](https://phaster.ca/), [Vibrant](https://github.com/AnantharamanLab/VIBRANT) and [Cenote Taker 3](https://github.com/mtisza1/Cenote-Taker3/tree/main).

Quality of retrieved prophages was checked by [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/).

Clustering at 100% and 95% similarity was performed using [MMseqs2](https://github.com/soedinglab/MMseqs2).

High-quality prophage genomes were annotated using [DRAM-v](https://github.com/WrightonLabCSU/DRAM/tree/master?tab=readme-ov-file).

Taxonomic classification was assigned by using [vConTACT3](https://bitbucket.org/MAVERICLab/vcontact3/src/master/), [PhageScope](https://phagescope.deepomics.org/) and BLAST.

To investigate the prevalence of identified prophages in the gut metaviromes we used [coverM](https://github.com/wwood/CoverM).


## Data Analysis
R script and all required data are within the folder __DataAnalysis__.

We used R (version 4.0.3) in RStudio (version 2022.12.0+353).




All tools have to be installed prior using any code in this repository.
