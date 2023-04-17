# Installation

Details of the analysis, figure and table generation are set out in Jupyter Notebooks. To run the scripts within those Notebooks, it is necessary to either install the following software and 
dependencies locally or we have provided a Dockerfile and container.

## Software dependencies

### R

Analyses were run using [R version 4.2.3](https://cran.r-project.org/) with the following package dependencies:
 
* `curl`
* `optparse`
* `stringi`
* `tidyverse`
* `reshape2`
* `vroom`
* `data.table`
* `ggpubr`
* `ggsci`
* `ggrepel`
* `patchwork`
* `GGally`
* `ggridges`
* `pROC`

### MAGeCK

[MAGeCK version 0.5.9.3](https://sourceforge.net/projects/mageck/files/0.5/) was used in the single guide analsis to identify enriched and depleted guides and genes.

### BAGEL2

[BAGEL2 version 2.0 build 114 (commit reference f9eedca)](https://github.com/hart-lab/bagel/tree/f9eedca7dc16299943dd1fd499bc1df4350ce8ef) was used in the single guide analysis to identify depleted 
guides and genes.

## Docker

A [Dockerfile](./Dockerfile) has been provided which contains all software dependencies required for running scripts and analyses in this repository.

To build:

```
docker build -t paralog_dgcrispr .
```

To run R:

```
docker run -it paralog_dgcrispr R --help
```

To run Rscript:

```
docker run -it paralog_dgcrispr Rscript --help
```

To run MAGeCK:

```
docker run -it paralog_dgcrispr mageck -h
```

To run BAGEL2:

```
docker run -it paralog_dgcrispr BAGEL.py -h
```
