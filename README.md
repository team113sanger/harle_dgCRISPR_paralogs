# A compendium of synthetic lethal gene pairs defined by extensive combinatorial pan-cancer CRISPR screening

[![DOI](https://zenodo.org/badge/628953602.svg)](https://doi.org/10.5281/zenodo.15827200)


## Installing dependencies 

Please see [INSTALL_README.md](INSTALL_README.md) for instructions on installing dependencies.

## CRISPR analysis

To download the datasets and metadata used in the analysis of the combinatorial CRISPR screen.

1. Change into the directory which already contains the Jupyter Notebooks and scripts that were used for the data analysis.

```
cd combinatorial_crispr_screen_analysis
```

2. Download the datasets from [Figshare](https://doi.org/10.6084/m9.figshare.25954027.v4).

```
curl -k -o dnld.zip -O https://figshare.com/ndownloader/articles/25954027/versions/4?folder_path=dnld
unzip -j dnld.zip 'DATA.tar.gz' 'METADATA.tar.gz' && rm dnld.zip
find . -name '*.tar.gz' -exec sh -c 'tar -xzvf "$1" -C "$(dirname "$1")" && rm "$1"' _ {} \;
find DATA -name '*.tar.gz' -exec sh -c 'tar -xzvf "$1" -C "$(dirname "$1")" && rm "$1"' _ {} \;
```

Initial data required for the analysis of the CRISPR screen can be found in the `DATA` and `METADATA` directories. Python Notebooks which detail the analyses can be found in `JupyterNotebooks`. These rely on R and Bash scripts in the `SCRIPTS` directory.

Please see [combinatorial_crispr_screen_analysis/README.md](combinatorial_crispr_screen_analysis/README.md).

## Imaging analysis

To download the datasets used in the analysis of the imaging dataset.

1. Change into the relevant directory.

```
cd imaging_analysis
```

2. Download the datasets from [Figshare](https://doi.org/10.6084/m9.figshare.25954027.v4).

```
curl -k -o dnld.zip -O https://figshare.com/ndownloader/articles/25954027/versions/4?folder_path=imaging_data
unzip -j dnld.zip 'DATA.tar.gz' && rm dnld.zip
find . -name '*.tar.gz' -exec sh -c 'tar -xzvf "$1" -C "$(dirname "$1")" && rm "$1"' _ {} \;
find DATA -name '*.tar.gz' -exec sh -c 'tar -xzvf "$1" -C "$(dirname "$1")" && rm "$1"' _ {} \;
```

## Manuscript figures

Please see [MANUSCRIPT/README.md](MANUSCRIPT/README.md).

