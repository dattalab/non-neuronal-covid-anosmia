# non-neuronal-covid-anosmia
==============================
# Non-neuronal expression of SARS-CoV-2 entry genes in the olfactory system suggests mechanisms underlying COVID-19-associated anosmia

## Authors
David H Brann\*, Tatsuya Tsukahara\*, Caleb Weinreb\*, Marcela Lipovsek, Koen Van den Berge, Boying Gong, Rebecca Chance, Iain C Macaulay, Hsin-Jung Chou, Russell B Fletcher, Diya Das, Kelly Street, Hector Roux de Bezieux, Yoon-Gi Choi, Davide Risso, Sandrine Dudoit, Elizabeth Purdom, Jonathan Mill, Ralph Abi Hachem, Hiroaki Matsunami, Darren W Logan, Bradley J Goldstein, Matthew S Grubb, John Ngai, Sandeep Robert Datta†

\*These authors contributed equally to this work.

†Corresponding author. Email: srdatta@hms.harvard.edu

# Abstract
Altered olfactory function is a common symptom of COVID-19, but its etiology is unknown. A key question is whether SARS-CoV-2 (CoV-2) – the causal agent in COVID-19 – affects olfaction directly, by infecting olfactory sensory neurons or their targets in the olfactory bulb, or indirectly, through perturbation of supporting cells. Here we identify cell types in the olfactory epithelium and olfactory bulb that express SARS-CoV-2 cell entry molecules. Bulk sequencing demonstrated that mouse, non-human primate and human olfactory mucosa expresses two key genes involved in CoV-2 entry, ACE2 and TMPRSS2. However, single cell sequencing revealed that ACE2 is expressed in support cells, stem cells, and perivascular cells, rather than in neurons. Immunostaining confirmed these results and revealed pervasive expression of ACE2 protein in dorsally-located olfactory epithelial sustentacular cells and olfactory bulb pericytes in the mouse. These findings suggest that CoV-2 infection of non-neuronal cell types leads to anosmia and related disturbances in odor perception in COVID-19 patients.

## Manuscript

For more details, please see our Open Access manuscript in Science Advances [here](https://doi.org/10.1126/sciadv.abc5801).

# Installation

## Requirements
1. Make a new conda env, e.g. `conda create -n cov python=3.8`
2. Activate that env `conda activate cov`.
3. Clone and enter this repo: `git clone git@github.com:dattalab/non-neuronal-covid-anosmia.git && cd non-neuronal-covid-anosmia`
4. To install the specific versions of packages used when testing the scripts in this repo you can do `pip install -r requirements.txt`. Alternatively, the minimal requirements for running the scripts and notebooks in this repo are:
```
pip install scanpy loompy
pip install jupyter notebook
```
5. Install the code in this directory from the `setup.py` file via `pip install -e .`

## Data

1. Raw and processed datasets generated as part of this study are available from the NCBI GEO at accessions: [GSE151346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151346) (Mouse olfactory epithelium), [GSE153730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153730) (Mouse olfactory epithelium injury), [GSE148360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148360) (Mouse olfactory bulb), and [GSE151709](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151709) (Mouse OB DA neurons). 
2. As described in the notebooks, the notebooks in this repository require the respective raw counts files for each dataset to be downloaded to the [data/raw](./data/raw) folders. These data can be found at the above GEO links, or at [GSE139522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522) for the reanalysis of the human olfactory epithelium data from [Durante et al. 2020](https://doi.org/10.1038/s41593-020-0587-9).
3. Additional metadata is saved in the [data/tables](./data/tables) folder.

# Examples
Code to replicate analyses in [Brann, Tsukahara, Weinreb et al. 2020](https://doi.org/10.1126/sciadv.abc5801). 

1. Open a new jupyter notebook with `jupyter notebook`.
2. Run the [notebooks](./notebooks). These notebooks should load the raw counts, combine them with the respective metadata and cell type annotations, and generate the main plots in the manuscript.

# Contact
For more details on additional preprocessing steps or analyses, please consult the methods in our manuscript, post an issue here, or contact the authors.