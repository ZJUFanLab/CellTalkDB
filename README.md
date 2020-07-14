# CellTalkDB

### A manually curated database of ligand-receptor interactions in human and mouse

<img src='https://img.shields.io/badge/ligand--receptor-database-brightgreen'> <img src='https://img.shields.io/badge/human-3%2C398-orange'> <img src='https://img.shields.io/badge/mouse-2%2C033-blue'> 

Cell-cell communications via secreting and receiving ligands frequently occur in multicellular organisms, which is a vital feature involving numerous biological processes.[[1]](https://pubmed.ncbi.nlm.nih.gov/32435978/) Recent advancements in single-cell RNA sequencing (scRNA-seq) have effectively resolving cellular phenotype heterogeneity and cell-type composition of complex tissues, which enables systematic investigation of cell-cell communications at a single-cell resolution. However, the common practice to study chemical signal-dependent cell-cell communications with scRNA-seq relies heavily on the prior knowledge of ligand-receptor (LR) interaction pairs. Here, we introduce CellTalkDB, a comprehensive database of LR interaction pairs in human and mouse by text mining and manual verification of known protein-protein interactions (PPIs) in STRING.

|Species  |LR pairs|Ligands|Receptors|Evidences|
|:---:    |:---:   |:---:  | :---:   | :---:   |
|__Human__| __3,398__  |815    |780      |3,735    |
|__Mouse__| __2,033__  |651     |588     |2,601    |

[1] Xin Shao et al. New Avenues for Systematically Inferring Cell-Cell Communication: Through Single-Cell Transcriptomics Data.Protein Cell.2020 May 21.[doi: 10.1007/s13238-020-00727-5](https://link.springer.com/article/10.1007/s13238-020-00727-5)

# Workflow
<img src='https://github.com/ZJUFanLab/CellTalkDB/blob/master/img/curation.svg'>

# Download
[![R data](https://img.shields.io/badge/R-data-blueviolet)](https://github.com/ZJUFanLab/CellTalkDB/tree/master/database) [![txt data](https://img.shields.io/badge/txt-data-ff69b4.svg)](http://tcm.zju.edu.cn/celltalkdb) [![Ensembl v99](https://img.shields.io/badge/Ensembl-v99-brightgreen)](http://www.ensembl.org) [![UniProt](https://img.shields.io/badge/UniProt-2020__03-yellowgreen)](https://www.uniprot.org/) [![NCBI](https://img.shields.io/badge/NCBI-2020--04--28-orange)](https://www.ncbi.nlm.nih.gov/)

- LR pairs for human and mouse can be download in[`database/`](https://github.com/ZJUFanLab/CellTalkDB/tree/master/database) 
- Annotation data for LR pairs can be downloaded in[`data/`](https://github.com/ZJUFanLab/CellTalkDB/tree/master/data)
- Raw data for reproduction of our results can be downloaded in the [release](https://github.com/ZJUFanLab/CellTalkDB/releases) page.

# Usage
Users can download the LR pairs in CellTalkDB and replace the underlying database in [SoptSC](https://github.com/mkarikom/RSoptSC), [SingleCellSignalR](https://github.com/SCA-IRCM/SingleCellSignalR_v1) and [CellPhoneDB](https://github.com/Teichlab/cellphonedb), etc. to identify significantly enriched LR pairs and to infer cell-cell communications. 

To help users use CellTalkDB conveniently, we have developed an integrated R packge by only replacing the underlying database in SingleCellSignalR and keeping the other functions unchanged. Source package of `scsrctdb-1.0` can be downloaded in the [release](https://github.com/ZJUFanLab/CellTalkDB/releases) page.

### Install
[![R >3.6](https://img.shields.io/badge/R-%3E3.6-brightgreen)](https://github.com/ZJUFanLab/CellTalkDB/releases) [![source package scsrctdb-1.0.tar.gz](https://img.shields.io/badge/source%20package-scsrctdb--1.0.tar.gz-blue)](https://github.com/ZJUFanLab/CellTalkDB/releases)

```
# download the source package of scsrctdb-1.0.tar.gz and install it
# ensure the right directory for scsrctdb-1.0.tar.gz
install.packages(pkgs = 'scsrctdb-1.0.tar.gz',repos = NULL, type = "source")
```

### Examples

- For human scRNA-seq datasets
```
# based on 3,398 human LR pairs
library(scsrctdb)
cell_signal <- cell_signaling(data = data,
                              genes = genes,
                              cluster = cluster,
                              gene_resive = T,
                              species = 'homo sapiens')
```

- For mouse scRNA-seq datasets
```
# based on 2,033 mouse LR pairs
library(scsrctdb)
cell_signal <- cell_signaling(data = data,
                              genes = genes,
                              cluster = cluster,
                              gene_resive = T,
                              species = 'mus musculus')
```

- Use `visualize()` for plotting
```
visualize(cell_signal)
```
<img src='https://github.com/ZJUFanLab/CellTalkDB/blob/master/img/Fig1.png'>

```
visualize(cell_signal,show.in = 1)
```
<img src='https://github.com/ZJUFanLab/CellTalkDB/blob/master/img/Fig2.png'>

__Note__: we have added an extra parameter`gene_resive`to revise gene symbols according to [NCBI Gene database](https://www.ncbi.nlm.nih.gov/gene) (updated in April 28,2020) as CellTalkDB has been revised with it. For more information about how to use SingleCellSignalR, please refer to [wiki page](https://github.com/ZJUFanLab/CellTalkDB/wiki/SCSR_UsersGuide) or  [SCA-IRCM/SingleCellSignalR](https://github.com/SCA-IRCM/SingleCellSignalR_v1)




# About

CellTalkDB manuscript is under submitting. For more information, please visit our website [tcm.zju.edu.cn/celltalkdb](http://tcm.zju.edu.cn/celltalkdb/). Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn
