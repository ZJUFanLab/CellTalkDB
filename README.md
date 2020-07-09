# CellTalkDB

### A manually curated database of ligand-receptor interactions in human and mouse

<img src='https://img.shields.io/badge/ligand--receptor-database-brightgreen'> <img src='https://img.shields.io/badge/human-3%2C398-orange'> <img src='https://img.shields.io/badge/mouse-1%2C988-blue'> 

Cell-cell communications via secreting and receiving ligands frequently occur in multicellular organisms, which is a vital feature involving numerous biological processes.[[1]](https://pubmed.ncbi.nlm.nih.gov/32435978/) Recent advancements in single-cell RNA sequencing (scRNA-seq) have effectively resolving cellular phenotype heterogeneity and cell-type composition of complex tissues, which enables systematic investigation of cell-cell communications at a single-cell resolution. However, the common practice to study chemical signal-dependent cell-cell communications with scRNA-seq relies heavily on the prior knowledge of ligand-receptor (LR) interaction pairs. Here, we introduce CellTalkDB, a comprehensive database of LR interaction pairs in human and mouse by text mining and manual verification of known protein-protein interactions (PPIs) in STRING.

|Species  |LR pairs|Ligands|Receptors|Evidences|
|:---:    |:---:   |:---:  | :---:   | :---:   |
|__Human__| 3,398  |815    |780      |1,557    |
|__Mouse__| 1,988  |NA     |NA       |NA       |

[1] Xin Shao et al. New Avenues for Systematically Inferring Cell-Cell Communication: Through Single-Cell Transcriptomics Data.Protein Cell.2020 May 21.[doi: 10.1007/s13238-020-00727-5](https://link.springer.com/article/10.1007/s13238-020-00727-5)

# Data curation

<img src='https://github.com/ZJUFanLab/CellTalkDB/blob/master/img/curation.png'>

# Download
<img src='https://img.shields.io/badge/R-data-blueviolet'> [![txt data](https://img.shields.io/badge/txt-data-ff69b4.svg)](http://tcm.zju.edu.cn/celltalkdb)

- Raw data for reproduction of our results can be downloaded in the release page.
- LR pairs for human and mouse can be download in `database/lr_pairs.rds` 
- Annotation data for LR pairs can be downloaded in `data/` and `database/lr_pairs_ann.rds`

# About
For more information, please visit our website.Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn
