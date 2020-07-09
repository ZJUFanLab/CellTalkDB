library(RCurl)
library(stringr)

# Data processing for human
# download protein annotation from STRING v11.0 (9606.protein.info.v11.0.txt)
pinfo<- read.table(file = '9606.protein.info.v11.0.txt',
                   stringsAsFactors = F,
                   sep = '\t',
                   header = T)

pinfo$protein_external_id<- substr(pinfo$protein_external_id,start = 6,stop = 24)

# download ppi from STRING v11.0 (9606.protein.links.v11.0.txt)
human_ppi<- read.table(file = '9606.protein.links.v11.0.txt',
                       stringsAsFactors = F,
                       header = T)
# obtain 11,759,454 ppi

# download protein and gene annotation from Ensembl v99 (Homo_sapiens.GRCh38.99.entrez.tsv)
pinfo_ginfo<- read.table(file = 'Homo_sapiens.GRCh38.99.entrez.tsv',
                         header = T,
                         sep = '\t',
                         stringsAsFactors = F)
# exclude duplicated ppi
human_ppi<- human_ppi[,c("protein1","protein2")]
human_ppi<- unique(human_ppi)

# revise the Ensembl protein ID
d1<- substr(human_ppi$protein1,start = 6,stop = 24)
human_ppi$protein1<- d1
d1<- substr(human_ppi$protein2,start = 6,stop = 24)
human_ppi$protein2<- d1

# create a sort function
p_sort<- function(object){
  object<- object[order(object)]
  object<- paste(object[1],object[2],sep = '_')
  return(object)
}

# exclude directed ppi in STRING protein network data
human_ppi2<- apply(human_ppi, 1, p_sort)
human_ppi$protein<- human_ppi2
human_ppi1<- unique(human_ppi$protein)
human_ppi$protein_c<- paste(human_ppi$protein1,human_ppi$protein2,sep = '_')
d1<- which(human_ppi$protein == human_ppi$protein_c)
human_ppi<- human_ppi[d1,]
human_ppi<- human_ppi[,c("protein1","protein2")]
# obtain 5,879,727 ppi

# annotate ppi with protein and gene information in STRING and Ensembl
p_all<- c(human_ppi$protein1,human_ppi$protein2)
p_all<- unique(p_all)
p_all<- data.frame(protein_id = p_all,
                   gene_id_length = 0,
                   gene_id = 'NA',
                   db_name = 'NA',
                   EntrezGene_id = 'NA',
                   pinfo_name = 'NA',
                   stringsAsFactors = F)

for (i in 1:nrow(p_all)) {
  print(i)
  d1<- p_all$protein_id[i]
  if(d1 %in% pinfo_ginfo$protein_stable_id){
    pinfo_ginfo1<- pinfo_ginfo[pinfo_ginfo$protein_stable_id == d1,]
    pinfo_ginfo1<- pinfo_ginfo1[,c("gene_stable_id","xref","db_name")]
    pinfo_ginfo1<- unique(pinfo_ginfo1)
    p_all$gene_id_length[i]<-  nrow(pinfo_ginfo1)
    if(nrow(pinfo_ginfo1) == 1){
      p_all$gene_id[i]<- pinfo_ginfo1$gene_stable_id
      p_all$db_name[i]<- pinfo_ginfo1$db_name
      p_all$EntrezGene_id[i]<- pinfo_ginfo1$xref
    }
    if(nrow(pinfo_ginfo1) > 1){
      if(length(unique(pinfo_ginfo1$gene_stable_id)) == 1){
        p_all$gene_id[i]<- unique(pinfo_ginfo1$gene_stable_id)
      }
      if(length(unique(pinfo_ginfo1$db_name)) == 1){
        p_all$db_name[i]<- unique(pinfo_ginfo1$db_name)
      }
    }
  }
  if(d1 %in% pinfo$protein_external_id){
    pinfo1<- pinfo[pinfo$protein_external_id == d1,]
    p_all$pinfo_name[i]<- pinfo1$preferred_name
  }
}

# load processed gene annotation from NCBI Gene database(gene2ensembl,updated in 2020.04.28)
human_gene2ensembl<- readRDS(file = '../data/human_gene2ensembl.rds')
human_gene_info<- readRDS(file = '../data/human_gene_info.rds')

# annotate ppi with gene information in NCBI (gene2ensembl)
p_all1<- p_all[p_all$protein_id %in% human_gene2ensembl$Ensembl_protein,]
p_all1<- p_all1[,c("protein_id","gene_id","EntrezGene_id")]
human_gene2ensembl1<- human_gene2ensembl[human_gene2ensembl$Ensembl_protein %in% p_all1$protein_id,]
rownames(human_gene2ensembl1)<- human_gene2ensembl1$Ensembl_protein
human_gene2ensembl1<- human_gene2ensembl1[p_all1$protein_id,]
all(p_all1$protein_id == human_gene2ensembl1$Ensembl_protein)
p_all1<- human_gene2ensembl1
p_all1<- p_all1[,c(3,2,1)]
rownames(p_all1)<- 1:nrow(p_all1)

p_all2<- p_all[!p_all$protein_id %in% human_gene2ensembl$Ensembl_protein,]
p_all3<- p_all2[!(p_all2$EntrezGene_id == 'NA' & p_all2$pinfo_name == 'NA'),]
p_all2<- p_all2[!p_all2$protein_id %in% p_all3$protein_id,]

p_all4<- p_all3[p_all3$EntrezGene_id == 'NA',]
p_all3<- p_all3[!p_all3$protein_id %in% p_all4$protein_id,]
p_all5<- p_all4[!p_all4$pinfo_name %in% human_gene_info$Symbol,]
p_all4<- p_all4[!p_all4$protein_id %in% p_all5$protein_id,]

human_gene_info1<- human_gene_info[human_gene_info$Symbol %in% p_all4$pinfo_name,c("GeneID","Symbol")]
human_gene_info1<- unique(human_gene_info1)
rownames(human_gene_info1)<- human_gene_info1$Symbol
human_gene_info1<- human_gene_info1[p_all4$pinfo_name,]
p_all4$EntrezGene_id<- human_gene_info1$GeneID

# load gene annotation 
human_ann_manual<- readRDS(file = '../data/human_ann_manual.rds')
p_all2$pinfo_name<- human_ann_manual$EntrezGene_id
p_all2_1<- p_all2[!p_all2$pinfo_name %in% human_gene_info$Symbol,]
p_all2<- p_all2[!p_all2$protein_id %in% p_all2_1$protein_id,]
p_all2_1<- p_all2_1[p_all2_1$pinfo_name %in% human_gene_info$Synonyms,]
p_all2_1<- p_all2_1[p_all2_1$pinfo_name %in% c('C10orf113','KIAA1257'),]
p_all2_1[p_all2_1$pinfo_name == 'KIAA1257',]$EntrezGene_id<- '57501'
p_all2_1[p_all2_1$pinfo_name == 'C10orf113',]$EntrezGene_id<- '10529'

human_gene_info1<- human_gene_info[human_gene_info$Symbol %in% p_all2$pinfo_name,c("GeneID","Symbol")]
human_gene_info1<- unique(human_gene_info1)
human_gene_info1<- human_gene_info1[-which(human_gene_info1$GeneID == '100187828'),]
rownames(human_gene_info1)<- human_gene_info1$Symbol
human_gene_info1<- human_gene_info1[p_all2$pinfo_name,]
p_all2$EntrezGene_id<- human_gene_info1$GeneID

p_all2<- rbind(p_all2,p_all3,p_all4)

p_all2<- p_all2[,c(1,3,5)]
colnames(p_all2)<- colnames(p_all1)

p_all<- rbind(p_all1,p_all2)

# obtain 19,101 unique proteins

# match human_ppi
human_ppi<- human_ppi[human_ppi$protein1 %in% p_all$Ensembl_protein,]
human_ppi<- human_ppi[human_ppi$protein2 %in% p_all$Ensembl_protein,]
# obtain 5,796,319 ppi

# annotate
p_all$gene_symbol<- 'NA'
p_all$des<- 'NA'

for (i in 1:nrow(p_all)) {
  print(i)
  d1<- p_all[i,]
  d2<- human_gene_info[human_gene_info$GeneID == d1$Gene_id,]
  d2<- d2[,c("Symbol","description")]
  d2<- unique(d2)
  p_all[i,"gene_symbol"]<- d2$Symbol
  p_all[i,"des"]<- d2$description
}

p_all1<- p_all
rownames(p_all1)<- p_all1$Ensembl_protein
p_all1<- p_all1[human_ppi$protein1,]
human_ppi<- cbind(human_ppi,p_all1)
human_ppi<- human_ppi[,-3]
colnames(human_ppi)[3:6]<- c('pro_gene1','gene1_id','gene1_symbol','gene1_description')

p_all1<- p_all
rownames(p_all1)<- p_all1$Ensembl_protein
p_all1<- p_all1[human_ppi$protein2,]
human_ppi<- cbind(human_ppi,p_all1)
human_ppi<- human_ppi[,-7]
colnames(human_ppi)[7:10]<- c('pro_gene2','gene2_id','gene2_symbol','gene2_description')
human_ppi<- human_ppi[,c(1,2,3,7,4,8,5,9,6,10)]

# classify proteins into poteintial ligands and receptors manually
# load classified protein data
human_potential_lr<- readRDS(file = '../data/human_potential_lr.rds')

# obtain 1,959 potential ligands and 3,157 potential receptors
human_ligand<- human_potential_lr[human_potential_lr$Non_lr_gene_manual == 'ligand',]
human_ligand<- human_ligand[human_ligand$type_gene == 'protein-coding',]

human_receptor<- human_potential_lr[human_potential_lr$Non_lr_gene_manual == 'receptor',]
human_receptor<- human_receptor[human_receptor$type_gene == 'protein-coding',]

# match ligands and receptors
human_ppi1<- human_ppi[human_ppi$protein1 %in% human_ligand$Ensembl_protein & human_ppi$protein2 %in% human_receptor$Ensembl_protein,]
human_ppi2<- human_ppi[human_ppi$protein1 %in% human_receptor$Ensembl_protein & human_ppi$protein2 %in% human_ligand$Ensembl_protein,]

# remove duplicayed LR pairs
p_sort<- function(object){
  object<- object[order(object)]
  object<- paste(object[1],object[2],sep = '_')
  return(object)
}

d1<- human_ppi1[,c("gene1_symbol","gene2_symbol")]
d2<- apply(d1, 1, p_sort)
human_ppi1$com_gene<- d2

d1<- human_ppi2[,c("gene1_symbol","gene2_symbol")]
d2<- apply(d1, 1, p_sort)
human_ppi2$com_gene<- d2

rownames(human_ppi1)<- 1:nrow(human_ppi1)
d1<- human_ppi1[,c("gene1_symbol","gene2_symbol","com_gene")]
d1<- unique(d1)
human_ppi1<- human_ppi1[rownames(d1),]

rownames(human_ppi2)<- 1:nrow(human_ppi2)
d1<- human_ppi2[,c("gene1_symbol","gene2_symbol","com_gene")]
d1<- unique(d1)
human_ppi2<- human_ppi2[rownames(d1),]

human_ppi3<- human_ppi2[!human_ppi2$com_gene %in% human_ppi1$com_gene,]
human_ppi3<- human_ppi3[,c(2,1,4,3,6,5,8,7,10,9,11)]
colnames(human_ppi3)<- colnames(human_ppi1)
human_ppi1<- rbind(human_ppi1,human_ppi3)
colnames(human_ppi1)<- c('ligand','receptor',
                         'ligand_ensembl','receptor_ensembl',
                         'ligand_gene_id','receptor_gene_id',
                         'ligand_gene_symbol','receptor_gene_symbol',
                         'ligand_description','receptor_description','com_gene')
# obtain 265,203 potential LR pairs

# load uniprot protein knowledegbase
human_uniprot<- readRDS('../data/human_uniprot.rds')

# geneating keyword in searching term for API
human_ppi1$search_term<- 'NA'

for (i in 1:nrow(human_ppi1)) {
  print(i)
  gene1<- human_ppi1$ligand_gene_symbol[i]
  gene1_name<- unique(human_gene_info[human_gene_info$Symbol == gene1,]$Synonyms)
  gene1_name<- gene1_name[which(gene1_name != '-')]
  
  if (gene1 %in% human_uniprot$gene) {
    gene1_pro<- unique(human_uniprot[human_uniprot$gene == gene1,]$protein)
    gene1<- c(gene1,gene1_pro)
  }
  
  gene1<- c(gene1,gene1_name)
  gene1<- unique(gene1)
  gene1_search_API<- paste0(gene1,'%5BTitle%2FAbstract%5D')
  gene1_name<- gene1[1]
  gene1_name_search_API<- gene1_search_API[1]
  
  if (length(gene1) > 1) {
    for (j in 2:length(gene1)) {
      gene1_name<- paste(gene1_name,gene1[j],sep = ',')
      gene1_name_search_API<- paste(gene1_name_search_API,'OR',gene1_search_API[j],sep = '+')
    }
  }
  
  
  gene2<- human_ppi1$receptor_gene_symbol[i]
  gene2_name<- unique(human_gene_info[human_gene_info$Symbol == gene2,]$Synonyms)
  gene2_name<- gene2_name[which(gene2_name != '-')]
  
  if (gene2 %in% human_uniprot$gene) {
    gene2_pro<- unique(human_uniprot[human_uniprot$gene == gene2,]$protein)
    gene2<- c(gene2,gene2_pro)
  }
  
  gene2<- c(gene2,gene2_name)
  gene2<- unique(gene2)
  gene2_search_API<- paste0(gene2,'%5BTitle%2FAbstract%5D')
  gene2_name<- gene2[1]
  gene2_name_search_API<- gene2_search_API[1]
  
  if (length(gene2) > 1) {
    for (j in 2:length(gene2)) {
      gene2_name<- paste(gene2_name,gene2[j],sep = ',')
      gene2_name_search_API<- paste(gene2_name_search_API,'OR',gene2_search_API[j],sep = '+')
    }
  }
  
  gene1_name_search_API<- paste0('%28',gene1_name_search_API,'%29')
  gene2_name_search_API<- paste0('%28',gene2_name_search_API,'%29')
  
  gene_name_search_API<- paste(gene1_name_search_API,'AND',gene2_name_search_API,sep = '+')
  
  human_ppi1[i,"search_term"]<- gene_name_search_API

}


# Exclude LR pairs without matched articles with Pubmed E-utilities

# Warning: please read the rule of NCBI E-utilities usage carefully before running the codes below.
human_ppi1$count<- '-1'

for (i in 1:nrow(human_ppi1)) {
  print(i)
  d1<- human_ppi1[i,]
  d1_term<- d1$search_term
  # API key is removed
  d1_url<- paste('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=',d1_term,sep = '')
  d1_res<- getURL(url = d1_url)
  d1_1<- str_locate(string = d1_res,pattern = '<Count>')
  d1_2<- str_locate(string = d1_res,pattern = '</Count>')
  d1_res1<- str_sub(string = d1_res,start = d1_1[2]+1,end = d1_2[1]-1)
  human_ppi1$count[i]<- d1_res1
  # Sys.sleep is removed
}

# Remove LR pairs without matched articles

human_ppi1$count<- as.numeric(human_ppi1$count)
human_ppi1<- human_ppi1[human_ppi1$count > 0,]

# obtain 25,317 potential LR pairs for manual verification.
