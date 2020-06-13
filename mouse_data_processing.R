library(RCurl)
library(stringr)

# Data processing for mouse
# download protein annotation from STRING v11.0 (10090.protein.info.v11.0.txt)
pinfo<- read.table('D:/workspace/Rstudio/STRING_v11/2020NAR/10090.protein.info.v11.0.txt',
                   stringsAsFactors = F,
                   sep = '\t',
                   header = T)

pinfo$protein_external_id<- substr(pinfo$protein_external_id,start = 7,stop = 24)

# download ppi from STRING v11.0 (10090.protein.links.v11.0.txt)
mouse_ppi<- read.table('D:/workspace/Rstudio/STRING_v11/2020NAR/10090.protein.links.v11.0.txt',
                 stringsAsFactors = F,
                 header = T)
# obtain 11,944,806 ppi

# download protein and gene annotation from Ensembl v99 (Mus_musculus.GRCm38.99.entrez.tsv)
pinfo_ginfo<- read.table('D:/workspace/Rstudio/STRING_v11/Mus_musculus.GRCm38.99.entrez.tsv',
                         header = T,
                         sep = '\t',
                         stringsAsFactors = F)

# exclude duplicated ppi
mouse_ppi<- mouse_ppi[,c("protein1","protein2")]
mouse_ppi<- unique(mouse_ppi)

# revise the Ensembl protein ID
d1<- substr(mouse_ppi$protein1,start = 7,stop = 24)
mouse_ppi$protein1<- d1
d1<- substr(mouse_ppi$protein2,start = 7,stop = 24)
mouse_ppi$protein2<- d1

# create a sort function
p_sort<- function(object){
  object<- object[order(object)]
  object<- paste(object[1],object[2],sep = '_')
  return(object)
}

# exclude directed ppi in STRING protein network data
mouse_ppi2<- apply(mouse_ppi, 1, p_sort)
mouse_ppi$protein<- mouse_ppi2
mouse_ppi1<- unique(mouse_ppi$protein)
mouse_ppi$protein_c<- paste(mouse_ppi$protein1,mouse_ppi$protein2,sep = '_')
d1<- which(mouse_ppi$protein == mouse_ppi$protein_c)
mouse_ppi<- mouse_ppi[d1,]
mouse_ppi<- mouse_ppi[,c("protein1","protein2")]
# obtain 5,972,403 ppi

# annotate ppi with protein and gene information in STRING and Ensembl
p_all<- c(mouse_ppi$protein1,mouse_ppi$protein2)
p_all<- unique(p_all)
p_all<- data.frame(protein_id = p_all,gene_id_length = 0,gene_id = 'NA',db_name = 'NA',EntrezGene_id = 'NA',pinfo_name = 'NA',stringsAsFactors = F)

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
mouse_gene2ensembl<- readRDS(file = 'data/mouse_gene2ensembl.rds')

# annotate ppi with gene information in NCBI (gene2ensembl)
p_all1<- p_all[p_all$protein_id %in% mouse_gene2ensembl$Ensembl_protein,]
p_all1<- p_all1[,c("protein_id","gene_id","EntrezGene_id")]
mouse_gene2ensembl1<- mouse_gene2ensembl[mouse_gene2ensembl$Ensembl_protein %in% p_all1$protein_id,]
rownames(mouse_gene2ensembl1)<- mouse_gene2ensembl1$Ensembl_protein
mouse_gene2ensembl1<- mouse_gene2ensembl1[p_all1$protein_id,]
all(p_all1$protein_id == mouse_gene2ensembl1$Ensembl_protein)
p_all1<- mouse_gene2ensembl1
p_all1<- p_all1[,c(3,2,1)]
rownames(p_all1)<- 1:nrow(p_all1)


p_all2<- p_all[!p_all$protein_id %in% mouse_gene2ensembl$Ensembl_protein,]
p_all2<- p_all2[!p_all2$gene_id == 'NA',]
p_all3<- p_all2[!p_all2$gene_id %in% mouse_gene2ensembl$Ensembl_gene_identifier,]
p_all2<- p_all2[p_all2$gene_id %in% mouse_gene2ensembl$Ensembl_gene_identifier,]
p_all2<- p_all2[,c(1,3)]
mouse_gene2ensembl1<- mouse_gene2ensembl[mouse_gene2ensembl$Ensembl_gene_identifier %in% p_all2$gene_id,]
mouse_gene2ensembl1<- mouse_gene2ensembl1[,c(2,1)]
mouse_gene2ensembl1<- unique(mouse_gene2ensembl1)

rownames(mouse_gene2ensembl1)<- 1:nrow(mouse_gene2ensembl1)
d1<- mouse_gene2ensembl1[mouse_gene2ensembl1$Ensembl_gene_identifier == 'ENSMUSG00000072694' & mouse_gene2ensembl1$Gene_id == '330173',]
mouse_gene2ensembl1<- mouse_gene2ensembl1[-as.integer(rownames(d1)),]

rownames(mouse_gene2ensembl1)<- 1:nrow(mouse_gene2ensembl1)
d1<- mouse_gene2ensembl1[mouse_gene2ensembl1$Ensembl_gene_identifier == 'ENSMUSG00000079033' & mouse_gene2ensembl1$Gene_id == '105980076',]
mouse_gene2ensembl1<- mouse_gene2ensembl1[-as.integer(rownames(d1)),]

p_all2$gene_id1<- 'NA'
colnames(p_all2)<- colnames(p_all1)

for (i in 1:nrow(p_all2)) {
  print(i)
  d1<- p_all2[i,]
  d2<- mouse_gene2ensembl1[mouse_gene2ensembl1$Ensembl_gene_identifier == d1$Ensembl_gene_identifier,]
  p_all2[i,"Gene_id"]<- d2$Gene_id
}

p_all3<- p_all3[,c(1,3,5)]
colnames(p_all3)<- colnames(p_all1)

p_all1<- rbind(p_all1,p_all2,p_all3)
rm(p_all2,p_all3,d1,d2,mouse_gene2ensembl1,i)
p_all<- p_all1
# obtain 20,140 unique proteins

# match mouse_ppi
mouse_ppi<- mouse_ppi[mouse_ppi$protein1 %in% p_all$Ensembl_protein,]
mouse_ppi<- mouse_ppi[mouse_ppi$protein2 %in% p_all$Ensembl_protein,]
# obtain 5,602,309 ppi

# load processed gene information from NCBI Gene database(Homo_sapiens.gene_info, updated in 2020.04.28)
mouse_gene_info<- readRDS(file = 'data/mouse_gene_info.rds')

# annotate
p_all$gene_symbol<- 'NA'
p_all$des<- 'NA'

for (i in 1:nrow(p_all)) {
  print(i)
  d1<- p_all[i,]
  d2<- mouse_gene_info[mouse_gene_info$GeneID == d1$Gene_id,]
  d2<- d2[,c("Symbol","description")]
  d2<- unique(d2)
  p_all[i,"gene_symbol"]<- d2$Symbol
  p_all[i,"des"]<- d2$description
}

p_all1<- p_all
rownames(p_all1)<- p_all1$Ensembl_protein
p_all1<- p_all1[mouse_ppi$protein1,]
mouse_ppi<- cbind(mouse_ppi,p_all1)
mouse_ppi<- mouse_ppi[,-3]
colnames(mouse_ppi)[3:6]<- c('pro_gene1','gene1_id','gene1_symbol','gene1_description')

p_all1<- p_all
rownames(p_all1)<- p_all1$Ensembl_protein
p_all1<- p_all1[mouse_ppi$protein2,]
mouse_ppi<- cbind(mouse_ppi,p_all1)
mouse_ppi<- mouse_ppi[,-7]
colnames(mouse_ppi)[7:10]<- c('pro_gene2','gene2_id','gene2_symbol','gene2_description')
mouse_ppi<- mouse_ppi[,c(1,2,3,7,4,8,5,9,6,10)]

# classify proteins into poteintial ligands and receptors manually
# load classified protein data
mouse_potential_lr<- readRDS(file = 'data/mouse_potential_lr.rds')

# obtain 2,005 potential ligands and 4,124 potential receptors
mouse_ligand<- mouse_potential_lr[mouse_potential_lr$Non_lr_gene_manual == 'ligand',]
mouse_receptor<- mouse_potential_lr[mouse_potential_lr$Non_lr_gene_manual == 'receptor',]
mouse_receptor<- mouse_receptor[mouse_receptor$type_of_gene == 'protein-coding',]

# match ligands and receptors
mouse_ppi1<- mouse_ppi[mouse_ppi$protein1 %in% mouse_ligand$Ensembl_protein & mouse_ppi$protein2 %in% mouse_receptor$Ensembl_protein,]
mouse_ppi2<- mouse_ppi[mouse_ppi$protein1 %in% mouse_receptor$Ensembl_protein & mouse_ppi$protein2 %in% mouse_ligand$Ensembl_protein,]

# remove duplicayed LR pairs
p_sort<- function(object){
  object<- object[order(object)]
  object<- paste(object[1],object[2],sep = '_')
  return(object)
}

d1<- mouse_ppi1[,c("gene1_symbol","gene2_symbol")]
d2<- apply(d1, 1, p_sort)
mouse_ppi1$com_gene<- d2

d1<- mouse_ppi2[,c("gene1_symbol","gene2_symbol")]
d2<- apply(d1, 1, p_sort)
mouse_ppi2$com_gene<- d2

rownames(mouse_ppi1)<- 1:nrow(mouse_ppi1)
d1<- mouse_ppi1[,c("gene1_symbol","gene2_symbol","com_gene")]
d1<- unique(d1)
mouse_ppi1<- mouse_ppi1[rownames(d1),]

rownames(mouse_ppi2)<- 1:nrow(mouse_ppi2)
d1<- mouse_ppi2[,c("gene1_symbol","gene2_symbol","com_gene")]
d1<- unique(d1)
mouse_ppi2<- mouse_ppi2[rownames(d1),]

mouse_ppi3<- mouse_ppi2[!mouse_ppi2$com_gene %in% mouse_ppi1$com_gene,]
mouse_ppi3<- mouse_ppi3[,c(2,1,4,3,6,5,8,7,10,9,11)]
colnames(mouse_ppi3)<- colnames(mouse_ppi1)
mouse_ppi1<- rbind(mouse_ppi1,mouse_ppi3)
colnames(mouse_ppi1)<- c('ligand','receptor',
                         'ligand_ensembl','receptor_ensembl',
                         'ligand_gene_id','receptor_gene_id',
                         'ligand_gene_symbol','receptor_gene_symbol',
                         'ligand_description','receptor_description','com_gene')
# obtain 255,584 potential LR pairs

# load uniprot protein knowledegbase
mouse_uniprot<- readRDS('data/mouse_uniprot.rds')





