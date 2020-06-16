library(readr)
library(biomaRt)
#devtools::install_github("MaayanLab/genesetr")
library(genesetr)
library(OmnipathR)


#download biogrid - https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-3.5.186/BIOGRID-ALL-3.5.186.mitab.zip

#read in the biogrid network
BIOGRID_ALL_3_5_186_mitab <- read_delim("BIOGRID-ALL-3.5.186.mitab.txt", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)

#keep human interactions
BIOGRID_ALL_3_5_186_mitab <- BIOGRID_ALL_3_5_186_mitab %>% filter(`Taxid Interactor A` == "taxid:9606" & `Taxid Interactor B` == "taxid:9606")

#get the gene symbols for the interactions
BIOGRID_ALL_3_5_186_mitab$A <- sapply(strsplit(BIOGRID_ALL_3_5_186_mitab$`Alt IDs Interactor A`,"|",fixed=T),"[[",2)
BIOGRID_ALL_3_5_186_mitab$B <- sapply(strsplit(BIOGRID_ALL_3_5_186_mitab$`Alt IDs Interactor B`,"|",fixed=T),"[[",2) 
BIOGRID_ALL_3_5_186_mitab$A <- gsub("entrez gene/locuslink:","",BIOGRID_ALL_3_5_186_mitab$A)
BIOGRID_ALL_3_5_186_mitab$B <- gsub("entrez gene/locuslink:","",BIOGRID_ALL_3_5_186_mitab$B)
biogrid <- data.frame(A=BIOGRID_ALL_3_5_186_mitab$A,B=BIOGRID_ALL_3_5_186_mitab$B)
biogrid$A <- genesetr::HGNCapproved(as.character(biogrid$A))
biogrid$B <- genesetr::HGNCapproved(as.character(biogrid$B))

#get the omnipath interactions excluding the Wang dataset as is noisy.
resources <- get_interaction_databases()
resources <- resources[ resources!="Wang"]
omnipath <- import_AllInteractions( filter_databases=resources)[,3:4]
colnames(biogrid) <- colnames(omnipath)
network <- rbind(omnipath,biogrid)

#make sure the gene symbols are the current official ones
network$source_genesymbol <- genesetr::HGNCapproved(network$source_genesymbol)
network$target_genesymbol <- genesetr::HGNCapproved(network$target_genesymbol)

#filter the network to keep just the protein coding genes
ensembl=useMart("ensembl")
ensemblHuman = useDataset("hsapiens_gene_ensembl",mart=ensembl)
humanProteinCodingGenes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), filters='biotype', values=c('protein_coding'), mart=ensemblHuman)

network <- network[ network$source_genesymbol %in% humanProteinCodingGenes$external_gene_name,]
network <- network[ network$target_genesymbol %in% humanProteinCodingGenes$external_gene_name,]

write.table(network, file="network.txt", sep = "\t",row.names=F,quote = F,col.names = T)


