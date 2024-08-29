# got this script from thomas
library(Seurat)
library(tidyverse)
library(Matrix)
library (ggpubr)
library(harmony)
library(scCustomize)

setwd("~/Documents/seurat_thomas_scripts")

### STEP1: Create Seurat Object from Parse DGE unfiltered outputs ##############
#...............................................................................

mat <- ReadParseBio("DGE_unfiltered")
mat

# check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
# outputs FALSE 28307, all rows have names..only run the command below if returns TRUE
# rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv("DGE_unfiltered/cell_metadata.csv")


# create seurat object
seu_obj <- CreateSeuratObject(mat, meta.data = cell_meta)

#...............................................................................


########### STEP2: Explore Seurat Object, plus some filtering ##############
#...............................................................................

seu_obj
# An object of class Seurat 
# 28307 features across 114349 samples within 1 assay 
# Active assay: RNA (28307 features, 0 variable features)
# 1 layer present: counts

#make qc plots to see number of genes and counts 
VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA"))

#dimensions in the object, number of genes and number of cells
dim(seu_obj)
# 28307 114349

# filtering by number of genes
seu_obj<-  subset(seu_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

# check what left after filtering
dim(seu_obj)
# 28307 22208

#make qc plots to see number of genes and counts 
VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA"))

mt.genes = c("ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "CYTB")

#Remove mt.genes that are not present in the seurat seu_obj i.e. are expressed in too few cells so was filtered out in "CreateSeuratseu_obj" function above
mt.genes.present <- mt.genes[mt.genes %in% rownames(seu_obj)]
#print(paste0(length(mt.genes.present),"/",length(mt.genes), " mito genes present in the count matrix"))
seu_obj <- PercentageFeatureSet(seu_obj,features = mt.genes.present, col.name = "percent.mt")
VlnPlot(seu_obj,features = c("percent.mt"),pt.size = 0) + geom_hline(yintercept = 20) + NoLegend()

# VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#...............................................................................


########### STEP3: Data normalisation, scaling, PCA #####################
#...............................................................................

#standard seurat workflow until elbow plot

# log-transform (normalisation)
seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu_obj), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# linear transformation (‘scaling’) 
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes)

# perform linear dimensional reduction (PCA)
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
DimPlot(seu_obj, reduction = "pca") + NoLegend()
DimHeatmap(seu_obj, dims = 1:25, cells = 500, balanced = TRUE)

#use this to assess number of PC's to use in downstream analysis
ElbowPlot(seu_obj, ndims=50)


# going ahead with dims = 30 and res 0.6, this generates me late infection cluster10

seu_obj<- FindNeighbors(seu_obj, dims = 1:30)
seu_obj<- FindClusters(seu_obj, resolution = 0.6)


#plot the UMAP
seu_obj <- RunUMAP(seu_obj, dims = 1:30)
DimPlot(seu_obj, reduction = "umap", label=T)

#save rds
saveRDS(seu_obj, file = "seu_obj_umap30x6.rds")
#......................................








# Get number of cells per cluster per group ID
#seu_obj <- Seurat::StashIdent(object = seu_obj, save.name = "my.clusters")
#bygroup<-table(seu_obj@meta.data$my.clusters, seu_obj@meta.data$sample)
#bysample<- table(hk@meta.data$my.clusters, seu_obj@meta.data$sample)

#to write the tables to an excel
#write.csv(bygroup, "bygroups.csv", row.names = T)
#write.csv(bysample, "M:/seu_obj/Single Cell Nuclei/bysample.csv", row.names = T)


# seu_objharmony <- seu_obj
# #with harmony  on samples. Harmony used to deal with batch effects and helps integrate the variation seen with various factors. 
# seu_objharmony <- NormalizeData(seu_objharmony, normalization.method = "LogNormalize", scale.factor = 10000)
# seu_objharmony <- FindVariableFeatures(seu_objharmony, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(seu_objharmony)
# seu_objharmony <- ScaleData(seu_objharmony, features = all.genes)
# seu_objharmony <- RunPCA(seu_objharmony, features = VariableFeatures(object = seu_objharmony))
# seu_objharmony<- RunHarmony(seu_objharmony, group.by.vars = "sample",plot_convergence = TRUE)
# ElbowPlot(seu_objharmony, ndims=50)
# seu_objharmony <- RunUMAP(seu_objharmony,reduction = "harmony",dims = 1:30)
# seu_objharmony<- FindNeighbors(seu_objharmony,reduction = "harmony", dims = 1:30)
# seu_objharmony<- FindClusters(seu_objharmony, resolution = 0.6)
# DimPlot(seu_objharmony, reduction = "umap", label=T)


# seu_objharmony <- Seurat::StashIdent(object = seu_objharmony, save.name = "my.clusters")
# bygroup<-table(seu_objharmony@meta.data$my.clusters, seu_objharmony@meta.data$groups)
# bysample<- table(seu_objharmony@meta.data$my.clusters, seu_objharmony@meta.data$sample)
# 
# write.csv(bygroup, "M:/seu_obj/Single Cell Nuclei/bygroups_harmonysample.csv", row.names = T)
# write.csv(bysample, "M:/seu_obj/Single Cell Nuclei/bysample_harmonysample.csv", row.names = T)
# 
# 
# p1<-DimPlot(seu_objharmony, reduction = "umap", label=T)
# p2<-DimPlot(seu_objharmony, reduction = "umap", label=F,group.by="sample")
# p3<-DimPlot(seu_objharmony, reduction = "umap", label=F,group.by="groups")
# ggarrange(p1,p2,p3, ncol = 3, nrow = 1)
# 
# 
# seu_obj1<-seu_objharmony
# new.cluster.ids <- c("Muscle ? 1","Cuticle 1","2","3","4","Cuticle 2","6","Cuticle 3","Muscle 1","Muscle 2","Cuticle 4","Chitins 1","12","Muscle ? 2","14","15","Chromatophores","Chitins 2","18","19","20","21")
# names(new.cluster.ids) <- levels(seu_obj)
# seu_obj <- RenameIdents(seu_obj, new.cluster.ids)
# p<-DimPlot(seu_obj, reduction = "umap", label=T,label.size = 4) + NoLegend()
# ggsave("M:/seu_obj/Single Cell Nuclei/analysis/Figures/umap.pdf", p, width=27, height=15, 
#    device = "pdf", units = "cm")  
# 
# 
# saveRDS(seu_objharmony, file = "M:/seu_obj/Single Cell Nuclei/harmony_by_sample.rds")



#Dotplot code
# moultdb<-c("LOC113804419",
# "LOC113826024",
# "LOC113820103",
# "LOC113800214",
# "LOC113814499",
# "LOC113814856",
# "LOC113827159",
# "LOC113801169",
# "LOC113828208",
# "LOC113804417",
# "LOC113800218",
# "LOC113823986",
# "LOC113820037",
# "LOC113824971")
# DotPlot(seu_obj,cols = c("white","darkorchid4"), col.min = 0.05, col.max = 3,dot.min = 0.05,  dot.scale = 8, features=moultdb)  + 
# theme(axis.text.x = element_text(angle = 45,vjust=0.8,hjust=0.8))  + coord_flip()
# 
# 
# 
# #Featureplot
# FeaturePlot(seu_obj, features = "LOC113804419")


###############................. find marker genes ################################

seu_obj <- readRDS(file = "seu_obj_umap30x6.rds")


#Finding markers for clusters and then plotting the top 5 with hierachial clustering
#all_markers <- FindAllMarkers(object = seu_obj, min.pct = 0.25, log2fc.threshold = 0.5)
all_markers <- FindAllMarkers(object = seu_obj, min.pct = 0.25, log2fc.threshold = 0.5)

# filter out expressed >1.5
all_markers_filt <- all_markers %>%
  dplyr::filter(avg_log2FC >= 1.5) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  rownames_to_column(var = "gene_id")

# count marker genes per cluster
all_markers_filt %>% 
  count(cluster)

# marker for cluster 10
all_markers_cluster_10 <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1.5 & cluster == 10)


VlnPlot(seu_obj, features = c("G20535", "G6643"))
FeaturePlot(seu_obj, features = c("G20535", "G6643"))

#write the marker list to file
write_tsv(all_markers_filt, "all_markers_filt.tsv")

#............................................................................



#get gene symbols from ncbi
# https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/29159/
# select all columns, check box next to Gene ID to select all
# download "ncbi_dataset_taxon_29159.tsv"

#load ncbi gene id table
features_ncbi <- read_tsv("ncbi_dataset_taxon_29159.tsv", col_names = T) %>% 
  dplyr::rename(gene_id = "Nomenclature ID" )

# format chr to match gtf/gff3 files, need later to getFasta
features_ncbi$Chromosomes <- case_when(
  features_ncbi$Chromosomes == "1" ~ "LR761634.1",
  features_ncbi$Chromosomes == "2" ~ "LR761635.1",
  features_ncbi$Chromosomes == "3" ~ "LR761636.1",
  features_ncbi$Chromosomes == "4" ~ "LR761637.1",
  features_ncbi$Chromosomes == "5" ~ "LR761638.1",
  features_ncbi$Chromosomes == "6" ~ "LR761639.1",
  features_ncbi$Chromosomes == "7" ~ "LR761640.1",
  features_ncbi$Chromosomes == "8" ~ "LR761641.1",
  features_ncbi$Chromosomes == "9" ~ "LR761642.1",
  features_ncbi$Chromosomes == "10" ~ "LR761643.1",
  features_ncbi$Chromosomes == "MT" ~ "MZ497416.1"
)
# load all_genes.csv from Parse output
#genes_parse <- read_csv("DGE_unfiltered/all_genes.csv")


all_markers_filt_geneID <- left_join(all_markers_filt, features_ncbi, by = "gene_id")

#write the marker list to file
write_tsv(all_markers_filt_geneID, "all_markers_filt_geneID.tsv")

###############################............................................


all_markers_filt_geneID %>% 
  select(1, 18, 21,22)


# wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-58/gff3/crassostrea_gigas/Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz
# gunzip Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3.gz 






# get lit of top 10 markers for each cluster
top_markers_20 <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 20, named_vector = FALSE,
                                  make_unique = TRUE, data_frame = TRUE)

write_tsv(top_markers_20, "top_markers_20.tsv")

Clustered_DotPlot(seurat_object = seu_obj, features = top_markers)



### get bed file coordinates for top_markers_20
top_genes <- top_markers_20$gene

read_table("Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3")

gff <- as.data.frame(rtracklayer::import("Crassostrea_gigas.cgigas_uk_roslin_v1.58.chr.gff3"))

gff_v2 %>% 
  filter(end - start >3)

gff_v2 <- gff %>% 
  select(1:3, 5, 7, 10) %>% 
  filter(type == "CDS") %>% 
  separate(col = "ID", sep = ":", into = c("a", "gene_ID", "c"), remove = FALSE) %>% 
  select(1:6,8) %>% 
  separate(col = "gene_ID", sep = "\\.", into = c("gene", "b"), remove = FALSE) %>% 
  select(1:4,6,8) %>% 
  filter(end - start >30) 

top_markers_20_bed <- inner_join(gff_v2, top_markers_20, by = "gene", relationship = "many-to-many") %>% 
  select(1:3,5,8,4) %>% 
  dplyr::rename(name = "ID")

write_tsv(top_markers_20_bed, "top_markers_20_bed.bed", col_names = FALSE)

## getFasta for top markes 20 bed, bash programming
#bedtools getfasta -name -s -fi Crassostrea_gigas_uk_roslin_v1.dna_sm.primary_assembly.fa -bed top_markers_20_bed.bed -fo top_markers_20.fasta


#######################....... translate all three frames using biostrings
library(Biostrings)


# import fasta file into biostring
fa_seqs  <- readDNAStringSet("top_markers_20.fasta")

# get all three position (1,2,3) as start
fa_subseqs <- lapply(1:3, function(pos) 
  subseq(fa_seqs, start=pos))

# translate
aa_seqs_list <- lapply(fa_subseqs, translate)


# Combine translated sequences into a single AAStringSet
# Initialize an empty AAStringSet
combined_aa_seqs <- AAStringSet()

# Loop to combine sequences from each frame and filter out those with stop codons
for (i in seq_along(aa_seqs_list)) {
  # Name each sequence
  names(aa_seqs_list[[i]]) <- paste0(names(fa_seqs), "_Frame_", i)
  
  # Filter sequences to remove those containing stop codons
  valid_seqs <- aa_seqs_list[[i]][!grepl("\\*", aa_seqs_list[[i]])]
  
  # Combine filtered sequences
  combined_aa_seqs <- c(combined_aa_seqs, valid_seqs)
}

# Step 5: Export the filtered amino acid sequences to a FASTA file
writeXStringSet(combined_aa_seqs, filepath = "filtered_translated_sequences.fasta")

# upload to egnog database for orhtology finding
#http://eggnog-mapper.embl.de/


##########.............................................................................

# load the top 20 marker genes file
top_markers_20 <- read_tsv ("top_markers_20.tsv")

# load eggnog output file
eggnog <- read_tsv("eggnog_out_from_aa/MM_tnx2xt46.emapper.annotations.tsv", comment = "##")


eggnog <- eggnog %>%
  separate(`#query`, sep = "\\.", into = c("a", "b"), remove = FALSE) %>% 
  separate(`a`, sep = ":", into = c("c", "gene")) 

# left join
top_markers_20_ann <- dplyr::left_join(top_markers_20, eggnog, 
                                       by = "gene", relationship = "many-to-many")

write_tsv(top_markers_20_ann, "top_markers_20_ann.tsv")

# vlnplot for each cluster

top_markers_10 <- top_markers %>% 
  filter(cluster == 10)

VlnPlot(seu_obj, features = top_markers_10$gene)
FeaturePlot(seu_obj, features = top_markers_10$gene)



# #I have my genome annotation file that I then use to add gene names and other annotations self done from papers to the gene IDs
# de_file <- read.csv("M:/seu_obj/Single Cell Nuclei/analysis/Figures/Moultdb/all_markers.csv", header = TRUE)
# ensemble_attributes <- read.csv("M:/seu_obj/Single Cell Nuclei/analysis/Gao2017/seu_obj_genome_withGao_Function.csv", header = TRUE, ",")
# final_output <- de_file %>% left_join(ensemble_attributes, by =c("gene"="Gene_ID"))
# write.csv(final_output,"M:/seu_obj/Single Cell Nuclei/analysis/Figures/Moultdb/all_markers.csv", , na="") 

