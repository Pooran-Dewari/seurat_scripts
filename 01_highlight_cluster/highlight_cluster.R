# MOTIVATION: create an output dir, plot and output umaps with every single cluster highlighted rest grayed out

# ASSUMPTIONS: seurat.filtered object is accessible that has umap assay in it, for example, seurat object created with bits of Seurat code below
#seurat.filtered <- FindNeighbors(seurat.filtered, dims = 1:25)
#seurat.filtered <- FindClusters(seurat.filtered, resolution = 0.5)
#seurat.filtered <- RunUMAP(seurat.filtered, dims = 1:25)
#seurat.filtered.markers <- FindAllMarkers(seurat.filtered, only.pos = TRUE)

# WHAT THESE FUNCTIONS DO?
# the functions below will first check if output_dir exists
# if it doesn't, it will create the dir and output results in there
# if it exists, it will create output_dir1, and output results in there
# if output_dir1 exists, it will create output_dir2, output results in there, so on and so forth...

# code STARTS ..........................................

library(Seurat)
library(tidyverse)

# part 1: Function to create the output directory and return its name
create_output_directory <- function(base_name) {
  counter <- 0
  dir_name <- as.character(base_name)
  
  # Check if the directory exists, increment the counter if it does
  while (file.exists(dir_name)) {
    counter <- counter + 1
    dir_name <- paste0(base_name, counter)
  }
  
  # Create the directory
  dir.create(dir_name)
  cat("Directory created:", dir_name, "\n")
  
  return(dir_name)
}


# part 2: Function to highlight a cluster and save a plot in the specified directory
plot_highlight_cluster <- function(x, output_dir) {
  
  # define variables to categorise based on cluster of interest, also for outputting file names
  query <- x
  cl <- paste0("cluster", query)
  ot <- "others"
  
  # add cluster information metadata in a new column, this will categorise clusters
  seurat.filtered$cluster_info <- ifelse(seurat.filtered$seurat_clusters == query, cl, ot)
  
  # Create a PDF in the specified directory
  pdf(file.path(output_dir, paste0(cl, ".pdf")))
  
  # Plotting (assumes appropriate packages and data are available)
  plot <- DimPlot(seurat.filtered, reduction = "umap", group.by = "cluster_info", split.by = "SampleGroup",
                  cols = c("red", "gray"), order = c(ot, cl), ncol = 2) +
    ggtitle(paste(cl, " across all samples"))
  
  # Print the plot to the PDF
  print(plot)
  
  # Close the PDF device
  dev.off()
}


# part 3: Main function to create the directory and generate plots
main_function <- function() {
  
  # create the output directory once, EDIT argument if needs be
  output_dir <- create_output_directory("cluster_out_dim25")
  
  # cluster count for the lapply loop
  cluster_count = length(levels(Idents(seurat.filtered)))-1
  
  # generate individual plots for all clusters, save them in the output_dir
  lapply(0:cluster_count, plot_highlight_cluster, output_dir)
}

# part 4: Run the main function, this will create the directory and output plots as pdf
main_function()

# code ENDS ..........................................

