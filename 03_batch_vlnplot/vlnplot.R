#Motivation: create vlnplots for all marker genes for a Seurat cluster


#................... code STARTS ....................

library(Seurat) # Ensure Seurat object is loaded
library(ggplot2)
library(glue)
library(dplyr)

# edit clust_num below as you wish!!
clust_num <- 3

# Function to create the output directory and return its name
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

# Step 1: Create the directory for saving plots
output_dir <- create_output_directory("vlnplots")


# Step 2: get marker genes for cluster of interest
feat <- seurat.filtered.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1.5, cluster == clust_num, p_val_adj < 0.05)

# Function to generate plots and save them as PNGs
vlnplot_fun <- function(x, output_dir) {
  plot <- VlnPlot(seurat.filtered, features = x) + 
    ggtitle(glue("{x}")) +
    NoLegend()
  
  # Define the filename for the PNG
  png_filename <- file.path(output_dir, glue("{x}.png"))
  
  # Open a PNG device with high resolution settings
  png(png_filename, width = 10, height = 7, units = "in", res = 300)
  
  # Print the plot to the PNG device
  print(plot)
  
  # Close the PNG device
  dev.off()
}

# Step 3: Apply the plotting function to each gene in feat
lapply(feat$gene, function(x) {
  vlnplot_fun(x, output_dir)
})

#................... code ENDS ....................
