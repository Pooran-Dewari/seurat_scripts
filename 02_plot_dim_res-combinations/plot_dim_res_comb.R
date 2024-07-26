#Motivation: It is quite a challenge to optimise dims and resolution parameters for clustering and UMAP plotting
# this script can create all possible combinations of user-supplied dims and resolutions parameters
# for example, the script below will create plots using all possible combinations of dims (20,25,30) and res (0.5, 0.6, 0.7, 0.8, 0.9) values


# code STARTS ...........................................

library(Seurat) # Ensure Seurat is loaded
library(ggplot2)
library(glue)

dim_lo <- 1
dim_up <- c(20, 25, 30)
res <- c(0.5, 0.6, 0.7, 0.8, 0.9)

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

# step 1: Create the directory for saving plots
output_dir <- create_output_directory("plot_dim_res_combinations")

# Function to generate plots and save them as PDFs
generate_plots <- function(dim_up_value, res_value, output_dir) {
  # Adjust the Seurat object according to the current parameters
  seurat.filtered <- FindNeighbors(seurat.filtered, dims = dim_lo:dim_up_value)
  seurat.filtered <- FindClusters(seurat.filtered, resolution = res_value)
  seurat.filtered <- RunUMAP(seurat.filtered, dims = dim_lo:dim_up_value)
  
  # Create the plot
  plot <- DimPlot(seurat.filtered, reduction = "umap", label = TRUE,
                  label.size = 4,
                  label.color = "black",
                  label.box = TRUE,
                  repel = TRUE) + 
    NoLegend() +
    ggtitle(glue("dims = {dim_up_value}, resolution={res_value}"))
  
  # Open a PDF device
  pdf_filename <- file.path(output_dir, glue("dims_{dim_up_value}_res_{res_value}.pdf"))
  pdf(pdf_filename)
  
  # Print the plot to the PDF
  print(plot)
  
  # Close the PDF device
  dev.off()
}

# step 2: Apply the plotting function to each combination
lapply(1:nrow(combinations), function(i) {
  generate_plots(combinations$dim_up[i], combinations$res[i], output_dir)
})


# code ENDS ...........................................

