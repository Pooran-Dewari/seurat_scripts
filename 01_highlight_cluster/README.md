Motivation: Create an output dir, plot and output umaps with every single cluster highlighted rest grayed out  

Assumptions: Assumes that seurat.filtered object is accessible with umap assay in it  

The R functions attached will 
- first check if output_dir exists
- if it doesn't, it will create the dir and output results in there
- if it exists, it will create output_dir1, and output results in there
- if output_dir1 exists, it will create output_dir2, output results in there, so on and so forth


