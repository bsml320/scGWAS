All parameters are provided through the configure file. The general rules of this file include:
- One parameter per row
- The parameters and their values are case sensitive. Thus, r_include and R_include would be two different parameters.
- No space. If the file names of your input files have space, please rename it.
- The format is: the row starts with a keyword followed exactly by its values without space, e.g., r_include=0.1
- A full list of all acceptable parameters is explained below. However, many parameters can be leave it with the default values. Thus, in each application, users may only provide the following "essential" parameters:
  -- gwas_node_file
  -- scrn_expr_file
  -- network_file
  -- outfile

 
Explanation of all parameters:

- run_model: always set as "node", indicating node-weighted module search
- normalization_model: default: calibration; available values: calibration, scale; suggestion: always set as calibration
- module_socre: default: penalty; available values: penalty, tw. When module_score=penalty, it indicates m=mg+mv-sd(mg,mv). When module_score=tw, it indicates m = mg+mv. 
- gwas_node_file: the file with boxcox-transformed z-scores per gene from GWAS. Two-column file (column 1: gene symbole; column 2: z-score); tab-seperated; no header line. An example file has been provided in the example folder
- scrn_expr_file: the file with cellular gene expression per cell type from scRNA-seq data. Each row represents a gene while each column represents one cell typpe. The expression value should be the average log-transformed CPM for each gene in each cell type (average across all cells assigned to the corresponding cell type). Note: All genes that pass a pre-defined threshold (e.g., with the zero-value in <95% cells) should be included for the estimation of the null distribution. This is different from most scRNA-seq related analyses where only the highly variably expressed genes are selected. Here we will need all expressed genes to estimate the parameters for the null distribution.
- network_file: the file including gene-gene relationships. We used the PathwayCommons collection as this collection contains various types of relationships such as catalysis, chemical effect, regulation of expression or phosphorylation, react, and interacts-with, among others. Each row represents a pair of genes. The interactions that involves MHC genes have been excluded.
- r_include (suggested value: 0.1), r_exclude (0.05), remove_zero (true) are parameters for model running and should be left as it is.
- outfile: the folder name (note: this is a folder name, not a file name) where the output files (as well as the intermediate files) will be saved. This folder should be created before running of scGWAS.
- exclude_genes_file: genes in this file will be excluded from the analysis, e.g., those annotated as house-keeping genes.
- permutation: whether to conduct the virtual search process or not. Suggest to be true.
- verbose: whether to print the intermediate results on screen.

Regarding the input files: Java can be quite confusing when looking for files. If not sure, use the full path to each of the input files.
