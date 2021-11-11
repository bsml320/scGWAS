# scGWAS: scRNA-seq assisted GWAS analysis

scGWAS leverages scRNA-seq data to identify the genetically mediated associations between traits and cell types. It requires user-provided files such as one with gene-based z-scores from GWAS, one with cellular expression matrix, one with gene-gene relationships, as well as several other parameters. All the input information is provided as a configure file to scGWAS.

Running example:

java -jar scGWAS_v3.jar CAD-Resting-Heart-Rate_Eppinga_2016.DER20.configure.txt > running.log

Example configure file can be found here: example/CAD-Resting-Heart-Rate_Eppinga_2016.DER20.configure.txt

Explanation of the items in the configure file:
- run_model, normalization_model (best using calibration), and module_socre (best using penalty) are the parameters for model running and should be provided with the values as in the example file.
- gwas_node_file: the file with boxcox-transformed z-scores per gene from GWAS. An example file has been provided in the example folder
- scrn_expr_file: the file with cellular gene expression per cell type from scRNA-seq data. Each row represents a gene while each column represents one cell typpe. The expression value should be the average log-transformed CPM for each gene in each cell type (average across all cells assigned to the corresponding cell type). Note: All genes that pass a pre-defined threshold (e.g., with the zero-value in <95% cells) should be included for the estimation of the null distribution. This is different from most scRNA-seq related analyses where only the highly variably expressed genes are selected. Here we will need all expressed genes to estimate the parameters for the null distribution.
- network_file: the file including gene-gene relationships. We used the PathwayCommons collection as this collection contains various types of relationships such as catalysis, chemical effect, regulation of expression or phosphorylation, react, and interacts-with, among others. Each row represents a pair of genes. The interactions that involves MHC genes have been excluded.
- r_include (suggested value: 0.1), r_exclude (0.05), remove_zero (true) are parameters for model running and should be left as it is.
- outfile: the folder name (note: this is a folder name, not a file name) where the output files (as well as the intermediate files) will be saved. This folder should be created before running of scGWAS.
- exclude_genes_file: genes in this file will be excluded from the analysis, e.g., those annotated as house-keeping genes.
- permutation: whether to conduct the virtual search process or not. Suggest to be true.
- verbose: whether to print the intermediate results on screen.
