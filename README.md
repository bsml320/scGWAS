# scGWAS: scRNA-seq assisted GWAS analysis

scGWAS leverages scRNA-seq data to identify the genetically mediated associations between traits and cell types. It requires user-provided files such as one with gene-based z-scores from GWAS, one with cellular expression matrix, one with gene-gene relationships, as well as several other parameters. All the input information is provided as a configure file to scGWAS.

Running example:

java -jar scGWAS_v3.jar CAD-Resting-Heart-Rate_Eppinga_2016.DER20.configure.txt > running.log

Example configure file:

run_model=node
normalization_model=calibration
module_score=penalty
gwas_node_file=/path/to/CAD-Resting-Heart-Rate_Eppinga_2016.z_magma.boxcox.txt
scrn_expr_file=/path/to/DER20.avg.tsv
network_file=/path/to/PathwayCommons12.All.hgnc.exPCDHA.MHC.tsv
r_include=0.1
r_exclude=0.05
remove_zero=true
outfile=/path/to/CAD-Resting-Heart-Rate_Eppinga_2016.DER20
permutation=true
verbose=true
exclude_genes_file=/path/to/HSIAO_HOUSEKEEPING_GENES.txt
