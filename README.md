# scGWAS: scRNA-seq assisted GWAS analysis

scGWAS leverages scRNA-seq data to identify the genetically mediated associations between traits and cell types. It requires user-provided files such as one with gene-based z-scores from GWAS, one with cellular expression matrix, one with gene-gene relationships, as well as several other parameters. All the input information is provided as a configure file to scGWAS.

Running example:

java -jar scGWAS_v3.jar CAD-Resting-Heart-Rate_Eppinga_2016.DER20.configure.txt > running.log

Example configure file can be found here

