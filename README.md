# scGWAS: scRNA-seq assisted GWAS analysis

scGWAS leverages scRNA-seq data to identify the genetically mediated associations between traits and cell types. 

## Install

scGWAS is a java package and does not need to install. Users only need to set up the Java Running Environment, which can be downloaded here (https://www.oracle.com/java/technologies/downloads/).

## Running example:

java -jar scGWAS_v4.jar CAD-Resting-Heart-Rate_Eppinga_2016.DER20.configure.txt > running.log

## Notes

- The JAR package is in the folder code
- A runnable example can be found in the folder example
- The folder analysis includes all codes for preparing the input files, post-processing including the proportional test, exploration of random modules and NES calculaiton, and figure preparation.
- We provide 18 scRNA-seq panels that were already pre-processed and ready for applications. These panels were collected for 9 major tissues and can be sufficient for the majority of complex diseases and traits. If you still need to process your own panel of scRNA-seq data, please following the code in the folder analysis.
- We also provide other supporting files, including the reference network and house-keeping genes that were generally suggested to be excluded.

## Citation

The underlying method is described in Jia P et al. Landscape of trait-cell type associations by integrating single-cell transcriptomics-wide and genome-wide association studies (scGWAS)

## Contact

Please contact Peilin Jia (peilin.jia@gmail.com) or Zhongming Zhao (zhongming.zhao@uth.tmc.edu) if you have any questions.
