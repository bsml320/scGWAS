# scGWAS: scRNA-seq assisted GWAS analysis

scGWAS leverages scRNA-seq data to identify the genetically mediated associations between traits and cell types. 

## Install

scGWAS is a java package and does not need to install. Users only need to set up the Java Running Environment, which can be downloaded here (https://www.oracle.com/java/technologies/downloads/).

## Running example:

java -jar scGWAS_v3.jar CAD-Resting-Heart-Rate_Eppinga_2016.DER20.configure.txt > running.log

## Notes

- A runnable example can be found in the folder example
- The code to prepare the input files is provided in the folder analysis
- The code for post-process can also be found in analysis, including the proportional test, exploration of random modules and NES calculaiton.
- 

## Citation

The underlying method is described in Jia P et al. Landscape of trait-cell type associations by integrating single-cell transcriptomics-wide and genome-wide association studies (scGWAS)

## Contact

Please contact Peilin Jia (peilin.jia@gmail.com) or Zhongming Zhao (zhongming.zhao@uth.tmc.edu) if you have any questions.
