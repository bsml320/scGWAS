#Example command line:

# Step 1:
magma.exe --annotate window=50,35 --snp-loc ..\path\to\snp_loc.txt --gene-loc ..\path\to\NCBI37.3\NCBI37.3.gene.loc --out clozuk_pgc2

# Step 2:
magma.exe --bfile ..\g1000_eur\g1000_eur --pval clozuk_pgc2.meta.sumstats.2.txt N=105318 --gene-annot clozuk_pgc2.genes.annot.txt --out clozuk_pgc2.Magma

# The output file: XXX.genes.out.txt will be used in the next step.
