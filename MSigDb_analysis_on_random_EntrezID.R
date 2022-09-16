# 1. Generate NCBI entrez gene ID for human

## download gene_info.gz from NCBI website (https://ftp.ncbi.nih.gov/gene/DATA/)

awk ' $1 ~ 9606 {print $1, $2, $3}' gene_info > human_gene_info_geneID.txt
# column 1 tax_id for human (9606). Select Homo sapiens (human), species, primates. 

# 2. Import the df in R
