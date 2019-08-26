## Ensembl_BED

#### Input   
- FASTA for chromosomes  
- GFF3 for annotations

#### Script action  
This script will:  
- Split up an Ensembl GFF3 by locus type (e.g. CDS, UTR, lncRNA, etc.)  
- Generate BED files for unannotated loci (e.g. introns, intergenic regions, noncoding genes)

##### Script editing  
- Currently, have to hard code in FA,GFF3, and desired annotation types  
