python implementation of the screening exome analysis pipeline for NBSeq1200 paper

## Requirements
- Python 2.7
- pandas 0.24.2
- numpy 1.13.1

## Usage
```
python PipelineMain.py <input_file> <param_file> <output_directory>
```

## input_file format
Tab-separated file of variants for each sample with the following columns.

Column name | Description
----|-----
sample| Sample ID
chrom | Chromosome
pos | Position
ref | Reference Allele
alt | alternate Allele
gene | Gene name
transcript| Gene transcript ID
mutation | Mutation type according to Varant
canon | Whether transcript is canonical transcript in Gencode V19
splice | Splice annotation according to Varant
kgaf | Population MAF of variant in 1000 Genomes Phase 3
espaf | Population MAF of variant in Exome Sequencing Project (ESP)
exacaf | Population MAF of variant in ExAC database
cadd | CADD score of the variant
meta_svm | meta_svm score of the variant in dbNSFP v2.0
hgmdvar | HGMD category for the variant 
clnsg170907 | ClinVar clinical significance of the variant (0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other)
clinstars170907 | Number of ClinVar review stars for the variant
appris | Whether transcript is a principal transcript according to Appris
rf_score | Random forest based integrated splicing effect score in dbscSNV
ug.gt | Genotype call of the variant by GATK UnifiedGenotyper
ug.gq |  Genotype quality of the variant by GATK UnifiedGenotyper
hc.gt | Genotype call of the variant by GATK HaplotypeCaller
hc.gq |  Genotype quality of the variant by GATK HaplotypeCaller
pp.gt |  Genotype call of the variant by Platypus
pp.gq |  Genotype quality of the variant by Platypus


## param_file format
Tab-separated file with the following columns. 

Column name | Description
----|----
id	| Descriptor for the parameter combination. 
group| Additional descriptor. 
cluster | Additional descriptor. 
name	 |  Additional descriptor. 
genelist | Single column file with list of gene names| 
transcript|  Transcript | 	
gqthres | Minimum allowable GQ 	
caller | Genotype calls from vcf
disease_maf_db | Databases for  	
disease_maf_thres| Population MAF threshold for disease databases
hgmd| HGMD variant classification
clnvar| 
clnstar| Number of ClinVar review stars to include
maf_db| Choice of population database for 
maf_thres| Population MAF threshold for population database
pa_list| List of terms for protein altering 
pathogen1| 	Optional pathogenicity tool 1
pathogen1_score| Optional pathogenicity tool 1
pathogen2| Optional pathogenicity tool 2
pathogen2_score| Optional pathogenicity tool 1
pathogen3| Optional pathogenicity tool 3
pathogen3_score| Optional pathogenicity tool 1
loftee| 
includefile| File path for list of NBSeq-reviewed variants to be included
excludefile| File path for list of NBSeq reviewed variants to be excluded
cnvfile| File path for 



