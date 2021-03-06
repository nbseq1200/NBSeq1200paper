python implementation of the screening exome analysis pipeline for NBSeq1200 paper

## Requirements
- Python 2.7
- pandas 0.24.2
- numpy 1.13.1

## Usage
```
python PipelineMain.py <input_file> <param_file> <output_directory>
```


### input_file format

The input file to the script is a tab-separated file of variants for each sample with the following columns. The file is generated by the tsv output module of the Varant variant annotation software (http://compbio.berkeley.edu/proj/varant/Home.html). 

Column name | Description
----|-----
sample| Sample ID
chrom | Chromosome
pos | Position
ref | Reference Allele
alt | Alternate Allele
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
hgmdvar | HGMD category for the variant (DM - Disease-causing mutation, DM? - Likely Disease-causing mutation)
clnsg170907 | ClinVar clinical significance of the variant (0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic)
clinstars170907 | Number of ClinVar review stars for the variant
appris | Whether transcript is a principal transcript according to Appris
rf_score | Random forest based integrated splicing effect score in dbscSNV
ug.gt | Genotype call of the variant by GATK UnifiedGenotyper
ug.gq |  Genotype quality of the variant by GATK UnifiedGenotyper
hc.gt | Genotype call of the variant by GATK HaplotypeCaller
hc.gq |  Genotype quality of the variant by GATK HaplotypeCaller
pp.gt |  Genotype call of the variant by Platypus
pp.gq |  Genotype quality of the variant by Platypus



### param_file format

Tab-separated file with the following columns. Each row defines a unique set of parameters. The file [NBSeq1200.params](param/NBSeq1200.params) in the 'param' folder defines the parameter sets used in the manuscript.

Column name | Description
----|----
id	| Descriptor for the parameter set. 
group| Additional descriptor. 
cluster | Additional descriptor. 
name	 |  Additional descriptor. 
genelist | Single column file with list of gene names
transcript|  Which transcript to select for the gene 
gqthres | Minimum allowable GQ 	
caller | Choice of genotype caller 
disease_maf_db | Source databases for population MAF for curation arm
disease_maf_thres| Population MAF threshold for database curated variants
hgmd| HGMD variant categories to include in curation arm
clnvar| Clinvar clinical significance categories to include in curation arm
clnstar| Number of ClinVar review stars to include
maf_db|  Source databases for population MAF for impact arm
maf_thres| Population MAF threshold for database impact arm variants
pa_list| List of Varant-annotated categories to include in defining protein-altering variants in impact arm  
pathogen1| 	Optional pathogenicity tool 1
pathogen1_score| Pathogenicity tool 1 pathogenicity score threshold
pathogen2| Optional pathogenicity tool 2
pathogen2_score| Pathogenicity tool 2 pathogenicity score threshold
pathogen3| Optional pathogenicity tool 3
pathogen3_score|  Pathogenicity tool 3 pathogenicity score threshold
loftee| Whether to include LOF variants as predicted by LOFTEE (Requires 'loftee' column in input_file) 
includefile| File path for list of NBSeq-curated variants to be included
excludefile| File path for list of NBSeq-curated variants to be excluded
cnvfile| File path for CNV calls (3 column file - with sample ID, gene name, and zygosity of the CNV calls)



