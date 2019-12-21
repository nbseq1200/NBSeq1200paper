# NBSeq assessment module


# 1. Download and install

### Download and unzip the attached package:

evaluation_package.zip

# 2. Usage

## Option 1: run the batch script 

### In the unzipped folder run the script, the script takes one parameters: 

   <answer key>

```
cd evaluation_package
./run_assessment.sh data/example_answer_key.tsv
```

### The script generates two folders, one blinded and one unblinded:

- result_**example_answer_key**_UCB_blinded
- result_**example_answer_key**_unblinded

## Option 2: run the perl script

```
cd evaluation_package
perl src/new_assess_180910.pl --dic data/dictionary.v10.studyID2.tsv --key data/example_answer_key.tsv --data example_pipeline --out result
```

| Required parameters:         | Descriptions                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| --dictionary or dic \<file\> | dictionary file in tsv format<br>columns: DiseaseID, Screen, GeneList(comma separated) |
| --key \<file\>               | answer key file in tsv format<br>columns: sample, DiseaseID  |
| --data \<dir\>               | directory which contains the pipeline prediction results     |
| --out \<file\>               | output directory name<br>(see details in 3. Output Files)    |

| Optional parameters: | Descriptions                                                 |
| -------------------- | ------------------------------------------------------------ |
| --list \<file\>      | sample list<br>(if it is not provided, all samples appearing in both key file and prediction result file will be caculated) |
| --screen \<file\>    | MS/MS screen result file in tsv format<br>columns: sample, MS/MS screen result |

## Option 3: Caculate the specificity on affected set

```
perl src/spec.pl data/dictionary.v10.studyID2.tsv <path to the summary_disorder.tsv> ./pipeline <path to the var_anno> <out_dir>
```

# 3. Output Files

- assess_log.txt	[file] log file
- summary_disorder.tsv	[file] report sample counts for TP/FN/FN+FP cases and FN cases with heterozygous for each disorder and pipeline
- summary_gene.tsv	[file] report sample counts for TP/FP/FN+FP cases and FN cases with heterozygous for each gene and pipeline
- summary_pipeline.tsv	[file] report total sample counts for TP/FP/FN+FP/FN cases on disorder/screen/affected status levels for each pipeline
- summary_screen.tsv	[file] report sample counts for TP/FN/FN+FP cases and FN cases with heterozygous for each screen and pipeline
- sens	[dir] contains sensitivity file for each pipeline
  - sens_test_<pipeline>_<filter|batch>.tsv	[file] sensitivity on test dataset
  - sens_val_<pipeline>_<filter|batch>.tsv	[file] sensitivity on validation dataset
  - sens_<pipeline>_<filter|batch>.tsv	[file] sensitivity on total dataset
- spec	[dir] contains file reporting fraction of false positive rate on unaffected samples for each pipeline
  - spec_test_<pipeline>_<filter/batch>.tsv	[file] false positive rate on unaffected samples from test dataset
  - spec_val_<pipeline>_<filter/batch>.tsv	[file] false positive rate on unaffected samples from validation dataset
  - spec_<pipeline>_<filter/batch>.tsv	[file] false positive rate on unaffected samples from total dataset
