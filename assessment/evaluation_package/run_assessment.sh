key=$1
name=$(basename $key)
name=${name%.*}
screen=0
perl src/new_assess_180910.pl --dictionary data/dictionary.v10.studyID2.tsv --data example_pipeline --out result_$name --key $key
