#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Time::localtime;
use File::Spec;
use File::Basename;
use Cwd;
my $VERSION = "v3.0.2";
#Yaqiong Wang
#yqw@berkeley.edu

my $usage = qq^
 assess.pl $VERSION

 usage: 
        perl assess.pl --dic dictionary.tsv --key key.tsv --data pipelines_dir --out result
        

Required:

 --dictionary <file>            dictionary file in tsv format, 
                                columns: DiseaseID, Screen, GeneList(comma separated)

 --key <file>                   answer key file in tsv format, 
                                columns: sample, DiseaseID, MS/MS results(required with --enable_screen_compare)

 --data <dir>                   directory which contains the pipeline prediction results

 --out <file>                   output file

Optional:

 --list <file>                  sample list
                                (if it is not provided, all samples appearing in both key file and prediction result file will be caculated)

^;

# option list:
my ($dictionary,$key,$datadir,$result0,$list,$fp_screen);
# defaults:
my $help_flag = 0;

unless(@ARGV){
        die "$usage\n";
        }
my @arg = ($0,@ARGV);
&GetOptions(
        'h|help' => \$help_flag,
        "dic|dictionary=s" => \$dictionary,
        "key=s" => \$key,
        "data=s" => \$datadir,
        "out=s" => \$result0,
        "list=s" =>\$list,
	"screen=s" => \$fp_screen,
);

my @_all_assess_params = qw(

dic
dictionary
key
data
list
screen
);

my %acceptable_opts = map {+ $_ => 1} @_all_assess_params;

my $opts_not_recognized_flag = 0;
for my $opt (@ARGV){
        if($opt =~ /^-+(\S+)/){
                my $opt_name = $1;
                unless ($acceptable_opts{$opt_name}){
                        print STDERR "ERROR, don't recognize parameter: $opt\n";
                        $opts_not_recognized_flag = 1;
                }
        }
}
if($opts_not_recognized_flag){
        die"Please review usage info for accepted parameters.\n"
}

if($help_flag){
        die "$usage\n";
        }

if(@ARGV){
        die"Error, do not understand options: @ARGV\n";
}


#output files
my $cwd = cwd();
my $result1 = $cwd."/$result0";

#unblinded results
my $result = $result1."_unblinded";


if(-d $result){
        die "Error:a folder called $result already exists\n";
}
elsif(-e $result){
        die "Error: a file named $result already exists\n";
}
else{
`mkdir $result`;
`mkdir $result/var_anno`;
`mkdir $result/sens`;
`mkdir $result/spec`;
`mkdir $result/MS_recode`;
}

#blind results
my $blind_result = $result1."_UCB_blinded";
my $blind_zip = $blind_result.".zip";
my $af = "$blind_result/AF";
if(-d $blind_result){
        die "Error:a folder called $blind_result already exists\n";
}
elsif(-e $blind_result){
        die "Error: a file named $blind_result already exists\n";
}
else{
`mkdir $blind_result`;
`mkdir $af`;

}


#output file names
my $sum_pipeline = $result."/summary_pipeline.tsv";
my $sum_disorder = $result."/summary_disorder.tsv";
my $sum_screen = $result."/summary_screen.tsv";
my $sum_gene = $result."/summary_gene.tsv";
my $logfile = $result."/assess_log.txt";
#my $casebycase_dir = $result."/manual_review";
#`mkdir $casebycase_dir`;

open(SUM,">$sum_pipeline") or die "can not open $sum_pipeline\n";
open(SUM_DIS,">$sum_disorder") or die "can not open $sum_disorder\n";

open(SUM_SCREEN,">$sum_screen") or die "can not open $sum_screen\n";
open(SUM_Gene,">$sum_gene") or die "can not open $sum_gene\n";
my $zip = "$result.zip";


#Log command line parameters
my $assess_arguments = "";
foreach (@arg){
        $assess_arguments = $assess_arguments . " " . $_;
}
open(LOG,">$logfile") or die "can not open $logfile\n";
my $now = ctime();
print LOG "$now\n";
print LOG "cwd:$cwd\n";
print LOG "perl $assess_arguments\n";

#######################################


my @dictionary_values = &ReadDictionary($dictionary);
# return in ReadDictionary:
# return (\$disorder2gene,\$disorder2screen,\$screen2disorder,\$gene2disorder);
my %disorder2gene = %{$dictionary_values[0]};
my %disorder2screen = %{$dictionary_values[1]};
my %screen2disorder = %{$dictionary_values[2]};
my %gene2disorder = %{$dictionary_values[3]};
my %gene2screen = %{$dictionary_values[4]};
my %screen2gene = %{$dictionary_values[5]}; 
my %id2name = %{$dictionary_values[6]};
my %sample2ID = &ReadKey($key);
my %sample2screen = &ConvertDisorder2Screen(\%sample2ID);
# map sample to screen



##################################
#check the answer key
foreach my $s (sort keys %sample2ID){		
	foreach my $id (@{$sample2ID{$s}}){
		my @valid_id = (keys %disorder2gene);
		if(!&Belong($id,@valid_id)){	
		die "Error:$s has invalid disorder ID $id\n";
		}
	}
}

################################
my @sample_subset=();
if($list){
        open(LIST,"<$list") or die "can not open $list\n";
        @sample_subset = <LIST>;
        chomp @sample_subset;
        foreach my $s (@sample_subset){
                if(!$sample2ID{$s}){
                        die "$list contains samples id $s which not in the key file: $key\n";
                }
        }

	close LIST;
}
else{
        @sample_subset = (sort (keys %sample2ID));
}


#################################
# Read pipeline results from pipeline dir
opendir(DIR,$datadir) or die "can not open dir: $datadir\n";

print SUM "Pipeline\tDisease.TP\tDisease.TN\tDisease.FP\tDisease.FN\tDisease.TP+FP\tDisease.TP+FN\tDisease.FN+FP\tScreen.TP\tScreen.TN\tScreen.FP\tScreen.FN\tScreen.TP+FP\tScreen.TP+FN\tScreen.FN+FP\tFNwHet\tFNwCNVHet\tFNw_Multi_HetGene\tFNFPwHet\tFNFPwCNVHet\tFNFPw_Multi_HetGene\n";
print SUM_DIS "Pipeline\tDisorderID\tDisorderName\tTP\tFN\tFN+FP\taffected\tFNwHet\tFNwCNVHet\tFNw_Multi_HetGene\tFNFPwHet\tFNFPwCNVHet\tFNFPw_Multi_HetGene\n";
print SUM_SCREEN "Pipeline\tScreen\tTP\tFN\tFN+FP\taffected\n";
print SUM_Gene "Pipeline\tGene\tDisorderID\tTP\tFN+FP\tFP\tFNwHet\tFNwCNVHet\tFNFPwHet\tFNFPwCNVHet\n";


while(my $file = readdir(DIR)){
	if($file =~ /^\./){next;}
	#test print
        my @assess_samples = @sample_subset;
	print LOG "$file\n";	
	&Assess("$datadir/$file",\@assess_samples);

}

close SUM_DIS;
close SUM_SCREEN;
close SUM_Gene;
close SUM;

#####################################

use lib 'src';
use MyModule qw(GetName NameRecode ReportCase ReportVariant RecodeFileInDir File2Hash);

if(-e "var"){
	opendir(VarDir,"var");

	my @sample_all = &GetName("data/all_sample",1,"sample");
	my %recode_hash = &NameRecode(@sample_all);

	while(my $file = readdir(VarDir)){
		if($file=~/(\S+)\.var$/){
			my $case = "$result/manual_review_$1/1_all_cases_$1.tsv";
			my $varoutfile = "$result/var_anno/$1.var";
			my $caseoutfile = "$result/var_anno/$1.case.tsv";
			my $var = "var/$file";
			if(-e $case){
        	        	&ReportVariant($case,1,"sample",$var,$dictionary,$varoutfile,$fp_screen,\%recode_hash);
				&ReportCase($case,$caseoutfile,$fp_screen,\%recode_hash);
			}
		}

	}
	closedir(VarDir);

	&RecodeFileInDir("data/MS","$result/MS_recode",\%recode_hash);

}
#####################################

##copy summarized files to blind folder
###
`cp $sum_pipeline $blind_result`;
`cp $sum_disorder $blind_result`;
`cp $sum_screen $blind_result`;
`cp $sum_gene $blind_result`;
`cp $logfile $blind_result`;
`cp -r $result/var_anno $blind_result`;
`cp -r $result/sens $blind_result`;
`cp -r $result/spec $blind_result`;
`cp -r $result/MS_recode $blind_result`;
# zip folders
#`cd $result/`;
`zip -r $zip $result/*`;
#`cd $blind_result/`;
`zip -r $blind_zip $blind_result/*`;



##################################
sub ConvertDisorder2Screen{
	my %mysample2ID = %{$_[0]};
	my %mysample2screen = ();
	foreach my $sample (sort keys %mysample2ID){
		foreach my $id (@{$mysample2ID{$sample}}){
			push (@{$mysample2screen{$sample}}, $disorder2screen{$id}[0]);
		}
	}
	return %mysample2screen;
}

##################################
# Read in the answer key file
sub ReadKey{
	my $mykey = shift @_;
	my %mysample2ID = ();
	open(IN, "$mykey") or die "can not open $key\n";
	while(my $line = <IN>){
		#skip the comment(#) lines and empty lines
		if(($line =~ /^#/) or ($line =~/^\n/) or ($line =~ /sample/)){
			next;
		}
		chomp $line;

		my @aa = split /\t/,$line;
		my $sample = $aa[0];

		#clean the sample names
		$sample =~ s/\s+//g;
		$sample =substr $sample,0,9;
		
		my @IDs = ();
		# single disorder ID
		$aa[1] =~ s/\s+//g;
		if($aa[1] =~ /,/ ){
			@IDs = split /,/,$aa[1];
			
		}
		else{
			$IDs[0] = $aa[1];
		}
		
		#check if there are duplicated sample names
		if(defined($mysample2ID{$sample})){
			die "Error: duplicated sample names in answer key!\n";
		}
		
		$mysample2ID{$sample} = \@IDs;
	}
	close IN;
	return %mysample2ID;
	
}
##################################
# Read in the dictionary
sub ReadDictionary{
	my $dic = shift @_;
	
	my %gene2disorder = ();
	my %disorder2gene = ();
	my %disorder2screen = ();
	my %screen2disorder = ();
	my %gene2screen = ();
	my %screen2gene = ();
	my %id2name = ();	
	open(IN,"<$dic") or die "can not open $dic\n";
	my $n=1;
	while(my $line = <IN>){
		#skip the comment(#) lines,empty lines and header lines
                if(($line =~ /^#/) or ($line =~/^\n/) or ($line =~ /^DiseaseID/)){
                        next;
                }
                chomp $line;
		my @aa = split /\t/,$line;
                my $disorder = $aa[0];
		my $name = $aa[1];
		$id2name{$disorder} = $name;
		my $screen = $aa[2];
		if(! $aa[3]){
			next;
		}
		$aa[3] =~ s/\s+//g;
		#$n++;
		#print "$n:$aa[3]\n";

		my @bb = ();
		#multiple genes
		if($aa[3] =~ /,/){
			@bb = split /,/,$aa[3];
		}
		else{
		#single gene
			$bb[0] = $aa[3];
		}

		$disorder2gene{$disorder} = [@bb];
	
		#if screen exists
		if($screen){
			push (@{$disorder2screen{$disorder}},$screen);
			push (@{$screen2disorder{$screen}},$disorder);
		}

		#gene to disorder map
		foreach my $gene (@bb){
			push (@{$gene2disorder{$gene}},$disorder);
			if($screen){
				push (@{$gene2screen{$gene}},$screen);
				push (@{$screen2gene{$screen}},$gene);
			}
		} 

	}

	foreach my $gene (keys %gene2screen){
		@{$gene2screen{$gene}} = &Unique(@{$gene2screen{$gene}});
		@{$gene2disorder{$gene}} = &Unique(@{$gene2disorder{$gene}});
	}
	foreach my $screen (keys %screen2gene){
		@{$screen2gene{$screen}} = &Unique(@{$screen2gene{$screen}});
		@{$screen2disorder{$screen}} = &Unique(@{$screen2disorder{$screen}});
	}
	foreach my $disorder (keys %disorder2gene){
		@{$disorder2gene{$disorder}} = &Unique(@{$disorder2gene{$disorder}});
		@{$disorder2screen{$disorder}} = &Unique(@{$disorder2screen{$disorder}});
	}
	close IN;
	#return four hashes that store the disorder, screen and gene relationships.
	return (\%disorder2gene,\%disorder2screen,\%screen2disorder,\%gene2disorder,\%gene2screen,\%screen2gene,\%id2name);
}
#################################
# Assessment on one pipeline
sub Assess{
	my $file = $_[0];
	my @assess_samples = @{$_[1]};
	my $pipeline = fileparse($file); 
	$pipeline =~ s/\.tsv//g;
	my $casebycase_dir = $result."/manual_review_".$pipeline;
	`mkdir $casebycase_dir`;
	my $pipeline_out = "$casebycase_dir/"."1_all_cases_".$pipeline.".tsv";
	my $fp_out = "$casebycase_dir/"."2_exome_false_positive_cases_".$pipeline.".tsv";
	my $fn_out = "$casebycase_dir/"."3_exome_false_negative_cases_".$pipeline.".tsv";
	my $fnfp_out = "$casebycase_dir/"."4_exome_wrong_gene_affected_cases_".$pipeline.".tsv";
	my $other_out = "$casebycase_dir/"."5_other_disorder_".$pipeline.".tsv";
	my $sens_out = "$result/"."sens/sens_".$pipeline.".tsv";
	my $spec_out = "$result/"."spec/spec_".$pipeline.".tsv";
	my $fp_af = "$af/"."af_FP_".$pipeline.".tsv";
	my $tp_af = "$af/"."af_TP_".$pipeline.".tsv";
	my $fnfp_af = "$af/"."af_FNFP_".$pipeline.".tsv";
	my $fnhet = "$af/"."FN_MultiHetGene_".$pipeline.".txt";
	my $fnfphet = "$af/"."FNFP_MultiHetGene_".$pipeline.".txt";	
	open (CASE,">$pipeline_out") or die "can not open $pipeline_out\n";
	open (FP,">$fp_out") or die "can not open $fp_out\n";
	open (FN,">$fn_out") or die "can not open $fn_out\n";
	open (FNFP,">$fnfp_out") or die "can not open $fnfp_out\n";
	open (OTHERD,">$other_out") or die "can not open $other_out\n";
	open (FP_AF,">$fp_af") or die "can not open $fp_af\n";
	open (TP_AF,">$tp_af") or die "can not open $fp_af\n";
	open (FNFP_AF,">$fnfp_af") or die "can not open $fp_af\n";
	open (FNHet,">$fnhet") or die "can not open $fnhet\n";
	open (FNFPHet,">$fnfphet") or die "can not open $fnfphet\n";	
	my $headerline = "sample\tDisorder\tScreen\tState\tKey_Disorder_ID\tKey_Disorder_Name\tPredicted_Gene\tPredicted_Disorder_ID\tPredict_Disorder_Name\tSCAD_variant_c.625G>A\tNum_of_SCAD_rare_variant\tKey_Disorder_ID_Correct?\tConfidence[0-1]\tDelete_Sample?\tChange ID\n";

	print CASE $headerline; 
	print FP $headerline; 
	my $headerline_FN = "sample\tDisorder\tScreen\tState\tKey_Disorder_ID\tKey_Disorder_Name\tPredicted_Gene\tPredicted_Disorder_ID\tPredict_Disorder_Name\tSCAD_variant_c.625G>A\tNum_of_SCAD_rare_variant\tHetGene_Correct\tHetGene\tCNV_HetGene_Correct\tCNV_HetGene\tMultiple_HetGene_Correct\tKey_Disorder_ID_Correct?\tConfidence[0-1]\tDelete_Sample?\tChange ID\n"; 
	print FN $headerline_FN;
	print FNFP $headerline_FN; 
	print OTHERD $headerline;

	#return values from ReadPipelineResult \%sample2gene,\%sample2SCAD,\%sample2variant,\%sample2variant_flag,\%sample2af,\%sample2Het,\%sample2CNVHet; 
	my ($ref_sample2gene,$ref_sample2SCAD,$ref_sample2variant,
		$ref_sample2variant_flag, $ref_sample2af, $ref_sample2Het, $ref_sample2CNVHet) = &ReadPipelineResult($file);
	my %sample2gene = %{$ref_sample2gene};
	my %sample2SCAD = %{$ref_sample2SCAD};
	my %sample2variant = %{$ref_sample2variant};
	my %sample2variant_flag = %{$ref_sample2variant_flag};
	my %sample2af = %{$ref_sample2af};
	my %sample2Het = %{$ref_sample2Het};
	my %sample2CNVHet = %{$ref_sample2CNVHet};
	my %d_sum = ();
	my %s_sum = ();
	my %het_sum = ();
	my %cnvhet_sum = ();
	my %d_hash = ();
	my %g_hash = ();
	my %s_hash = ();
	my %het_hash = ();
	my %het_hash_gene = ();
	my %het_cnvhash = ();
	my %het_cnvhash_gene = ();
	my $hetgenecount = 0;
	my %MultiHetGene = ();
	my %fnfphet_sum = ();
	my %fnfpcnvhet_sum = ();
	my %fnfphet_hash = ();
	my %fnfphet_hash_gene = ();
	my %fnfphet_cnvhash = ();
	my %fnfphet_cnvhash_gene = ();
	my $fnfphetgenecount = 0;
	my %fnfpMultiHetGene = ();

	my %affected = ();
	my %affected_screen = ();
	my @sample2 = keys %sample2gene;
	my @subset_samples = &Intersect(\@assess_samples,\@sample2);
	my $other_sum = 0;
	foreach my $sample (@assess_samples){
		if(!&Belong($sample,@sample2)){
		print LOG "$sample\tNot pipeline predicted sample\n";
		}
	}

	foreach my $sample (@subset_samples){
		#count affected samples
		my @disnames = ();
		foreach my $dis (@{$sample2ID{$sample}}){
			$affected{$dis}++;
			push(@disnames,$id2name{$dis});
		}

		foreach my $screen (@{$sample2screen{$sample}}){
			$affected_screen{$screen}++;
		}
		my @predict_dis_id = ();	
		my @predict_dis_name = ();
		foreach my $gene(@{$sample2gene{$sample}}){
			foreach my $id (@{$gene2disorder{$gene}}){
			push(@predict_dis_name,$id2name{$id});
				}
			my $dis_id = join(",",@{$gene2disorder{$gene}});
			push (@predict_dis_id,$dis_id);
		
		}
		#999 other disorders
		my ($tag,$d_ref,$g_ref) = ("NA","","");
		my ($stag,$s_ref,$sg_ref) = ("NA","",""); 		
		#non-999
		if(${$sample2ID{$sample}}[0] ne "999"){
			($tag,$d_ref,$g_ref) = &CompareDisorder($sample2ID{$sample},$sample2gene{$sample},\%gene2disorder,\%disorder2gene);
			($stag,$s_ref,$sg_ref) = &CompareDisorder($sample2screen{$sample},$sample2gene{$sample},\%gene2screen,\%screen2gene);
		}
		my $diseaselist = join("; ",@{$sample2ID{$sample}});
		my $screenlist = join("; ",@{$sample2screen{$sample}});
		my $genelist = join("; ",@{$sample2gene{$sample}});
		my $diseasenamelist = join("; ",@disnames);
		my $predict_dis_id_list = join("; ",@predict_dis_id);
		my $predict_dis_name_list = join("; ",@predict_dis_name);
		my $state = $tag;
		my $scadvariant = $sample2SCAD{$sample};
		my $scad_rare_variant = $sample2variant{$sample};
		if($tag eq "FN+FP"){
			$state = "TP";
		}
		elsif($tag eq "TP+FP"){
			$state = "TP";
		}
		elsif($tag eq "TP+FN"){
			$state = "TP";
		}
		
		my $case_printout_line = "$sample\t$tag\t$stag\t$state\t$diseaselist\t$diseasenamelist\t$genelist\t$predict_dis_id_list\t$predict_dis_name_list\t$scadvariant\t$scad_rare_variant\n";
		print CASE $case_printout_line; 
		$d_sum{$tag}++;
		$s_sum{$stag}++;

		#print out FN cases
		if(($tag eq "FN") or ($tag eq "TP+FN")){
			chomp $case_printout_line;
			print FN $case_printout_line;
			
			#### Het genes
			my ($het_tag,$ref_het,$ref_het_gene) = &CompareDisorder($sample2ID{$sample},$sample2Het{$sample},\%gene2disorder,\%disorder2gene);
			my ($cnvhet_tag,$ref_cnvhet,$ref_cnvhet_gene) = &CompareDisorder($sample2ID{$sample},$sample2CNVHet{$sample},\%gene2disorder,\%disorder2gene);
			
			### case counts
			$het_sum{$het_tag}++;
			$cnvhet_sum{$cnvhet_tag}++;
			if(($het_tag eq "TP") | ($het_tag eq "TP+FP")){
				$het_tag = "Y";
			}
			else{$het_tag = "N";}
			
			if(($cnvhet_tag eq "TP") | ($cnvhet_tag eq "TP+FP")){
                                $cnvhet_tag = "Y";
                        }
                        else{$cnvhet_tag = "N";}
			my $hetgene = join(",",@{$sample2Het{$sample}});
			my $cnvhetgene = join(",",@{$sample2CNVHet{$sample}});
			print FN "\t$het_tag\t$hetgene\t$cnvhet_tag\t$cnvhetgene";

			#mult Het genes
			my($mult_tag,$ref_hetgenelist) = &HetGeneCount($ref_het_gene,$ref_cnvhet_gene);
			
			print FN "\t$mult_tag\n";
			if($mult_tag eq "Y"){
				my $list = join (",",@{$ref_hetgenelist});
				print FNHet "$list\n";
				$hetgenecount++;

				
			}
			### gene counts begin ###
			foreach my $d (sort keys %{$ref_het}){
                                        my $result = $ref_het->{$d};
                                        ${$het_hash{$d}}{$result}++;
					if((($result eq "TP") or ($ref_cnvhet->{$d} eq "TP")) and ($mult_tag eq "Y")){
					$MultiHetGene{$d}++;
					}
                                }

			foreach my $g (sort keys %{$ref_het_gene}){
                                        my $result = $ref_het_gene->{$g};
                                        ${$het_hash_gene{$g}}{$result}++;
                                }

			foreach my $d (sort keys %{$ref_cnvhet}){
                                        my $result = $ref_cnvhet->{$d};
                                        ${$het_cnvhash{$d}}{$result}++;
                                }
			foreach my $g (sort keys %{$ref_cnvhet_gene}){
                                        my $result = $ref_cnvhet_gene->{$g};
                                        ${$het_cnvhash_gene{$g}}{$result}++;
                                }
				
			
		}

		#print out FP cases
		if(($tag eq "FP") or ($tag eq "TP+FP")){
                        print FP $case_printout_line;
			if($sample2variant_flag{$sample} ne "NoVariants"){
				print FP_AF "$sample2variant_flag{$sample}\t$sample2af{$sample}\n"; 
                	}
		}
		#TP cases
		if($tag eq "TP"){
			if($sample2variant_flag{$sample} ne "NoVariants"){
				print TP_AF "$sample2variant_flag{$sample}\t$sample2af{$sample}\n";
			}

		}
		#print out FN+FP cases
		if($tag eq "FN+FP"){
			chomp $case_printout_line;
			print FNFP $case_printout_line;
			# FP variants
			if($sample2variant_flag{$sample} ne "NoVariants"){
				print FNFP_AF "$sample2variant_flag{$sample}\t$sample2af{$sample}\n";
			}

			###FN Het variants
			my ($het_tag,$ref_het,$ref_het_gene) = &CompareDisorder($sample2ID{$sample},$sample2Het{$sample},\%gene2disorder,\%disorder2gene);
			my ($cnvhet_tag,$ref_cnvhet,$ref_cnvhet_gene) = &CompareDisorder($sample2ID{$sample},$sample2CNVHet{$sample},\%gene2disorder,\%disorder2gene);
			
			### case counts
			$fnfphet_sum{$het_tag}++;
			$fnfpcnvhet_sum{$cnvhet_tag}++;
			if(($het_tag eq "TP") | ($het_tag eq "TP+FP")){
				$het_tag = "Y";
			}
			else{$het_tag = "N";}
			
			if(($cnvhet_tag eq "TP") | ($cnvhet_tag eq "TP+FP")){
                                $cnvhet_tag = "Y";
                        }
                        else{$cnvhet_tag = "N";}
			my $hetgene = join(",",@{$sample2Het{$sample}});
			my $cnvhetgene = join(",",@{$sample2CNVHet{$sample}});
			print FNFP "\t$het_tag\t$hetgene\t$cnvhet_tag\t$cnvhetgene";

			#mult Het genes
			my($mult_tag,$ref_hetgenelist) = &HetGeneCount($ref_het_gene,$ref_cnvhet_gene);
			#test print
			#print $mult_tag;
			print FNFP "\t$mult_tag\n";
			if($mult_tag eq "Y"){
				my $list = join (",",@{$ref_hetgenelist});
				print FNFPHet "$list\n";
				$fnfphetgenecount++;

				
			}
			### gene counts begin ###
			foreach my $d (sort keys %{$ref_het}){
                                        my $result = $ref_het->{$d};
                                        ${$fnfphet_hash{$d}}{$result}++;
					if((($result eq "TP") or ($ref_cnvhet->{$d} eq "TP")) and ($mult_tag eq "Y")){
					$fnfpMultiHetGene{$d}++;
					}
                                }

			foreach my $g (sort keys %{$ref_het_gene}){
                                        my $result = $ref_het_gene->{$g};
                                        ${$fnfphet_hash_gene{$g}}{$result}++;
                                }

			foreach my $d (sort keys %{$ref_cnvhet}){
                                        my $result = $ref_cnvhet->{$d};
                                        ${$fnfphet_cnvhash{$d}}{$result}++;
                                }
			foreach my $g (sort keys %{$ref_cnvhet_gene}){
                                        my $result = $ref_cnvhet_gene->{$g};
                                        ${$fnfphet_cnvhash_gene{$g}}{$result}++;
                                }
			

		}

		#print out "other disorder" cases
		if($tag eq "NA"){
			$other_sum++;
			print OTHERD $case_printout_line;
		}
			#summarize non-999 cases
			if($tag ne "NA"){
				foreach my $d (sort keys %{$d_ref}){
					my $result = $d_ref->{$d};
					${$d_hash{$d}}{$result}++;
				}
				foreach my $s (sort keys %{$s_ref}){
					my $result = $s_ref->{$s};
					${$s_hash{$s}}{$result}++;
				}
				foreach my $g (sort keys %{$g_ref}){
                                	my $result = $g_ref->{$g};
                                	${$g_hash{$g}}{$result}++;
                        	}

			}
	}
	print LOG "Number of other disorder cases: $other_sum\n";
	close CASE;
	close FN;
	close FP;
	close FNFP;
	close OTHERD;
	close FP_AF;
	close TP_AF;
	close FNFP_AF;
	close FNHet;
	close FNFPHet;
	
	### calculate matrix:
	`perl src/caculate_171031.pl $pipeline_out data/prevalence.txt data/disorder_group.txt $sens_out $spec_out`;


	#### print pipeline summary

	# on disorder level
	print SUM "$pipeline";
	foreach my $i ("TP","TN","FP","FN","TP+FP","TP+FN","FN+FP"){
		my $value = 0;
		if(defined $d_sum{$i}){
			$value = $d_sum{$i};
		}
		print SUM "\t$value";
	}	
	# on screen level
	foreach my $i ("TP","TN","FP","FN","TP+FP","TP+FN","FN+FP"){
                my $value = 0;
                if(defined $s_sum{$i}){
                        $value = $s_sum{$i};
                }
                print SUM "\t$value";
        }

	# on affected state level
	# on het level

	if(!defined $het_sum{"TP"}){ $het_sum{"TP"} = 0;}
	if(!defined $het_sum{"TP+FP"}){ $het_sum{"TP+FP"} = 0;}
	if(!defined $cnvhet_sum{"TP"}){ $cnvhet_sum{"TP"} = 0;}
	if(!defined $cnvhet_sum{"TP+FP"}){ $cnvhet_sum{"TP+FP"} = 0;}
	my $hetvalue = $het_sum{"TP"} + $het_sum{"TP+FP"};
	my $cnvhetvalue = $cnvhet_sum{"TP"} + $cnvhet_sum{"TP+FP"};	
	
	if(!defined $fnfphet_sum{"TP"}){ $fnfphet_sum{"TP"} = 0;}
        if(!defined $fnfphet_sum{"TP+FP"}){ $fnfphet_sum{"TP+FP"} = 0;}
        if(!defined $fnfpcnvhet_sum{"TP"}){ $fnfpcnvhet_sum{"TP"} = 0;}
        if(!defined $fnfpcnvhet_sum{"TP+FP"}){ $fnfpcnvhet_sum{"TP+FP"} = 0;}
        my $fnfphetvalue = $fnfphet_sum{"TP"} + $fnfphet_sum{"TP+FP"};
        my $fnfpcnvhetvalue = $fnfpcnvhet_sum{"TP"} + $fnfpcnvhet_sum{"TP+FP"};

	print SUM "\t$hetvalue\t$cnvhetvalue\t$hetgenecount\t$fnfphetvalue\t$fnfpcnvhetvalue\t$fnfphetgenecount";       
	print SUM "\n";

	#### print disorder summary
	foreach my $d (sort keys %d_hash){
		print SUM_DIS "$pipeline\t$d\t$id2name{$d}";
		foreach my $i ("TP","FN","FN+FP"){
			my $value = 0;
			if($d_hash{$d}->{$i}){
				$value = $d_hash{$d}->{$i};
			}
			print SUM_DIS "\t$value";
		}
		print SUM_DIS "\t$affected{$d}";

		###het counts
		foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $het_hash{$d}->{$i}){
                                $value = $het_hash{$d}->{$i};
                        }
                        print SUM_DIS "\t$value";
                }
		foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $het_cnvhash{$d}->{$i}){
                                $value = $het_cnvhash{$d}->{$i};
                        }
			#multi-Hets
                        print SUM_DIS "\t$value";
			if(!defined $MultiHetGene{$d}){
				$MultiHetGene{$d} = 0;
			}
			print SUM_DIS "\t$MultiHetGene{$d}";
                }

		#FNFP cases
		foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $fnfphet_hash{$d}->{$i}){
                                $value = $fnfphet_hash{$d}->{$i};
                        }
                        print SUM_DIS "\t$value";
                }
		foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $fnfphet_cnvhash{$d}->{$i}){
                                $value = $fnfphet_cnvhash{$d}->{$i};
                        }
			#multi-Hets
                        print SUM_DIS "\t$value";
			if(!defined $fnfpMultiHetGene{$d}){
				$fnfpMultiHetGene{$d} = 0;
			}
			print SUM_DIS "\t$fnfpMultiHetGene{$d}";
                }

		print SUM_DIS "\n";

	
	}
	### print screen summary
	foreach my $s (sort keys %s_hash){
                print SUM_SCREEN "$pipeline\t$s\t";
                foreach my $i ("TP","FN","FN+FP"){
                        my $value = 0;
                        if($s_hash{$s}->{$i}){
                                $value = $s_hash{$s}->{$i};
                        }
                        print SUM_SCREEN "$value\t";
                }
                print SUM_SCREEN "$affected_screen{$s}\n";
        }

	#### print gene summary
        foreach my $g (sort keys %gene2disorder){
                my @dis_array = @{$gene2disorder{$g}};
		my $dis_list = join (",",@dis_array);
		print SUM_Gene "$pipeline\t$g\t$dis_list";
		
                foreach my $i ("TP","FN+FP","FP"){
                        my $value = 0;
                        if(defined $g_hash{$g}->{$i}){
                                $value = $g_hash{$g}->{$i};
                        }
                        print SUM_Gene "\t$value";
                }
        
		###het counts 
                foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $het_hash_gene{$g}->{$i}){
                                $value = $het_hash_gene{$g}->{$i};
                        }
                        print SUM_Gene "\t$value";
                }       
                foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $het_cnvhash_gene{$g}->{$i}){
                                $value = $het_cnvhash_gene{$g}->{$i};
                        }
                        print SUM_Gene "\t$value";
                }
	###het counts for FNFP 
                foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $fnfphet_hash_gene{$g}->{$i}){
                                $value = $fnfphet_hash_gene{$g}->{$i};
                        }
                        print SUM_Gene "\t$value";
                }       
                foreach my $i ("TP"){
                        my $value = 0;
                        if(defined $fnfphet_cnvhash_gene{$g}->{$i}){
                                $value = $fnfphet_cnvhash_gene{$g}->{$i};
                        }
                        print SUM_Gene "\t$value";
                }

               
                print SUM_Gene "\n";


	}


}
##################################
	




##################################
sub ReadPipelineResult{
        my $in = $_[0];
        my %sample2gene = ();
        my %sample2SCAD = ();
	my %sample2variant = ();
	my %sample2variant_flag = ();
	my %sample2af = ();
	my %sample2Het = ();
	my %sample2CNVHet = ();
	open(PIN,"<$in") or die "can not open $in\n";
	while(my $line = <PIN>){
		if($line =~ /^\n/){
			next;
		}
		chomp $line;
		my @aa = split /\t/,$line;
		$aa[0] = substr $aa[0],0,9;
		if($aa[1] eq "NoGene"){
			$sample2gene{$aa[0]}[0] = "None";
		}
		elsif($aa[1] =~ /,/){
			my @bb = split /,/,$aa[1];
			$sample2gene{$aa[0]} = \@bb;
		}
		else{
			$sample2gene{$aa[0]}[0] = $aa[1];
		}

		if(!$aa[2]){
			$aa[2]="";
		}

		$sample2SCAD{$aa[0]} = $aa[2];
		if(!$aa[3]){
			$aa[3]=0;
		}
		$sample2variant{$aa[0]} = $aa[3];
		$sample2variant_flag{$aa[0]} = $aa[4];
		$sample2af{$aa[0]} = $aa[5];
		$sample2Het{$aa[0]} = &String2Array($aa[6], ",", "NoHetGene");

		$sample2CNVHet{$aa[0]} = &String2Array($aa[7], ",", "NoHetCNVGene");
		
	}
	close PIN;
	return (\%sample2gene,\%sample2SCAD,\%sample2variant,\%sample2variant_flag,\%sample2af,\%sample2Het,\%sample2CNVHet);

}
###################################
sub String2Array{
        my ($string,$delm,$null) = @_; 
        my @bb;
	if(defined $string){
		if($string eq $null){
                	$bb[0] = "None";
        	}

		#remove space
        	elsif($string =~ "$delm"){
                	@bb = split "$delm",$string;
        	}
        	else{
                	$bb[0] = $string;
        	}
	}
	else{
		$bb[0] = "None";
		}
	return \@bb;
}

##################################
sub CompareDisorder{
	my @disease = @{$_[0]};
	my @gene = @{$_[1]};
	my %mygene2disorder = %{$_[2]};
	my %mydisorder2gene = %{$_[3]};
	my $tag = "";
	my %hash;
	my %ghash;
	#Unaffect individuals
	if(($disease[0] eq 0) or ($disease[0] eq "None")){
		if($gene[0] eq "None"){
			$tag = "TN";
		}
		else{	
			my $mark = 0;
			foreach my $g (@gene){
				if(&Belong($g,(keys %mygene2disorder))){
					$tag = "FP";
					$mark = 1;
					$ghash{$g} = "FP";	
				}
		
			}
	
			if ($mark == 0){
				$tag = "TN";
			}
			
		}

	}

        #Affected individuals
		#Gene call is None
	elsif($gene[0] eq "None"){
		$tag = "FN";
		foreach my $d (@disease){
			$hash{$d} = "FN";
		}	
	}
	
		#Has disease and gene calls
	else{
		my ($TP,$FP,$FN) = (0,0,0);
		#check every disorder
		foreach my $d (@disease){
			my $dmark = 0;
			foreach my $g (@gene){
				if(&Belong($d,@{$mygene2disorder{$g}})){
					$TP = 1;
					$dmark = 1;
					$hash{$d} = "TP";
					next;
				}
			}
			if($dmark == 0){
					$FN = 1;
					$hash{$d} = "FN+FP";
			}

		}
		#check every gene
		foreach my $g (@gene){
			my $gmark = 0;
			foreach my $d (@disease){
				if(&Belong($g,@{$mydisorder2gene{$d}})){
					$gmark = 1;
					$ghash{$g} = "TP";
				}
			}
			if($gmark == 0){
				$FP = 1;
				$ghash{$g} = "FN+FP";
			}

		}
		
		if($TP == 1 and $FN == 0 and $FP == 0){
			$tag = "TP";
		}
		if($TP == 1 and $FN == 1 and $FP == 0){
			$tag = "TP+FN";
		}
		if($TP == 1 and $FN == 0 and $FP == 1){
			$tag = "TP+FP";
		}
		if($TP == 1 and $FN == 1 and $FP == 1){
			$tag = "TP+FP+FN";
		}
		if($TP == 0 and $FN == 1 and $FP == 1){
			$tag = "FN+FP";
		}
		if($TP == 0 and $FN == 1 and $FP == 0){
			$tag = "FN";
		}


		
	}
	
	if($tag eq "TP+FN"){
		foreach my $i (keys %hash){
			if($hash{$i} eq "FN+FP"){
				$hash{$i} = "FN";
		}
		}
	}
	
	return($tag,\%hash,\%ghash);
}


##################################
#
#
sub HetGeneCount{
	my ($ref,$cnv_ref) = @_;
	my %hash = %{$ref};
	my %cnvhash = %{$cnv_ref};
	my @aa = ();
	my $tag = "N";
	foreach my $i (keys %hash){
		if($hash{$i} eq "TP"){
			#print "$i\t";
			push (@aa,$i);
		}
	}
	foreach my $i (keys %cnvhash){
		if($cnvhash{$i} eq "TP"){
			#print "$i\t";
			push (@aa,$i);
		}
	}

	my @bb = &Unique(@aa);
	if($#bb >= 1){
		$tag = "Y";
	}
	
	return($tag,\@bb);

}

##################################
# Export sample list
# disorder_list: array reference
# sample2ID: hash reference
# filename for exporting the sample list: should be under the sample path
sub RemoveDisorder{

	my @disorder_list = @{$_[0]};
	my %my2sample2ID = %{$_[1]};
	my $sample_list_output = $_[2];

	foreach my $sample (sort keys %my2sample2ID){
		foreach my $id (@{$my2sample2ID{$sample}}){
			if(&Belong($id,@disorder_list)){
				if($my2sample2ID{$sample}){
					delete $my2sample2ID{$sample};
				}
			}
		}
	}
	
	
	open(OUT,">$result/$sample_list_output") or die "can not open $result/$sample_list_output\n";
	print OUT "#remove samples has disorder:";
	foreach my $id (@disorder_list){
		print OUT "$id,";
		}
	print OUT "\n";

	foreach my $sample (sort keys %my2sample2ID){
		print OUT "$sample\n";
	}

	close OUT;
	return %my2sample2ID;
}

##################################
# Gene mapped to multiple disorders
# Get the gene and disorder IDs 
# argument: %gene2disorder
# return value: (\@mul_disorder,\@mul_gene)
sub GeneMappedMultipleDisorder{
	my @mul_disorder;
	my @mul_gene;
	my %gene2disorder = %{$_[0]};
	foreach my $gene (keys %gene2disorder){
		if($#{$gene2disorder{$gene}}>0){
			push (@mul_disorder,@{$gene2disorder{$gene}});
			push (@mul_gene,$gene);
			}
		
	}
	@mul_disorder = &Unique(@mul_disorder);
	@mul_gene = &Unique(@mul_gene);
	return (\@mul_disorder,\@mul_gene);
}
##################################
# Check if elment a is in @A
sub Belong{
        my ($gene,@list) = @_;
        my $value = 0;
        foreach my $i (@list){
		
                if ($i eq $gene){
                        $value = 1;
                        last;
                }
        }
        return $value;
}

##################################
# find unique elements in @A
# arguments: @A
# return: unique @A
sub Unique{
        my %seen = ();
        my @list = @_;
        foreach my $i (@list){
                $seen{$i}++;
        }
        my @uniq = sort keys %seen;
        return @uniq;
}

##################################
# find elements in both @A and @B
# arguments: \@A,\@B    
# return: intersection of @A and @B
sub Intersect{
	my @input = @_;
	my @aa = @{$input[0]};
	my @bb = @{$input[1]};
	@aa = &Unique(@aa);
	@bb = &Unique(@bb);
	my @both = ();
	foreach my $i (@aa){
		if(&Belong($i,@bb)){
			push(@both,$i);
		}
	}
	return @both;
}

##################################
# find elements in @A but not in @B
# arguments: \@A,\@B
# return: @A-@B
sub Diff{
	my @input = @_;
        my @aa = @{$input[0]};
        my @bb = @{$input[1]};
        @aa = &Unique(@aa);
        @bb = &Unique(@bb);
	my @diff = ();
	foreach my $i (@aa){
		unless(&Belong($i,@bb)){
			push(@diff,$i);
		}
	}
	return @diff;
}

exit;
