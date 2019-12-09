package MyModule;
use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(GetName NameRecode CaseResult ReportCase ReportVariant ReadDictionary RecodeFileInDir File2Hash Belong Unique Intersect Diff);

########################################
# GetName
# -input: file path, column of sample ids, tags for skipping lines
# -output: unique sample names(ARRAY)
########################################
sub GetName{
	my ($file,$col,$skipline) = @_;
	my @sample_name0 = `cut -f $col $file | uniq | grep $skipline -v`;
	my @sample_name = map {substr($_,0,9)} @sample_name0;
	chomp @sample_name;
	return @sample_name;
}


########################################
# NameRecode
# -input: sample names(ARRAY)
# -output: old name to new name (HASH)
########################################
sub NameRecode{
	my @sample_name = @_;
	#old_name => new_name store in %new;
	my %new = ();
	my $n = $#sample_name;
	my %seen = ();
	while(keys %seen < $n+2){
		my $j = int (rand 10000)+10000;
		$seen{$j} =1;
	}
	
	my @new_name_id = keys %seen;
	foreach my $i (0..$n){
		$new{$sample_name[$i]}=$new_name_id[$i];
	}

	return %new;

}

########################################
#CaseResult
# -input: case by case file name from assessment tool e.g. manual_review_<pipeline>/1_all_cases_<pipeline>
# -output: \%case2tag \%case2dis,\%case2gene \%case2scad(HASH references)
########################################
sub CaseResult{
	my $file = shift @_;
	chomp $file;
	my %case2tag = ();
	my %case2dis = ();
	my %case2gene = ();
	my %case2scad = ();
	open(IN,"<$file") or die "can not open $file\n";
	while(my $line = <IN>){
		if ($line =~ "sample"){
			next;
		}
		my @aa = split /\t/,$line;
		my $sample = substr $aa[0],0,9;
		$case2tag{$sample} = $aa[1];
		$case2dis{$sample} = $aa[4];
		$case2gene{$sample} = $aa[6];
		$case2scad{$sample} = $aa[9];
	}
	close IN;
	return (\%case2tag,\%case2dis,\%case2gene,\%case2scad);

}

########################################
#CombineValue
# -input: reference of Hash 1, reference of Hash 2
# -output: Hash 3: key=>"var1\tvar2"
########################################
sub CombineValue{
	my ($ref1,$ref2) = @_;
	my %hash1 = %{$ref1};
	my %hash2 = %{$ref2};
	my %hash3 = ();
	foreach my $key (keys %hash1){
		if(!defined($hash1{$key})){
			next;
		}
		if(!defined($hash2{$key})){
                	next;
		}
		$hash3{$key} = "$hash1{$key}\t$hash2{$key}";
	
	}
	return %hash3;

}

#######################################
# ReportCase
######################################
sub ReportCase{
	my ($case,$outfile,$screen_file,$hash_ref) = @_;
	my %new = %{$hash_ref};
	my %case2fp_screen = &File2Hash($screen_file);
        my ($ref1,$ref2,$ref3,$ref4) = &CaseResult($case);
        my %case2tag = %{$ref1};
        my %case2dis = %{$ref2};
        my %case2gene = %{$ref3};
        my %case2scad = %{$ref4};
        my %case2scadraw = ();

#Apply SCAD clinical standard
        foreach my $sample (keys %case2dis){
                if(($case2dis{$sample} eq "77") and ($case2tag{$sample} eq "FN") and ($case2scad{$sample} eq "Y")){
                        $case2dis{$sample} = "0";
                        $case2tag{$sample} = "TN";
                        $case2scadraw{$sample} = "SCAD_FN";
                }
                else{
                        $case2scadraw{$sample} = "-";
                }
        }

        my %case2newtag1 = &CombineValue(\%new,\%case2tag);
        my %case2newtag = &CombineValue(\%case2newtag1,\%case2scadraw);
	

	open(OUT,">$outfile") or die "can not open $outfile\n";
	print OUT "sample\tcode\tSCADraw\tDisorder\tFlagged_Gene\tFP_screen\n";
	foreach my $sample (keys %case2dis){
		my $fp_screen = "unknown";
                if(defined $case2fp_screen{$sample}){
                	$fp_screen = $case2fp_screen{$sample};
                }
	my $newline = join ("\t",$case2newtag{$sample},$case2dis{$sample},$case2gene{$sample},$fp_screen);
	print OUT "$newline\n";

	}

}



#######################################
#ReportVariant
# -input:case result file, col, skipline, variant file, dictionary, output file
# -output: result file
######################################
sub ReportVariant{
	my($case,$col,$skipline,$var,$dictionary,$outfile,$screen_file,$hash_ref) = @_;
	my @name = &GetName($case,$col,$skipline);
	my %new = %{$hash_ref};
	my %case2fp_screen = &File2Hash($screen_file);
	my ($ref1,$ref2,$ref3,$ref4) = &CaseResult($case);
	my %case2tag = %{$ref1};
	my %case2dis = %{$ref2};
	my %case2gene = %{$ref3};
	my %case2scad = %{$ref4};
	my %case2scadraw = ();
	
#Apply SCAD clinical standard
	foreach my $sample (keys %case2dis){
		if(($case2dis{$sample} eq "77") and ($case2tag{$sample} eq "FN") and ($case2scad{$sample} eq "Y")){
			$case2dis{$sample} = "0";
			$case2tag{$sample} = "TN";
			$case2scadraw{$sample} = "SCAD_FN";
		}
		else{
			$case2scadraw{$sample} = "-";
		}
	}	

	my %disorder2gene = &ReadDictionary($dictionary);
	my %case2newtag1 = &CombineValue(\%new,\%case2tag);
	my %case2newtag = &CombineValue(\%case2newtag1,\%case2scadraw);

	open(IN,"<$var") or die "can not open $var\n";
	open(OUT,">$outfile") or die "can not open $outfile\n";
	while(my $line = <IN>){
		my @aa = split /\t/,$line;
		if($aa[1] eq "sample"){
			shift @aa;
			shift @aa;
			my $newheader = join("\t","sample","code","SCADraw","CorrectDiseaseGene","FlaggedByPipeline","disorder","FP_screen",@aa);
			print OUT $newheader;
			next;
		}
		my $sample = substr $aa[1],0,9;
		if($sample and $case2newtag{$sample}){
			
			shift @aa;
			shift @aa;
			my $gene = $aa[6];
			unless(defined $case2dis{$sample}){
				next;
			}
			my $disorder = $case2dis{$sample};
			my $fp_screen = "unknown";
			if(defined $case2fp_screen{$sample}){
				$fp_screen = $case2fp_screen{$sample};
			}
			my @candidate = ();
			if($disorder =~ /; /){
				my @d = split /; /,$disorder;
				foreach my $i (@d){
					push (@candidate, @{$disorder2gene{$i}});

				}
			}

			else{
				@candidate = @{$disorder2gene{$disorder}};
			}
			my $is_dis_gene = "incorrect";
			my $is_pipeline_gene = "not_flagged";
			if($gene eq $case2gene{$sample}){
				$is_pipeline_gene = "flagged";
			}
			elsif($case2gene{$sample} =~ /; /){
				my @multigene = split /; /,$case2gene{$sample};
				#print "test:$case2gene{$sample}\n";
				if(&Belong($gene,@multigene)){
					$is_pipeline_gene = "flagged";
				}
			}

			if(&Belong($gene,@candidate)){
				$is_dis_gene = "correct";
			}

			my $newline = join ("\t",$case2newtag{$sample},$is_dis_gene,$is_pipeline_gene,$disorder,$fp_screen,@aa);
			print OUT $newline;
                }

	}
	
	if (-e "var/xhmm.cnv.auto" ){
	open(CNV,"var/xhmm.cnv.auto");
	while(my $line = <CNV>){
		chomp $line;
		my @aa = split /\t/,$line;
		my $sample = substr $aa[0],0,9;
		unless(defined $case2gene{$sample}){
			next;
		}
		if($sample and $case2newtag{$sample}){
			shift @aa;
			my $gene = $aa[0];
			my $disorder = $case2dis{$sample};
			my $fp_screen = "unknown";
                        if(defined $case2fp_screen{$sample}){
                                $fp_screen = $case2fp_screen{$sample};
                        }
			my @candidate = ();
			if($disorder =~ /; /){
				my @d = split /; /,$disorder;
				foreach my $i (@d){
					push (@candidate, @{$disorder2gene{$i}});

				}
			}

			else{
				@candidate = @{$disorder2gene{$disorder}};
			}
			my $is_dis_gene = "incorrect";
			my $is_pipeline_gene = "not_flagged";
			if($gene eq $case2gene{$sample}){
				$is_pipeline_gene = "flagged";
			}
			if(&Belong($gene,@candidate)){
				$is_dis_gene = "correct";
			}

			my $newline = join ("\t",$case2newtag{$sample},$is_dis_gene,$is_pipeline_gene,$disorder,$fp_screen,"CNV","CNV","CNV","CNV","CNV","CNV",@aa);
			print OUT "$newline\n";
                }

	}

	close CNV;

	}

	close IN;
	close OUT;




}

######################################
# ReadDictionary
# Read in the dictionary
#
#
######################################
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
	my @bb = @{$disorder2gene{0}};
	return %disorder2gene;
}


####################################
# RecodeFileInDir
# change sample IDs of each files in the folder
# Input: dir, output dir, hash reference from File2Hash
#####################################
sub RecodeFileInDir{
        my ($dir,$outdir,$hash_ref) = @_;
        opendir(DIR,$dir);
        my %hash = %{$hash_ref};
        while(my $file = readdir(DIR)){
		if($file =~ /^\./){
			next;
		}

                my $out = "$outdir/recode_"."$file";

                open(OUT,">$out") or die "can not open $out\n";

                open(IN,"<$dir/$file") or die "can not open $dir/$file\n";
                while(my $line = <IN>){
                        my $newname = "";
                        my @aa = split "\t",$line;
                        chomp $aa[0];
                        my $name = substr($aa[0],0,9);
                        if (defined $hash{$name}){
                                $newname = $hash{$name};
                        }
                        else{
                                $newname = $aa[0];
                        }
                        shift @aa;

                        my $newline = join("\t",$newname,@aa);
                        print OUT $newline;

                }

        close IN;
        close OUT;
        }
}



###################################
# File2Hash
##################################
sub File2Hash{
        my $file = shift @_;
        my %hash = ();
        open(IN,"<$file") or die "can not open $file\n";
        while(my $line = <IN>){
                chomp $line;
                my @aa = split /\t/,$line;
		my $sample = substr($aa[0],0,9);
                $hash{$sample} = $aa[1];
        }
        return %hash;

}


##################################
# Check if elment a is in @A
#
#
##################################

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
1
