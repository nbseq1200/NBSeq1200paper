#!usr/bin/perl -w
use strict;
use MyModule qw(Unique Belong Intersect Diff);

my ($dic,$dis_sum,$pipeline_dir,$varfilepath,$out_dir) = @ARGV;

my $usage = "spec.pl <dic> <dis_sum> <pipeline_dir> <varfilepath><out_dir>\n";

if($#ARGV != 4) {die "$usage";}

my ($dis_ref,$gene_ref) = &ReadDictionary($dic);
my %gene2unrelated = &Gene2UnrelatedDis($dis_ref,$gene_ref);
my %gene2dis = %{$gene_ref};
my %dis2gene = %{$dis_ref};
delete $dis2gene{"0"};
delete $dis2gene{"999"};
delete $gene2dis{"999"};
delete $gene2dis{"None"};
delete $gene2unrelated{"999"};
delete $gene2unrelated{"None"};

foreach my $gene (keys %gene2dis){
	print "$gene2dis{gene}\n";
}


opendir(DIR,$pipeline_dir);
while(my $file = readdir(DIR)){
if($file =~ /^\./){next;}

$file =~ /(\S+?)\.tsv/;
my $pipeline = $1;
my $spec_out_file = "$out_dir/$pipeline.spec.tsv";

my $varfile = "$varfilepath/$pipeline.var";
my %genecount = &VarGeneCount($varfile);
my %disordercount = &DisorderCount($dis_sum,$pipeline);

my %spec = ();

my @genelist = (keys %gene2dis);
my $overall_fp = &ValueSum(\@genelist,\%genecount);
my $overall_assay = 0;
foreach my $gene (keys %gene2unrelated){
	
		$overall_assay = $overall_assay + &ValueSum($gene2unrelated{$gene},\%disordercount);
}

if ($overall_assay == 0){
	next;
}
my $n = @genelist;
my $total_spec = 1-$n*$overall_fp/$overall_assay;

print "caculating specificity: $pipeline ...\n";
open (OUT,">$spec_out_file") or die "can not open $spec_out_file\n";

print OUT "Pipeline\tDisorder\tNumOfGenes\tFPs\tAssays\tSpecificity\n";
print OUT "$pipeline\toverall\t$n\t$overall_fp\t$overall_assay\t$total_spec\n";


foreach my $dis (sort keys %dis2gene){
	my $fp = &ValueSum($dis2gene{$dis},\%genecount);
	my $assay = 0;
	foreach my $gene (@{$dis2gene{$dis}}){
		$assay = $assay + &ValueSum($gene2unrelated{$gene},\%disordercount);
	}		
	my $n = @{$dis2gene{$dis}};	
	$spec{$dis} = 1-$n*$fp/$assay;
	print OUT "$pipeline\t$dis\t$n\t$fp\t$assay\t$spec{$dis}\n";
}

close OUT;

}



###################
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
        return (\%screen2gene,\%gene2screen);
}



sub Gene2UnrelatedDis{
	my ($dis_ref,$gene_ref) = @_;
	my %dis2gene = %{$dis_ref};
	my %gene2dis = %{$gene_ref};

	my @dis_all = &Unique(keys %dis2gene);
	@dis_all = &Diff(\@dis_all,[0,999]);
	foreach my $gene (keys %gene2dis){
		$gene2unrelated{$gene} = [&Diff(\@dis_all,$gene2dis{$gene})];
	
	}
	return %gene2unrelated;

}


sub DisorderCount{
	my ($file,$pipeline) = @_;
	my %disorder_count = ();
	chomp $pipeline;

	open(IN,"<$file") or die "can not open $file\n";
	while(my $line = <IN>){
		chomp $line;
		my @aa = split /\t/,$line;
		#1,2,7
		if($aa[0] eq $pipeline){
			$disorder_count{$aa[1]} = $aa[5];
		}
		

	}
	
	return %disorder_count;

}

sub ValueSum{
	my ($list_ref,$hash_ref) = @_;
	my @list = @{$list_ref};
	my %hash = %{$hash_ref};	
	my $sum = 0;
	
	foreach my $i  (@list){
		if(defined $hash{$i}){
			$sum = $sum + $hash{$i};
		}
	}	

	return $sum;
}


sub VarGeneCount{

	my $varfile = shift @_;
	
	my $grep1 = "grep -e FN+FP -e TP+FP $varfile |grep incorrect| grep not_flagged -v | cut -f 1,14 | uniq";
        my @var = `$grep1`;
	my $grep2 = "grep CNV $varfile |grep -e FN+FP -e TP+FP | grep incorrect | cut -f 1,14 | uniq";
	my @var2 = `$grep2`;
	push (@var,@var2);
	
	chomp @var;
	my %genecount = ();
	foreach my $l (@var){
		my @aa = split("\t",$l);
		if($aa[1]){
			$genecount{$aa[1]}++;
		}
	}
	return %genecount;
}


