#!usr/bin/perl -w
use strict;
my ($file,$prevalence,$dis_group_file,$sensitivity_file,$spec_file) = @ARGV;
open(IN,"<$file") or die "can not open $file\n";
my %DIS = ();
my %SC = ();
my %state = ();
my %PREV = ();
my %FP = ();
my %FP_screen = ();
my %FP_state = ();
my $unaffected = 0;


open(P,"$prevalence") or die "can not open $prevalence\n";
while(my $line = <P>){
        chomp $line;
        my @aa = split /\t/,$line;
        my $dis = $aa[0];
        $PREV{$dis}=$aa[1];
        
}       
close P;


while(my $line = <IN>){
	#skip header line
	if($line =~ /^sample/){
		next;
	}
	
	chomp $line;
	my @aa = split /\t/,$line;
	my $sample =  $aa[0];
	my $dis = $aa[4];
	my $predic_dis = $aa[7];
	if($aa[4] =~ /;/){
		my @bb = split /;/,$aa[4];
		$dis = $bb[0];
		$dis =~ s/\s+//g;
	}

	if(($dis eq "77") and ($aa[1] eq "FN") and ($aa[9] eq "Y")){
		$dis = 0;
	}
	
	#unaffected cases:
	if($dis eq "0"){
		$unaffected++;
		if($aa[1] eq "FP"){
			$FP{$predic_dis}++;
		}
		if($aa[2] eq "FP"){
			$FP_screen{$predic_dis}++;
		}
		if($aa[3] eq "FP"){
			$FP_state{$predic_dis}++;
		}
		next;
	}
	
	#affected cases:
	elsif($dis ne "999"){
		$DIS{$dis}{"affect"}++;
		$DIS{$dis}{"affect_w"} += $PREV{$dis};
		$SC{$dis}{"affect"}++;
		$state{$dis}{"affect"}++;
	
		if(($aa[1] eq "TP") or ($aa[1] eq "TP+FN")){
			$DIS{$dis}{"TP"} ++;
			$DIS{$dis}{"TP_w"} += $PREV{$dis};
		}

		if(($aa[2] eq "TP") or ($aa[2] eq "TP+FN")){
			$SC{$dis}{"TP"}++;
                	$SC{$dis}{"TP_w"} += $PREV{$dis};
		}
		if(($aa[3] eq "TP") or ($aa[3] eq "TP+FN")){
			$state{$dis}{"TP"}++;
			$state{$dis}{"TP_w"} += $PREV{$dis};
		}
	}
}
close IN;

open(SEN_FILE,">$sensitivity_file") or die "can not open $sensitivity_file\n";
print SEN_FILE "Sample set\tPrevalence weighted sensitivity(d)\tSensitivity(d)\tPrevalence weighted sensitivity(s)\tSensitivity(s)\tPrevalence weighted sensitivity(a)\tSensitivity(a)\tTP(d)\tTP(s)\tTP(a)\tTotal affected individuals\n";

my ($dis_group_ref,$group_list_ref) = &ReadGroup($dis_group_file);
	my %dis_group = %{$dis_group_ref};
	my @group_list = @{$group_list_ref};

	foreach my $group (@group_list){
		my $affected_sum = 0;
		my $affected_sum_w = 0;
		my $TP_sum = 0;
		my $TP_sum_w = 0;
		my $TP_screen_sum = 0;
		my $TP_screen_sum_w = 0;
		my $TP_state_sum = 0;
		my $TP_state_sum_w = 0;
		foreach my $dis (@{$dis_group{$group}}){
			if(!$DIS{$dis}{'affect'}){next;}
			$affected_sum += $DIS{$dis}{'affect'};
			$affected_sum_w += $DIS{$dis}{'affect_w'};
			if(!$DIS{$dis}{'TP'}){
				$DIS{$dis}{'TP'} = 0;
				$DIS{$dis}{'TP_w'} = 0;
			}
			if(!$SC{$dis}{'TP'}){
				$SC{$dis}{'TP'} = 0;
				$SC{$dis}{'TP_w'} = 0;

			}

			if(!$state{$dis}{'TP'}){
                                $state{$dis}{'TP'} = 0;
                                $state{$dis}{'TP_w'} = 0;

                        }
			$TP_sum += $DIS{$dis}{'TP'};
			$TP_sum_w += $DIS{$dis}{'TP_w'};
			$TP_screen_sum += $SC{$dis}{"TP"};
			$TP_screen_sum_w += $SC{$dis}{"TP_w"};
			$TP_state_sum += $state{$dis}{"TP"};
			$TP_state_sum_w += $state{$dis}{"TP_w"};
			}
				
	
	if($affected_sum != 0){
		
	my $SEN = &Percentage($TP_sum/$affected_sum);
	my $SENw = &Percentage($TP_sum_w/$affected_sum_w);	
	my $SEN_screen = &Percentage($TP_screen_sum/$affected_sum);
	my $SENw_screen = &Percentage($TP_screen_sum_w/$affected_sum_w);
	my $SEN_state = &Percentage($TP_state_sum/$affected_sum);
	my $SENw_state = &Percentage($TP_state_sum_w/$affected_sum_w);
	my $dis_list = join(",",@{$dis_group{$group}}); 

	
	print SEN_FILE "$group\t$SENw\t$SEN\t$SENw_screen\t$SEN_screen\t$SENw_state\t$SEN_state\t$TP_sum\t$TP_screen_sum\t$TP_state_sum\t$affected_sum\n";
	}
	else{
		my($SEN,$SENw,$SEN_screen,$SENw_screen,$SEN_state,$SENw_state) = ("NA","NA","NA","NA","NA","NA");
		my $dis_list = join(",",@{$dis_group{$group}}); 
		print SEN_FILE "$group\t$SENw\t$SEN\t$SENw_screen\t$SEN_screen\t$SENw_state\t$SEN_state\t$TP_sum\t$TP_screen_sum\t$TP_state_sum\t$affected_sum\n"; 
	}

}

##################################
#FP on unaffected cases
open(SPEC,">$spec_file") or die "can not open $spec_file\n";
print SPEC "Sample set\tFraction(d)\tFraction(s)\tFraction(a)\tExome flagged individuals(d)\tExome flagged individuals(s)\tExome flagged individuals(a)\tTotal unaffected individuals\n";
foreach my $group (sort keys %dis_group){
	my $FP_sum = 0;
	my $FP_screen_sum = 0;
	my $FP_state_sum = 0;
	foreach my $i (keys %FP){
		my @cc = split ",",$i;
			my @dd = &Intersect(\@cc,\@{$dis_group{$group}});
			if($dd[0]){
				$FP_sum += $FP{$i};
				$FP_screen_sum += $FP_screen{$i};
				$FP_state_sum += $FP_state{$i};		
			}
		}

	if($unaffected){
		my $SPEC = &Percentage($FP_sum/$unaffected);
		my $SPEC_screen = &Percentage($FP_screen_sum/$unaffected);
		my $SPEC_state = &Percentage($FP_state_sum/$unaffected);
		my $dis_list = join(",",@{$dis_group{$group}});
		print SPEC "$group\t$SPEC\t$SPEC_screen\t$SPEC_state\t$FP_sum\t$FP_screen_sum\t$FP_state_sum\t$unaffected\n";
	}

}



##################################
#raw SCAD key

my $affected_sum = 0;
my $affected_sum_w = 0;
my $TP_sum = 0;
my $TP_sum_w = 0;
my $TP_screen_sum = 0;
my $TP_screen_sum_w = 0;
my $TP_state_sum = 0;
my $TP_state_sum_w = 0;
$unaffected = 0;
my $FP_sum = 0;
my $FP_screen_sum = 0;
my $FP_state_sum = 0;
open(IN,"<$file") or die "can not open $file\n";
while(my $line = <IN>){
	#skip header line
	if($line =~ /^sample/){
		next;
	}
	
	chomp $line;
	my @aa = split /\t/,$line;
	my $sample =  $aa[0];
	my $dis = $aa[4];
	my $predic_dis = $aa[7];
	if($aa[4] =~ /;/){
		my @bb = split /;/,$aa[4];
		$dis = $bb[0];
		$dis =~ s/\s+//g;
	}

	
	#unaffected cases:
	if($dis eq "0"){
		$unaffected++;
		if($aa[1] eq "FP"){
			$FP_sum ++;
		}
		if($aa[2] eq "FP"){
			$FP_screen_sum ++;
		}
		if($aa[3] eq "FP"){
			$FP_state_sum ++;
		}
		next;
	}
	
	#affected cases:
	elsif($dis ne "999"){
		$affected_sum ++;
		$affected_sum_w += $PREV{$dis};
	
		if(($aa[1] eq "TP") or ($aa[1] eq "TP+FN")){
			$TP_sum ++;
			$TP_sum_w += $PREV{$dis};
		}

		if(($aa[2] eq "TP") or ($aa[2] eq "TP+FN")){
			$TP_screen_sum ++;
                	$TP_screen_sum_w += $PREV{$dis};
		}
		if(($aa[3] eq "TP") or ($aa[3] eq "TP+FN")){
			$TP_state_sum++;
			$TP_state_sum_w += $PREV{$dis};
		}
	}
}
close IN;


	if($affected_sum != 0){
		
	my $SEN = &Percentage($TP_sum/$affected_sum);
	my $SENw = &Percentage($TP_sum_w/$affected_sum_w);	
	my $SEN_screen = &Percentage($TP_screen_sum/$affected_sum);
	my $SENw_screen = &Percentage($TP_screen_sum_w/$affected_sum_w);
	my $SEN_state = &Percentage($TP_state_sum/$affected_sum);
	my $SENw_state = &Percentage($TP_state_sum_w/$affected_sum_w);
	my $dis_list = join(",",@{$dis_group{"Full set"}} );

	
	print SEN_FILE "Full set(raw SCAD)\t$SENw\t$SEN\t$SENw_screen\t$SEN_screen\t$SENw_state\t$SEN_state\t$TP_sum\t$TP_screen_sum\t$TP_state_sum\t$affected_sum\n";
	}
	else{
		my($SEN,$SENw,$SEN_screen,$SENw_screen,$SEN_state,$SENw_state) = ("NA","NA","NA","NA","NA","NA");
		my $dis_list = join(",",@{$dis_group{"Full set"}}); 
		print SEN_FILE "Full set(raw SCAD)\t$SENw\t$SEN\t$SENw_screen\t$SEN_screen\t$SENw_state\t$SEN_state\t$TP_sum\t$TP_screen_sum\t$TP_state_sum\t$affected_sum\n"; 
	}

close SEN_FILE;


 if($unaffected){
		my $SPEC = &Percentage($FP_sum/$unaffected);
		my $SPEC_screen = &Percentage($FP_screen_sum/$unaffected);
		my $SPEC_state = &Percentage($FP_state_sum/$unaffected);
	
		my $dis_list = join(",",@{$dis_group{"Full set"}});
        print SPEC "Full set(raw SCAD)\t$SPEC\t$SPEC_screen\t$SPEC_state\t$FP_sum\t$FP_screen_sum\t$FP_state_sum\t$unaffected\t$dis_list\n";
	}

close SPEC;



##################################
#read in disorder groups
sub ReadGroup{
	my $group = shift @_;
	chomp $group;
	my %dis_group = ();
	my @group_list = ();
	open(GROUP,"<$group") or die "can not open $group\n";
	while(my $line = <GROUP>){
		chomp $line;
		my @aa = split /\t/,$line;
		my @bb = split /,/,$aa[1];
		$dis_group{$aa[0]} = \@bb;
		push @group_list,$aa[0];
	}
close GROUP;
return (\%dis_group,\@group_list);
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
sub Percentage{
	my $i = shift @_;
	my $percent = sprintf("%.2f%%",$i*100);
	return $percent;
}
