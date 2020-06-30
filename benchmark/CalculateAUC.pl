#!/usr/bin/perl

use strict;
use warnings; 

my $infile = "TPR_FPR_Neg.txt";
my %hash=();

open(IN,"$infile") or die "Could not open the file:$!\n";
while(<IN>)
{
	chomp;
	my ($tpr,$fpr) = (0,0);
	my $type="";
	unless(/^Type	TPR	FPR/)
	{
		($type, $tpr, $fpr)=(split /\t/)[0,1,2];
		unless(exists $hash{$type}{$fpr})
		{
			$hash{$type}{$fpr}="$tpr";			
		}
		else
		{
			my $v = $hash{$type}{$fpr};
			if($tpr > $v)
			{
				$hash{$type}{$fpr}="$tpr";
			}	
		}
	}
}
close IN;

foreach my $k1 (sort(keys %hash))
{
	print "$k1";
	my @arr= sort(keys %{$hash{$k1}});
	my $auc=0;
	for(my $i=1; $i<(scalar(@arr));$i++)
	{
		#print "$arr[$i]\t$hash{$k1}{$arr[$i]}\n";
		my $s = (($arr[$i])-($arr[$i-1]))*((($hash{$k1}{$arr[$i]}) + ($hash{$k1}{$arr[$i-1]}))/2);
		$auc+=$s;
	}
	print "\t$auc\n";
	print "********\n";
}