#!/usr/bin/perl

use strict;
use warnings;

my %hashMono=();
open(MONO,"MonoIsotopicMass.txt") or die "Could not open the file:$!\n";
while(<MONO>)
{
	chomp;
	$_=~s/\r//g;
	unless(/PubChemId/)
	{
		my ($cha,$pubchem)=(split /\t/)[0,1];
		unless($pubchem eq "NA")
		{
			$hashMono{$cha}="$pubchem";
		}
		else
		{
			$hashMono{$cha}="$cha";			
		}
	}
}
close MONO;

open(IN,"PositiveScore0.05.txt") or die "Could not open the file:$!\n";
open(OUT,">TPR_FPR_Pos.txt") or die "Could not create the file:$!\n";
my $ref="";
my %hashSEN=();
my %hashSPE=();
my %hashRank=();
my (%rank, %rankShared, %Shared)=()x3;
my $highRank=0;
my ($NumTP,$NumTN)=(0,0);
my @scores=("DP","mS","HGS","FS");
while(<IN>)
{
	chomp;
	$_=~s/\r//;
	if(/^\#/)
	{
		$ref=$_;
		$ref=~s /\#//;
	}
	elsif(/^[A-Z|a-z|0-9]/)
	{
		my @line=(split /\t/);
		my $candidate=shift(@line);
		for(my $i=0;$i<scalar(@line);$i++)
		{
			if(($ref ne $candidate) && ($hashMono{$ref} eq $hashMono{$candidate}))
			{				
			}
			else
			{
				unless((exists $hashRank{$scores[$i]}) && (exists $hashRank{$scores[$i]}{$line[$i]}))
				{
					$hashRank{$scores[$i]}{$line[$i]} = "$candidate";
				}
				else
				{
					$hashRank{$scores[$i]}{$line[$i]} .= "|$candidate";
				}
			}	
		}
	}
	else
	{	
		foreach my $k1 (keys %hashRank)
		{
			my @k2All = sort{$b <=> $a}(keys %{$hashRank{$k1}});
			
			#print OUT "Santosh\t$k1\t@k2All";
			#print OUT "BEHERA @k2All\n";
			for(my $i=0;$i<scalar(@k2All);$i++)
			{
				my ($TP,$TN)=(0,0);
				my $NumCan=0;
				#print OUT "SANTOSH$i\n";
				if($hashRank{$k1}{$k2All[$i]}=~/\|/)
				{
					#$NumCan=scalar(split /\|/,$hashRank{$k1}{$k2All[$i]});
					#print OUT "##########\n";
					#print OUT "SANTOSH\n$hashRank{$k1}{$k2All[$i]}\n";
					$NumCan = $hashRank{$k1}{$k2All[$i]} =~ tr/\|//;
					#print OUT "$NumCan\n";
					if($hashRank{$k1}{$k2All[$i]}=~/$ref/)
					{
						$TP=1;
						$NumTP++;
						$TN+=($NumCan);
						#print OUT "$NumTP\t$NumTN\n"
					}
					else
					{
						$TN+=($NumCan+1);
					}
				}
				else
				{
					#$NumCan=1;
					if($hashRank{$k1}{$k2All[$i]}=~/$ref/)
					{
						$TP++;
						#$NumTN+=($NumCan-1);
					}
					else
					{
						$TN++;
					}					
				}
				my $j=$i+1;
				#print OUT "$j\n";
				#print OUT "$NumTP\t$NumTN\n";
				$highRank = $j if ($j >= $highRank);
				$hashSEN{$k1}{$j}+=$TP;
				$hashSPE{$k1}{$j}+=$TN;
				#print OUT "$hashRank{$k1}{$k2All[$i]}\n";
				#if($hashRank{$k1}{$k2All[$i]}=~/$ref/)
				#{
				#	my $j=$i+1;
				#	print OUT "B\t$j\n";
				#	$rank{$k1}{$j}++;
				#	$highRank = $j if ($j >= $highRank);
				#	if($hashRank{$k1}{$k2All[$i]}=~/\|/)
				#	{
				#		my $count=0;
				#		$count = $hashRank{$k1}{$k2All[$i]} =~ tr/\|//;
				#		$rankShared{$k1}{$j}+=($count);
				#		$Shared{$k1}{$j}++;# For NR shared
				#	}
				#}
			}	
		}
		%hashRank=();
	}
}
close IN;

#print "$highRank\n";

#my ($TotalTP,$TotalTN)=(0,0);
my (%hashTN, %hashTP)=()x2;
foreach my $k1 (@scores)
{
	my ($TotTP,$TotTN)=(0,0);
	#print OUT "#$k1\n";
	foreach my $k2 (1..$highRank)
	{
		$TotTP+=$hashSEN{$k1}{$k2};
		$TotTN+=$hashSPE{$k1}{$k2};
		#print OUT "$k2\t$hashSEN{$k1}{$k2}\t$hashSPE{$k1}{$k2}\n";
	}
	$hashTP{$k1}="$TotTP";
	$hashTN{$k1}="$TotTN";
}

print OUT "Type\tTPR\tFPR\n";
foreach my $k1 (@scores)
{
	#print OUT "#$k1\t$hashTP{$k1}\t$hashTN{$k1}\n";
	my ($tpr_backup, $fpr_backup)=(0,0);
	print OUT "$k1\t$tpr_backup\t$fpr_backup\n";	
	foreach my $k2 (1..$highRank)
	{
		my $tpr = ((($hashSEN{$k1}{$k2})/$hashTP{$k1}) + $tpr_backup);
		my $fpr = ((($hashSPE{$k1}{$k2})/$hashTN{$k1}) + $fpr_backup);
		$tpr_backup = $tpr;
		$fpr_backup = $fpr;
		#my ($tpr,$fpr)=(0,0);
		print OUT "$k1\t$tpr\t$fpr\n";
	}
}

close OUT;


exit;
