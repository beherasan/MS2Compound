#!/usr/bin/perl
#
## Usage: perl ScoreComparison_Pos.pl casmi2017_pos
###########

use strict;
use warnings;

my $Refdir="casmi2017_pos";
my $preTol=0.05;
my $fraTol=0.05;
my %hashMatchedData=();

my %hashPreMass=(); ## Precursor mass of Challenge data
open(IN,"PrecursorMass.txt") or die "Could not open the file:$!\n";
while(<IN>)
{
	chomp;
	$_=~s/\r//g;
	unless(/^challengename/)
	{
		my ($name,$pre)=(split /\t/)[0,1];
		$hashPreMass{$name}="$pre";
	}
}
close IN;

my %hashMolWeight=(); ## Molecular weight of the PubChem IDs to be matched.
open(MW,"MonoIsotopicMass.txt") or die "Could not open the file:$!\n";
while(<MW>)
{
	chomp;
	$_=~s/\r//g;
	unless(/PubChemId/)
	{
		my ($id,$mw) = (split /\t/)[0,2];
		$hashMolWeight{$id}="$mw";
	}	
}
close MW;

my %hashRef=(); ## Store all the predicted spectra to be match
my %hashRaw=(); ## Store the raw file, which is used during Fit score calculation
foreach my $fp (glob("$Refdir/*intensityRatio.log"))
{
	open(REF,"$fp") or die "Could not open the file:$!\n";
	my ($energy,$mz,$int)=("")x3;
	my $f=$fp;
	$f=~s/$Refdir\///;$f=~s/\_intensityRatio\.log//;
	while (<REF>)
	{
		chomp;
		$_=~s/\r//g;
		if(/^energy(.*)/)
		{
			$energy=$1;
		}
		else
		{
			($mz, $int)=(split /\s/)[0,1];
			if(defined $mz and defined $int)
			{
				if((exists $hashRef{$f}) and (exists $hashRef{$f}{$energy}))
				{
					$hashRef{$f}{$energy}.=",$mz-$int";
				}
				else
				{
					$hashRef{$f}{$energy}="$mz-$int";
				}
			}
		}
	}
	close REF;
}

my @parent=();
my (%hash_N,%hash_N1,%hash_K,%hash_K1)=()x4;
foreach my $cha (keys %hashPreMass)
{
	my $mz=$hashPreMass{$cha};
	my $mz_up=$mz+$preTol;
	my $mz_down=$mz-$preTol;
	foreach my $mw (keys %hashMolWeight)
	{		
		if($hashMolWeight{$mw} >= $mz_down && $hashMolWeight{$mw} <= $mz_up)
		{
			push(@parent,$mw);
		}
		else
		{
			my $mz_new = $mz-1.007276; # M+H
			$mz_up=$mz_new+$preTol;
			$mz_down=$mz_new-$preTol;
			if($hashMolWeight{$mw} >= $mz_down && $hashMolWeight{$mw} <= $mz_up)
			{
				push(@parent,$mw);
			}
			$mz_new = $mz-22.989218; # M+Na
			$mz_up=$mz_new+$preTol;
			$mz_down=$mz_new-$preTol;
			if($hashMolWeight{$mw} >= $mz_down && $hashMolWeight{$mw} <= $mz_up)
			{
				push(@parent,$mw);
			}
			$mz_new = $mz-18.033823; # M+NH4
			$mz_up=$mz_new+$preTol;
			$mz_down=$mz_new-$preTol;
			if($hashMolWeight{$mw} >= $mz_down && $hashMolWeight{$mw} <= $mz_up)
			{
				push(@parent,$mw);
			}
			$mz_new = $mz-38.963158; # M+K
			$mz_up=$mz_new+$preTol;
			$mz_down=$mz_new-$preTol;
			if($hashMolWeight{$mw} >= $mz_down && $hashMolWeight{$mw} <= $mz_up)
			{
				push(@parent,$mw);
			}
			$mz_new = $mz-42.033823; #M+ACN+H
			$mz_up=$mz_new+$preTol;
			$mz_down=$mz_new-$preTol;
			if($hashMolWeight{$mw} >= $mz_down && $hashMolWeight{$mw} <= $mz_up)
			{
				push(@parent,$mw);
			}
			$mz_new = $mz-44.971160; #M+2Na-H
			$mz_up=$mz_new+$preTol;
			$mz_down=$mz_new-$preTol;
			if($hashMolWeight{$mw} >= $mz_down && $hashMolWeight{$mw} <= $mz_up)
			{
				push(@parent,$mw);
			}			
		}
	}
	my $rawFile="Positive/".$cha."-msms.mgf";
	if(scalar(@parent) >= 1 && -e $rawFile)
	{
		my $rawData="";
		open(RAW,"$rawFile") or die "Could not open the file:$!\n";
		while(<RAW>)
		{
			chomp;
			$_=~s/\r//g;
			if(/^[0-9]/)
			{
				my ($mz, $inten)=(split /\t/)[0,1];
				if(length($rawData) >= 1)
				{
					$rawData=$rawData.",".$mz."-".$inten;
				}
				else
				{
					$rawData=$mz."-".$inten;
				}		
			}
		}
		close RAW;
		my ($count_N,$count_N1,$count_K,$count_K1)=(0,0,0,0);
		foreach my $p (@parent)
		{
			my $pData=$hashRef{$p}{"1"};
			$hashRaw{$cha}="$rawData";
			($count_N1,$count_K1)=(0,0);
			my @Raw=(split /\,/,$rawData);
			my @Par=(split /\,/,$pData);
			$count_N1=scalar(@Par);
			$hash_N1{$cha}{$p}=$count_N1;
			$count_N+=(scalar(@Par));
			foreach my $i (@Raw)
			{
				my ($i_mz,$i_int)=(split /\-/,$i)[0,1];
				foreach my $j (@Par)
				{
					my ($j_mz,$j_int)=(split /\-/,$j)[0,1];
					if((abs($i_mz-$j_mz))<=$fraTol)
					{
						if(exists $hashMatchedData{$cha}{$p})
						{
							$hashMatchedData{$cha}{$p}.="|$i_mz-$i_int-$j_mz-$j_int";
						}
						else
						{
							$hashMatchedData{$cha}{$p}="$i_mz-$i_int-$j_mz-$j_int";							
						}
					}
				}	
			}
			if(exists $hashMatchedData{$cha} && exists $hashMatchedData{$cha}{$p})
			{
				my $rD = removeDuplicates($hashMatchedData{$cha}{$p});
				my $tmp_K=scalar(split /\|/,$rD);
				$hashMatchedData{$cha}{$p}=$rD;
				$count_K+=$tmp_K;
				$count_K1+=$tmp_K;
				$hash_K1{$cha}{$p}=$count_K1;
			}	
		}
		$hash_N{$cha}=$count_N;
		$hash_K{$cha}=$count_K;
	}	
	@parent=();
}

open(OUT,">PositiveScore0.05.txt") or die "Could not create the file:$!\n";
foreach my $k1 (sort (keys %hashMatchedData))
{
	print OUT "#$k1\n";
	my $N=$hash_N{$k1};
	my $K=$hash_K{$k1};
	foreach my $k2 (keys %{$hashMatchedData{$k1}})
	{
		my $N1=$hash_N1{$k1}{$k2};
		my $K1=$hash_K1{$k1}{$k2};
		my $dot = DotProduct($hashMatchedData{$k1}{$k2});
		my $S = Sscore($hashMatchedData{$k1}{$k2});
		my $H = HyperGeometric($N,$K,$N1,$K1); 
		my $rData=$hashRaw{$k1};
		my $FS = FitScore($hashMatchedData{$k1}{$k2},$rData);
		print OUT "$k2\t$dot\t$S\t$H\t$FS\n";
	}
	print OUT "\n";
}
close OUT;

### Subroutines

sub removeDuplicates{
	my $all=$_[0];
	my @ToScore=(split /\|/,$all);
	my %hashRmDup=();
	foreach my $i (@ToScore)
	{
		my ($rw_mz,$rw_int,$rf_mz,$rf_int)=(split /\-/,$i)[0,1,2,3];
		if(exists $hashRmDup{$rf_mz})
		{
			my $v=$hashRmDup{$rf_mz};
			my ($rw_int_tmp)=(split /\-/,$v)[1];
			if($rw_int > $rw_int_tmp)
			{
				$hashRmDup{$rf_mz}="$rw_mz-$rw_int-$rf_mz-$rf_int";
			}
		}
		else
		{
			$hashRmDup{$rf_mz}="$rw_mz-$rw_int-$rf_mz-$rf_int";
		}
	}
	my $new="";
	foreach my $k (sort(keys %hashRmDup))
	{
		if($new eq "")
		{
			$new="$hashRmDup{$k}";
		}
		else
		{
			$new.="|$hashRmDup{$k}";
		}
	}
	return $new;
}

sub DotProduct{
	my $all=$_[0];
	my @ToScore=(split /\|/,$all);
	my $den=0;
	my ($rw_sq,$rf_sq)=(0,0);
	foreach my $i (@ToScore)
	{
		my ($rw_mz,$rw_int,$rf_mz,$rf_int)=(split /\-/,$i)[0,1,2,3];
		$den+=($rw_int*$rf_int);
		$rw_sq+=($rw_int*$rw_int);
		$rf_sq+=($rf_int*$rf_int);
	}
	my $score=($den)/(sqrt($rw_sq*$rf_sq));
	return $score;
}

sub Sscore{
	my $all=$_[0];
	my @ToScore=(split /\|/,$all);
	my $score=0;
	my ($rw_sq,$rf_sq)=(0,0);
	my $l=scalar(@ToScore);
	foreach my $i (@ToScore)
	{
		my ($rw_mz,$rw_int,$rf_mz,$rf_int)=(split /\-/,$i)[0,1,2,3];
		$score+=((abs($rw_mz - $rf_mz))*(100/($rw_int))*(abs ($rw_int - $rf_int)));
	}
	my $s=((sqrt($score))/($l*$l));
	if($l == 1)
	{
		$s=-(log($s+1));
	}
	else
	{
		$s=-(log($s));		
	}
	return $s;
}

sub HyperGeometric{
	my ($n,$k,$n1,$k1)=@_;
	my ($one)=(Factorial($k)/(Factorial($k1)*Factorial($k-$k1)));
	my ($two)=(Factorial($n-$k)/(Factorial($n1-$k1)*Factorial(($n-$k)-($n1-$k1))));
	my ($three)=(Factorial($n)/(Factorial($n1)*Factorial($n-$n1)));
	my $h=(($one*$two)/$three);
	unless($h == 0)
	{
		my $hscore=-(log($h));
		return $hscore;
	}
	else
	{
		my $hscore=10;
		return $hscore;
	}
}

sub FitScore{
	my ($MatchedFragments,$AllFragments)=@_;
	my @MF=(split /\|/,$MatchedFragments);
	my $one=0;
	foreach my $i (@MF)
	{
		my ($rw_mz,$rw_int)=(split /\-/,$i)[0,1];
		$one+=($rw_mz*$rw_int);		
	}
	my @AF=(split /\,/,$AllFragments);
	my $two=0;
	foreach my $i (@AF)
	{
		my ($rw_mz,$rw_int)=(split /\-/,$i)[0,1];
		$two+=($rw_mz*$rw_int);		
	}
	my $ftscore=($one/$two);	
	return $ftscore;	
}

sub Factorial{
	my $num=$_[0];
	my $fac = 1;
	for(my $i=1;$i<=$num;$i++)
	{
		$fac = $fac*$i;
	}
	return $fac;
}


exit;


