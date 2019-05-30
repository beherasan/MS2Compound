#!/usr/bin/perl

## Usage ##
## C:\Users\santosh\Desktop\Test_delete\Test_delete\bin\Debug\perl.exe MS2Compound.pl D:\CSBMM\PhD\MSMS_compoundIdentifier\negative\NR_B_subtilis_neg_031119 D:\CSBMM\PhD\MSMS_compoundIdentifier\negative\pre_relativeIntensity.mgf D:\CSBMM\PhD\MSMS_compoundIdentifier\negative\B_sub_monoIsotopicmass_051519.txt 0 0.05 0.5 D:\CSBMM\PhD\MSMS_compoundIdentifier\negative
##

$db=$ARGV[0];
$query_mgf=$ARGV[1];
$monoMass_file= $ARGV[2];
$energy=$ARGV[3];
$error=$ARGV[4];
$frag_error=$ARGV[5];
$out_dir = $ARGV[6];
%hash=();


opendir(DIR, $db) or die "Could not open $db\n";
while ($fp1 = readdir(DIR)) {
	$fp = $db."\\".$fp1; 
	if($fp =~/\.log$/)
	{
		open my $fh, "<", $fp or die "can't read open '$fp': $OS_ERROR";
		$DBcompound= (split /\\/,$fp)[-1];
		$DBcompound=~s/\.log//;
		while (<$fh>)
		{
			chomp;
			unless(/^\s/)
			{
				if(/energy0/)
				{
					$switch=0;
				}
				elsif(/energy1/)
				{
					$switch=1;	
				}
				elsif(/energy2/)
				{
					$switch=2;	
				}
				else
				{
					($frag,$inten)=(split /\s/)[0,1];
					$hash{$switch}{$DBcompound}.="$frag-$inten;";
				}	
			}
		}
		close $fh or die "can't read close '$fp': $OS_ERROR";		
	}
}
closedir(DIR);

## Read MonoIsotopic parent mass ##
open(IN,"$monoMass_file") or die "Could not open the file:$!\n";
while(<IN>)
{
	chomp;
	($id,$monoMass)=(split /\t/)[0,1];
	$hashMono{$id}="$monoMass";
}
close IN;

##### MS/MS match #####
open(MGF,"$query_mgf") or die "Could not open the file:$!\n";
open(OUT,">$out_dir\\output.txt") or die "Could not create the file:$!\n";

print OUT "FEATURE_ID\tCompoundId:Mean_Deviation_Score:NumberOfMatch||\n";
while(<MGF>)
{
	chomp;
	if(/^BEGIN IONS/)
	{
		## inititae all the variables here for each instance of metabolite
		@parentArr=();
		$frag_mz_int_exp = "";
	}
	if(/^FEATURE_ID\=(.*)/)
	{
		$id=$1;
	}
	if(/^PEPMASS\=(.*)/)
	{
		$precursor_mass=$1;
		$pre_up = $precursor_mass + $error;
		$pre_down = $precursor_mass - $error;
		foreach	$m (keys%hashMono)
		{
			$mono = $hashMono{$m};
			if($mono >= $pre_down && $mono <= $pre_up)
			{
				push(@parentArr,$m);
			}	
		}
	}
	if(/^CHARGE\=(.*)/)
	{
		$charge_all=$1;
		if($charge_all =~/\+/)
		{
			$charge_num=$charge_all;$charge_num=~s/\+//g;
			$charge="+";
		}
		elsif($charge_all =~/\-/)
		{
			$charge_num=$charge_all;$charge_num=~s/\-//g;
			$charge="-";
		}
		else
		{
			print "Check the charge state in MGF file\n";
			exit;	
		}	
	}
	if(/^[0-9]/)
	{
		($frag_mz,$frag_int)=(split /\s/)[0,1];
		$frag_mz_int_exp .= "$frag_mz-$frag_int;";
	}
	if(/^END IONS/)
	{
		$sw_match=0;
		if (scalar(@parentArr) >= 1)
		{
			@frag_exp_arr = (split /\;/,$frag_mz_int_exp);
			foreach (@parentArr)
			{
				$raw_score=$match_count=0;
				$parent_id = $_;
				$frag_theo = $hash{$energy}{$_};
				@frag_theo_arr = (split /\;/,$frag_theo);
				foreach $exp (@frag_exp_arr)
				{
					($exp_mz,$exp_int)=(split /\-/,$exp)[0,1];
					
					$exp_mz_down = $exp_mz - $frag_error;
					$exp_mz_up = $exp_mz + $frag_error;
					foreach $theo (@frag_theo_arr)
					{
						($theo_mz,$theo_int) =(split /\-/,$theo)[0,1];
						if($theo_mz >= $exp_mz_down && $theo_mz <= $exp_mz_up)
						{
							$score=(abs($theo_mz - $exp_mz)) *(abs ($theo_int - $exp_int));
							$raw_score+=$score;
							$match_count++;
							$sw_match=1;
						}
					}			
				}
				unless ($match_count == 0)
				{
					$mean_dev = $raw_score/$match_count;
					$hash_rawScore{$id}.="$parent_id:$mean_dev:$match_count||";
				}
				else
				{
					$mean_dev = 0;
				}	
			}	
		}
		if($sw_match == 1)
		{
			print OUT "$id\t$hash_rawScore{$id}\n";
		}
	}
}
close OUT;
close MGF;
#############################
### End of the program ######
#############################
