#!/usr/bin/perl
use strict;
use warnings;

my $infile=$ARGV[0];
my $outfile=$ARGV[1];
open(IN,"$infile") or die "could not open the file:$!\n";
open(OUT,">$outfile") or die "Could not create the file:$!\n";

my $max_intensity;
my @mz_arr;my @int_arr;
while(<IN>)
{
	chomp;
	$_=~s/\r//g;
	if(/^[^0-9]/ && /^[^(END IONS)]/)
	{
		print OUT "$_\n";
		$max_intensity=0;
		@mz_arr=@int_arr=();
	}
	elsif(/^[0-9]/)
	{
		my ($mz,$int)=(split /\s/)[0,1];
		if($int >= $max_intensity)
		{
			$max_intensity = $int;
		}
		push(@mz_arr,$mz);
		push(@int_arr,$int);
	}
	elsif(/^END IONS/)
	{
		for (my $i=0;$i < scalar(@mz_arr);$i++)
		{
			print OUT "$mz_arr[$i]";
			my $relativeInt=($int_arr[$i]/$max_intensity)*100;
			print OUT " $relativeInt\n";
		}
		print OUT "END IONS\n\n";
	}
}

close IN;
close OUT;
