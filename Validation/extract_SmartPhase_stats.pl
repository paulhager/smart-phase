#!/bin/perl -w

#  Copyright (C) 2018 the SmartPhase contributors.
#  Website: https://github.com/paulhager/smart-phase
#
#  This file is part of the SmartPhase phasing tool.
#
#  The SmartPhase phasing tool is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

use strict;

# call: perl extract_SmartPhase_stats.pl dir <CEU|YRI> <NA..> <chr..> <trio|NOtrio> 10 0.1

my $dir = $ARGV[0];
my $fam = $ARGV[1];
my $sample = $ARGV[2];
my $chr = $ARGV[3];
my $trio = $ARGV[4];
my $maxVarPerGene = $ARGV[5];
my $minConf = $ARGV[6];

my $genes = 0;
my $variants = 0;
my $total_variants = 0;
my $potpairs = 0;
my $pc = 0; #phase connections
my $pcr = 0;
my $flag;
my $score;
my $tp;
my $innocuous = 0;
my $cis = 0;
my $trans = 0;
my $notphased = 0;
my $notfound = 0;
my @fields = ();
my %excludeIntervals = ();
my %pairToFlag = ();

my $invalcsv = "$dir/sim.$fam.trio.chr$chr.phased.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.SmartPhase_${trio}_VALIDATION.csv";
my $inresults = "$dir/sim.$fam.trio.chr$chr.phased.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.SmartPhase_$trio.tsv";
my $invalconftsv = "$dir/sim.$fam.trio.chr$chr.phased.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.SmartPhase_${trio}_VALIDATION_CONF.tsv";
my $out = $inresults;
$out =~ s/.tsv/.evaluation.txt/ig;

my $eval_pairs = "$dir/sim.$fam.trio.chr$chr.phased.intersect.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37.recode.evaluation_pairs.tsv";
#my $spout_good_intervals = $inresults;
#$spout_good_intervals =~ s/.tsv/_good_intervals.tsv/ig;

# open output file
open(my $outfh, '>', $out) or die "Could not open outfile ($out)!";


# read VALIDATION.csv
open(my $fh, '<', $invalcsv) or die "Could not open file ($invalcsv)!";

while(my $line = <$fh>) {
	chomp $line;
	@fields = split(',', $line);
	$total_variants+=$fields[1];
	if($fields[1] > $maxVarPerGene) {
		print $outfh "WARNING: Interval $fields[0] will be excluded because it contains $fields[1] variants!\n";
		$excludeIntervals{$fields[0]} = 1;
	} else {
		$variants+=$fields[1];
		$genes++;
		$potpairs+=(($fields[1]*($fields[1]-1))/2);
		$pc+=$fields[2];
	}
}

close($fh);

$pcr = $pc/($variants-$genes);

print $outfh "In total, there are $total_variants variants in all regions.\n";
print $outfh "There are $genes genes, a total of $variants variants resulting in $potpairs potential compound heterozygous variant pairs.\n";
print $outfh "The phase connnection ratio is: $pcr\n";

# write SP results of evaluation intervals
open(my $fhout2, '>', $eval_pairs) or die "Could not open file ($eval_pairs)!";

# read results.txt
open($fh, '<', $inresults) or die "Could not open file ($inresults)!";

while(my $line = <$fh>) {
	chomp $line;
	@fields = split('\t', $line);
	if(scalar(@fields) == 5 && ! $excludeIntervals{$fields[0]}) {
		print $fhout2 "$sample\t$fields[1]\t$fields[2]\n";
		#print $fhout2 "$line\n";
		$flag = $fields[3];
		$pairToFlag{"$fields[0]+$fields[1]+$fields[2]"} = $flag;
		# determine innocuous pairs
		if($flag == 9 || $flag == 10 || $flag == 12 || $flag == 25 || $flag == 26 || $flag == 28) {
			$innocuous++;
		} elsif($flag == 1 || $flag == 17 ) {
			$cis++;
		} elsif($flag == 2 || $flag == 18 ) {
			$trans++;
		} elsif($flag == 4) {
			$notphased++;
		} elsif($flag == 20)  {
			$notfound++;
		}
		
		
	}
}

close($fhout2);

print $outfh "Pairs by flag:\n";
print $outfh "Innocuous\t$innocuous\n";
print $outfh "Cis\t$cis\n";
print $outfh "Trans\t$trans\n";
print $outfh "Not phased\t$notphased\n";
print $outfh "Not found\t$notfound\n";

# read VALIDATION_CONF.tsv

my $cis_fp = 0;
my $trans_fp = 0;
my $cis_lowConf = 0;
my $trans_lowConf = 0;
my $cis_fp_lowConf = 0;
my $trans_fp_lowConf = 0;


open($fh, '<', $invalconftsv) or die "Could not open file ($invalconftsv)!";

while(my $line = <$fh>) {
        chomp $line;
        @fields = split('\t', $line);
        if(! $excludeIntervals{$fields[0]}) {
		$flag = $pairToFlag{"$fields[0]+$fields[1]+$fields[2]"};
		$score = $fields[3];
		$tp = $fields[4];
		if($flag == 1 || $flag == 17) {
			if(!$tp) {
				$cis_fp++;
			}
			if($score < $minConf) {
				$cis_lowConf++;
				if(!$tp) {
					$cis_fp_lowConf++;
				}
			}
		} elsif($flag ==2 || $flag == 18) {
			if(!$tp) {
				$trans_fp++;
			}
			if($score < $minConf) {
				$trans_lowConf++;
				if(!$tp) {
					$trans_fp_lowConf++;
				}
			}
		}
	}
}

close($fh);

print $outfh "For $cis cis predictions, there are $cis_fp false-positives and $cis_lowConf with score below $minConf (overlap $cis_fp_lowConf).\n";
print $outfh "For $trans trans predictions, there are $trans_fp false-positives and $trans_lowConf with score below $minConf (overlap $trans_fp_lowConf).\n";

close($outfh);
