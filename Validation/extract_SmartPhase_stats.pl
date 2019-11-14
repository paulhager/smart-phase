#!/usr/bin/perl -w

use strict;

# call: perl extract_SmartPhase_stats.pl dir <CEU|YRI> <NA...> <chr..> <trio|NOtrio> 10 0.1

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

my $invalcsv = "$dir/sim.$fam.trio.$chr.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37_SmartPhase_${trio}_VALIDATION.csv";
my $inresults = "$dir/sim.$fam.trio.$chr.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37_SmartPhase_$trio.tsv";
my $invalconftsv = "$dir/sim.$fam.trio.$chr.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37_SmartPhase_${trio}_VALIDATION_CONF.tsv";
my $out = $inresults;
$out =~ s/.tsv/.evaluation.txt/ig;

my $eval_pairs = "$dir/sim.$fam.trio.$chr.phased.AGV6UTR.allGeneRegionsCanonical.HG19.GRCh37_evaluation_pairs.tsv";

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

# write SP results of good intervals
open(my $fhout2, '>', $eval_pairs) or die "Could not open file ($eval_pairs)!";

# read results.txt
open($fh, '<', $inresults) or die "Could not open file ($inresults)!";

my $lastline = "";
while(my $line = <$fh>) {
	chomp $line;
	$lastline = $line;
	@fields = split('\t', $line);
	if(scalar(@fields) == 5 && ! $excludeIntervals{$fields[0]}) {
		print $fhout2 "$sample\t$fields[1]\t$fields[2]\n";
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

# extract runtime of last line
@fields = split(':', $lastline);
my $runtime = $fields[3];

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
				print $outfh "False positive cis call for pair $fields[1] and $fields[2] in $fields[0] has confidence score $score.\n";
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
				print $outfh "False positive trans call for pair $fields[1] and $fields[2] in $fields[0] has confidence score $score.\n";
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

# print line for table
my $trio_mode = "yes";
if ($trio eq "read-only") {
	$trio_mode = "no"
}

$chr =~ s/chr//;
my $phased = $cis + $trans;
my $solved = $innocuous + $phased;
my $lowConf = $cis_lowConf + $trans_lowConf;
my $fp = $cis_fp + $trans_fp;
my $fp_lowConf = $cis_fp_lowConf + $trans_fp_lowConf;
$runtime = sprintf("%.2f", $runtime/1000);
print "$fam\t$trio_mode\t$chr\t$runtime\t$genes\t$variants\t$potpairs\t$innocuous\t$cis\t$trans\t$phased\t$solved\t$notphased\t$notfound\t$lowConf\t$fp\t$fp_lowConf\n";
