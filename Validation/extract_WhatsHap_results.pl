#!/bin/perl -w

# Copyright (C) 2018 the SmartPhase contributors.
# Website: https://github.com/paulhager/smart-phase
# 
# This file is part of the SmartPhase phasing tool.
#
# The SmartPhase phasing tool is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


use strict;

my $in_vcf = $ARGV[0];
my $vars = $ARGV[1];
my $out = $in_vcf;
$out =~ s/.vcf.gz/_phased.tsv/g;

my @header;
my @parts;
my @format;
my @gt_array;
my %phased_gt_data;

my $fh;
if($in_vcf =~ /.gz$/) {
        open($fh, "gunzip -c $in_vcf |") or die "Can't open pipe to $in_vcf!";
} else {
        open($fh, "<", $in_vcf) or die "Can't open VCF ($in_vcf) file!";
}


while(my $line = <$fh>) {
	chomp $line;
	if($line =~ m/^#CHROM/){#Print all Header lines
		@header = split('\t',$line);     
        }

	if(!($line =~ m/^#/)){
		@parts = split('\t',$line);
		if($parts[8] =~ m/PS/) {
			#print $line."\n";
			my $chr = $parts[0];
                        $chr =~ s/chr//g;
                        my $variant = "$chr-$parts[1]-$parts[3]-$parts[4]";

			my $ps_idx = "-1";
			@format = split(':', $parts[8]);
			for(my $i = 0; $i < scalar(@format); $i++) {
				if($format[$i] =~ m/PS/) {
					$ps_idx = $i;
				}
			}
			#print $ps_idx."\n";

			for(my $i = 9; $i < scalar(@parts); $i++) {
				if($parts[$i] =~ m/:/) {
					@gt_array = split(':', $parts[$i]);
					if(scalar(@gt_array) == scalar(@format)) {
						my $ps = $gt_array[$ps_idx];
						if($ps ne ".") {
							my $value = "$gt_array[0]:$ps";
							#print "$header[$i]\t$variant\t$value\n";
							$phased_gt_data{$header[$i]}{$variant} = $value;
						}
					}
				}
			}
		}
	}

}
close($fh);

open($fh, '<', $vars) or die "Could not open $vars!";
open(my $fhout, '>', $out) or die "Could not open $out!";

print $fhout "sample_id\tpos\tpos (#1)\tCompound Het. Prediction\n";

#my $line = <$fh>;

while(my $line = <$fh>) {
	chomp $line;
	@parts = split('\t', $line);
	my $sample = $parts[0];
	my $var1 = $parts[1];
	my $var2 = $parts[2];

	if( exists $phased_gt_data{$sample}{$var1} && exists $phased_gt_data{$sample}{$var2} ) {
		my $phase1 = $phased_gt_data{$sample}{$var1};
		my $phase2 = $phased_gt_data{$sample}{$var2};

		my @fields1 = split(':',$phase1);
		my @fields2 = split(':',$phase2);

		if($fields1[1] eq $fields2[1]) {
			if(($fields1[0] eq "0|1" && $fields2[0] eq "1|0") || ($fields1[0] eq "1|0" && $fields2[0] eq "0|1")) {
				print $fhout "$sample\t$var1\t$var2\t2\n";
			} elsif (($fields1[0] eq "0|1" && $fields2[0] eq "0|1") || ($fields1[0] eq "1|0" && $fields2[0] eq "1|0")) {
				print $fhout "$sample\t$var1\t$var2\t1\n";
			}
		} else {
			 print $fhout "$sample\t$var1\t$var2\t4\n";
		}
	} else {
		print $fhout "$sample\t$var1\t$var2\t4\n";
	}
}
close($fh);
close($fhout);
