#~!/usr/bin/perl -w
##extracts an individual, creates a file for that individual for the locations in vcf file
#does not filter for anything
use strict;
use Getopt::Long;
use warnings;

die "MISSING INPUT OPTIONSE!\nperl format_vcf_file_readable.pl <vcf file> <id file> <output folder>\n" if @ARGV !=3;

#print "input files: $ARGV[0], $ARGV[1], $ARGV[2]\n";
my $vcf=$ARGV[0];
open (id_input_file, $ARGV[1]) || die "ID file cannot be opened\n";
my $ind;
my $output;
my $output_folder=$ARGV[2];


#------------
# Begin MAIN 
#------------

while (<id_input_file>){
	
	chomp;
	$ind=$_;
	$output=$ind;
	print "processing $ind\n";
	system("GATK SelectVariants -R /data/ClinSeq/reference/human_g1k_v37.fasta \\
	-V /data/hongcs2/DAVID_REDO/analysis_1929/$vcf \\
	-L /data/hongcs2/DAVID_REDO/analysis_1929/DMET_indel_location.convert.loc.interval_list \\
	-o $output_folder/$ind.recode.vcf");
	
	#reformats individual vcf file
	open (SNV, "$output_folder/$output.recode.vcf") || die "$output.SNV.recode.vcf cannot open file\n";
	extract("$output_folder/$output.recode.vcf");
}

#------------
# End MAIN
#------------


#------------
#Subroutine
#------------

sub extract {
	my $filename=shift;
	print "$filename!!\n";
	open (file, "$filename") || die "$filename cannot open file\n";

	while (<file>){
		if (/CHROM/){last;}
	}
	##breaks down each line and reformats
	open (out, ">$filename.reformatted") || die "cannot open output file $filename.reformatted\n";
	print out "CHROM\tPOS\tID\tREF\tALT\tGENOTYPE\tCOVERAGE\tQUALITY\n";
	my %hash;
	my @line;
	while (<file>){
		chomp;
		#print $line;
		@line=split /\t/, $_;	
		
		print out "$line[0]\t$line[1]\t$ind\t$line[3]\t$line[4]\t";

		my @format=split/:/,$line[8];
		my @values=split /:/,$line[9];

		my $len=@format;

		for my $i (0..$len-1){
			#print "$format[$i] $values[$i] ";
			if ($format[$i] eq "RGQ"){ $hash{'GQ'}=$values[$i];}
			else{
				$hash{$format[$i]}=$values[$i];
			}
		}
		if ($hash{'GT'} eq "./."){
			$hash{'DP'} = 0;
			$hash{'GQ'} =0;
		}
		print out "$hash{'GT'}\t$hash{'DP'}\t$hash{'GQ'}\n";
	}
}
