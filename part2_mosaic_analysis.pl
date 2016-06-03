#~!/usr/bin/perl
#this program will accept case/control files from previous 1st program and filter it by criteria
#filtering criteria: 
#1) base called for only one file
#2) if there are no alternate allele in both files
#3) an alternate alelle by 1 single read
#4) heterozygous in both samples

####input file: <file1.summary.bed> <file2.summary.bed>
####output file:  
#1) *.targets file for each input file
#2) *RESULT file that summarizes the allele comparisons

use warnings;
use strict;
use Getopt::Long;

my $control_hash;
my $test_hash;
my ($cA, $cG, $cT, $cC, $tA, $tG, $tT, $tC);
my ($rA, $rG, $rT, $rC, $r2A, $r2G, $r2T, $r2C);
my (@freqC, @freqT, @readC, @readT);
my @testbreak;
my ($one1, $one2);
my ($filter);
my $t=0;
my ($count);
my ($control_min, $test_min);


my $usage="perl 2_filter_locations_by_criteria.pl <file1.summary.bed> <file2.summary.bed>\n";
$#ARGV>0 || die $usage;

my $file1=$ARGV[0];
my $file2=$ARGV[1];

open (result, ">$file1.$file2.RESULT") || die "cannot open output file\n";

##intersect with target file
system("bedtools intersect -u -a $file1 -b PIK3CA_AKT1.targets > $file1.targets");
system("bedtools intersect -u -a $file2 -b PIK3CA_AKT1.targets > $file2.targets");

#print header
print result "CHR\tLoc\tLoc\t$file1.cov\t$file1.read\t$file1.freq\t$file2.cov\t$file2.read\t$file2.freq\tFILTERED\n";

my @control =`cat $file1.targets`;
my @test =`cat $file2.targets`;
my $len1=@control;
my $len2=@test;
my (@cLine, @tLine);

###storing data to hash
$control_hash=store_to_hash(\@control);
$test_hash=store_to_hash(\@test);

=pod
foreach my $key (keys %{$test_hash}){
        print $key."\t";
}
=cut

##-------------
# main
##-----------

foreach my $key (keys %{$control_hash}){
        print $key."\t";
#       print "in loop! $control_hash->{$key},".$test_hash->{$key}."\n";        
#       if (exists $test_hash->{$key}){

        if (exists  $test_hash->{$key}){
                #print "both have $key existing!\n";
                 ###coverage is >2 
                if ($test_hash->{$key}->{'cov'}>1 && $control_hash->{$key}->{'cov'}>1){
                        print  result "$control_hash->{$key}->{'chr'}\t$key\t$control_hash->{$key}->{'cov'}\t$control_hash->{$key}->{'count'}\t$control_hash->{$key}->{'freq'}\t$test_hash->{$key}->{'cov'}\t$test_hash->{$key}->{'count'}\t$test_hash->{$key}->{'freq'}\t";
                         #breaking down of variables
                         @freqC= split /\:|\,/, $control_hash->{$key}->{'freq'};
                        @freqT= split /\:|\,/, $test_hash->{$key}->{'freq'};
                        $cA= $freqC[1]; $cG= $freqC[3]; $cT=$freqC[5]; $cC=$freqC[7];
                        $tA= $freqT[1]; $tG= $freqT[3]; $tT=$freqT[5]; $tC=$freqT[7];
                        @readC= split /\:|\,/, $control_hash->{$key}->{'count'};
                        @readT= split /\:|\,/, $test_hash->{$key}->{'count'};
                        $rA= $readC[1]; $rG= $readC[3]; $rT=$readC[5]; $rC=$readC[7];
                        $r2A= $readT[1]; $r2G= $readT[3]; $r2T=$readT[5]; $r2C=$readT[7];
                        if ($control_hash->{$key}->{'cov'} < 200){
                                $control_min=1;
                        }elsif ($control_hash->{$key}->{'cov'} < 400){
                                $control_min=2;
                        }
                        elsif ($control_hash->{$key}->{'cov'} < 600){
                                $control_min=3;
                        }
                        elsif ($control_hash->{$key}->{'cov'} > 600){
                                $control_min=4;
                        }
                        if ($test_hash->{$key}->{'cov'} < 200){
                                $test_min=1;
                        }elsif ($test_hash->{$key}->{'cov'} < 400){
                                $test_min=2;
                       }
                         elsif ($test_hash->{$key}->{'cov'} < 600){
                                $test_min=3;
                        }
                        elsif ($test_hash->{$key}->{'cov'} > 600){
                                $test_min=4;
                        }
                        $count=0;

                        print "coverage $test_hash->{$key}->{'cov'}: test min $test_min control: $control_hash->{$key}->{'cov'} $control_min\n";


  ##checks if allele frequency of 1 exists
                        ##print "$control_hash->{$key}->{'freq'}:$test_hash->{$key}->{'freq'}\t";
                        my ($one1) = scalar(@{[($control_hash->{$key}->{'count'}) =~ /:0/gi]});
                        my ($one2) = scalar(@{[($test_hash->{$key}->{'count'}) =~ /:0/gi]});

                        $filter=0;
                        #print "$one1:$one2\t"; 
                         if ($one1==3 && $one2 ==3){
                                #print "both has full allele $cLine[5], $tLine[5]\n";
                                ###checks which one has 1 and checks the same allele for the others
                                ##print "frequency: $cA, $cG, $cT, $cC, \t $tA, $tG, $tT, $tC\n";
                                ###filter2: both are homozygous for the same allele
                                if (($rA>0 && $r2A>0) || ($rG>0 && $r2G>0) || ($rT>0  && $r2T >0) || ($rC>0 && $r2C>0)){
                                        print result "filter: both have same homozygous allele!";
                                         #print "filter: both have same homozygous allele!";
                                }else{
                                        print result "KEEP: different homozygous allele found!";
                                        #print  "KEEP: different homozygous allele found!";
                                }
                        }
                        elsif ($one1==3 && $one2!=3){
                                print result "single_heterozygous:test shows alternative allele: min=$test_min\t";
                                if ($r2A!=0 &&  ($r2A > $test_min)) { $count++;}
                                if ($r2G!=0 && ( $r2G > $test_min)) {$count++;}
                                if ($r2T!=0 && ($r2T > $test_min)) {$count++;}
                                if ($r2C!=0 && ($r2C > $test_min)){ $count++;}
                                if ($count>=2){
                                        print result "KEEP!: New variant found $r2A,$r2G, $r2T, $r2C!!";
                                }
                                else{
                                        print result "filter!! Min read found!!";
                                }
                        }
                        elsif ($one1!=3 && $one2==3){
                                print result "single_heterozygous:control shows alternative allele! min=$control_min\t";
                                if ($rA!=0 && ($rA > $control_min)){ $count++;}
                                if ($rG!=0 && ($rG > $control_min)){$count++;}
                                if ($rT!=0 && ($rT >  $control_min)){$count++;}
                                if ($rC!=0 && ($rC > $control_min)){$count++;}
                                if ($count>=2){
                                        print result "KEEP!: New Variant found!!";
                                }
                                else{
                                        print result "filter!! Min read found!";
                                }

                        }
                        else{  ##both heterozygous sites found
                                print result "both heterozygous found! ";
                                if ($rA > $control_min && $r2A > $test_min){
                                        #my $cAT=2*$rA; my $tAT=2*$r2A; 
                                        my $control_cutA=3*$tA; my $test_cutA=3*$cA;
                                        if ( ($cA > $control_cutA)){ print result "control greater! A $rA, $r2A"; print 2*$r2A; $filter++;}
                                        if ( ($tA > $test_cutA)){ print result "test greater!A $rA, $r2A "; print 2*$rA; $filter++;}
                                }
                                 if ($rC > $control_min && $r2C > $test_min){
                                        my $control_cutC=3*$tC; my $test_cutC=3*$cC;
                                        if ( ($cC > $control_cutC)){ print result "control greater! C $rC, $r2C"; print 2*$r2C; $filter++;}
                                        if ( ($tC > $test_cutC)){ print result "test greater!C $rC, $r2C "; print 2*$rC; $filter++;}
                                 }
                                if ($rG > $control_min && $r2G > $test_min){
                                        my $control_cutG=3*$tG; my $test_cutG=3*$cG;
                                        if ( ($cG > $control_cutG)){ print result "control greater! G $rG, $r2G"; print 2*$r2G; $filter++;}
                                        if ( ($tG > $test_cutG)){ print result "test greater!G $rG, $r2G "; print 2*$rG; $filter++;}
                                }

                                if ($rT > $control_min && $r2T > $test_min){
                                        my $control_cutT=3*$tT; my $test_cutT=3*$cT;
                                        if ( ($cT > $control_cutT)){ print result "control greater! T $rT, $r2T"; print 2*$r2T; $filter++;}
                                        if ( ($tT > $test_cutT)){ print result "test greater!T $rT, $r2T "; print 2*$rT; $filter++;}
                                 }
                                if ($filter >=1){
                                        print result "KEEP! heterozygous sites differ significantly in read counts!";
                                }
                                else {
                                        #print result "filter! does not meet min read/significance difference requirement!";

                                        if ( ($rA<$control_min &&  $r2A> $test_min) || ($r2A<$test_min && $rA > $control_min) ||  ($rG<$control_min &&  $r2G> $test_min) || ($r2G<$test_min && $rG > $control_min) ||($rC<$control_min &&  $r2C> $test_min) || ($r2C<$test_min && $rC > $control_min) || ($rT<$control_min &&  $r2T> $test_min) || ($r2T<$test_min && $rT > $control_min)){
                                                print result "filter! one nucleotide over the min threshold, one not";
                                                $filter=1;
                                        }

                                         if ($filter!=1)
                                        { print result "filter! heterozygous sites are does not reach threshold or significant difference filtered!";}
                                }
                        }



                }
                else {
                        print result "filtered by readcount=1 for at least one file!";
                        #print "filtered by readcount=2 for at least one file!";
                }
                print result "\n";
                }
#}      



}
#-------------------
#subroutine
#------------------
sub store_to_hash{
        my @data=@{$_[0]};
        my %store_hash;
        my @store_temp;

        foreach (@data){
                chomp;
                @store_temp=split /\t/ ,$_;
                $store_hash{$store_temp[1]}->{'chr'}=$store_temp[0];
                $store_hash{$store_temp[1]}->{'cov'}=$store_temp[3];
                $store_hash{$store_temp[1]}->{'count'}=$store_temp[4];
                $store_hash{$store_temp[1]}->{'freq'}=$store_temp[5];
        }
        return \%store_hash;
}

