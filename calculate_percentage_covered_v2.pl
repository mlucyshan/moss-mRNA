#!/usr/bin/perl -w
use strict;

#SET UP VARIABLES
my $geneLengthFile=shift @ARGV; #generate this file from the gtf file
my $bedFileFolder=shift @ARGV;
my $mRNAFile=shift @ARGV; 
my $outputFile=shift @ARGV;

my %PPSCoverage;
my %mRNARead;

if ($bedFileFolder !~/\/$/) {
    $bedFileFolder=$bedFileFolder."/";
}

opendir DIR, $bedFileFolder or die "Cannot open directory $bedFileFolder:$!";
my @bedFiles=grep(/\.bed$/,readdir DIR);
closedir DIR;
foreach (@bedFiles){
    print $_,"\n";
}
#Read geneLength into hash
my %geneSpan; 

open(GENESPAN, "<", $geneLengthFile) or die "Cannot open the gene length file.\n";
while (my $geneSpanLine=<GENESPAN>) {
    chomp($geneSpanLine);
    my @geneLine=split(/\t/, $geneSpanLine);
    my $geneIdentifier = (split(/=/, (split /;/, (split(/\t/, $geneSpanLine))[8])[1]))[1];
    $geneSpan{$geneIdentifier}=[$geneLine[3], $geneLine[4], $geneLine[5], $geneLine[4]-$geneLine[3]];
    $PPSCoverage{$geneIdentifier}=[0,0,0]; 
}

sub get_PPS_coverages{
    my %output={};
    return %output; 
}
while (my $bedFile=shift(@bedFiles)) {
    my $bedFilePath=$bedFileFolder.$bedFile; 
    open(BEDFILE, "<", $bedFilePath) or die "Cannot open the bedflile $bedFile.\n";
    while (my $bedFileLine=<BEDFILE>) {
        chomp($bedFileLine);
        my @bedLine=split(/\t/, $bedFileLine);
        my $bedGeneID=(split(/=/, ((split(/;/, $bedLine[6]))[1])))[-1];
        my $PPS_start=$bedLine[1];
        my $PPS_stop=$bedLine[2];
        my $strand=$bedLine[5]; 
        my $timePoints=$bedLine[3];
        if (index($timePoints, "00min") != -1) {
            @{$PPSCoverage{$bedGeneID}}[0]+=($PPS_stop-$PPS_start)/@{$geneSpan{$bedGeneID}}[3];
        }
        if (index($timePoints, "30min") != -1) {
            @{$PPSCoverage{$bedGeneID}}[1]+=($PPS_stop-$PPS_start)/@{$geneSpan{$bedGeneID}}[3];
        }
        if (index($timePoints, "60min") != -1) {
            @{$PPSCoverage{$bedGeneID}}[2]+=($PPS_stop-$PPS_start)/@{$geneSpan{$bedGeneID}}[3];
        } 
    }
    close(BEDFILE); 
}
close(GENESPAN);

my @testValue=@{$PPSCoverage{"Pp3c1_40"}};
foreach my $v (@testValue){
    print $v."\n";
    
}

foreach my $key (keys %PPSCoverage){
    my @values=@{$PPSCoverage{$key}};
    for my $i (0..scalar @values-1){
        if ($values[$i]*1 >1) {
            $values[$i]=1;
        }
    }
}
open(mRNA, "<", $mRNAFile) or die "Cannot open the mRNA read count file.\n";
while (my $mRNALine=<mRNA>) {
    next unless $.>1;
    chomp($mRNALine);
    my @mRNAValues=split(/\t/, $mRNALine); 
    $mRNARead{$mRNAValues[0]}=[$mRNAValues[1],$mRNAValues[2], $mRNAValues[3]]; 
}

close(mRNA); 
open(OUTPUT, ">",$outputFile) or die "Cannot create the output file.\n";
print OUTPUT "TxID\tmin00_mRNA\tmin30_mRNA\tmin60_mRNA\tmin00_PPS\tmin30_PPS\tmin60_PPS\n";
foreach my $key (keys %PPSCoverage){
    my @values=@{$PPSCoverage{$key}};
    my @mRNA=@{$mRNARead{$key}}; 
    print OUTPUT $key;
    for my $l (0..scalar @mRNA-1){
        print OUTPUT "\t$mRNA[$l]";
    }
    for my $i (0..scalar @values-1){
        if ($values[$i]*1 >1) {
            $values[$i]=1;
        }
        print OUTPUT "\t$values[$i]"; 
    }
    print OUTPUT "\n"; 
}
close(OUTPUT); 
exit;
