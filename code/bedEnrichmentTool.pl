#!/usr/bin/perl -w
##This script use BEDtools to annotate bed file with multiple genomic features, with percentage overlap, within same strand, reverse strand or both.
##it extend BEDtools to allow count from 5'end, 3'end, both border or center of bed. it also allows different upstream and downstream. 
##it will generate summary files that count how many regions are overlapped with each feature(p value using hypergenomitric test), 
##how many features are in this region(p value using hypergenomitric test)

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2012-7-25

#Usege:  perl bedEnrichmentTool.pl [Options] output.txt input_bed/gff_file hg19.chrom.size --anno feature_bed_file_dir1 --anno feature_bed_file_dir2 ...

use strict;
use Getopt::Long;
use File::Basename;

sub usage {

    print "\nUsage:\n";
    print "perl bedEnrichmentTool.pl [Options] outputPrefix input_bed/gff_file genome.chrom.size --anno feature_bed_file_dir1 --anno feature_bed_file_dir2 ...\n\n";

    print " This script use BEDtools to annotate bed file with multiple genomic features.\n";
   

	print "  outputPrefix     Output prefix about enrichment result.\n\n";
	print "  bed/gff_file     Input bed/gff file which requires the enrichment analysis.\n";
	print "  genome.chrom.size     genome chrom size file, used to generate random number in each of the chromosome.\n";
	
	print "  [Options]:\n\n";
    print "  --anno FILE : genomic features' bed files, which are used to annotate input bed files, allow multiple files (required at least one directory)\n\n"; 
    print "  --random_same_chr : random feature generated in the same chr or not (Default: Not enabled)\n\n"; 
    print "  --random_file_num : number of matched random files to generated and calculate SD (Default: 10)\n\n"; 
 print "  --bg FILE : use background file rather than random match (Default: Not enabled)\n\n";	
    print "  --bedtools_path : path to BEDtools/bin/.if not specified in PATH environment variable\n\n";
 	print "  --r_script : complete path to attached r script to generate figure and p value (Default: in the same dir)\n\n";


    exit(1);
}

##default option setting
my @anno_dir = ();
my $random_same_chr = "";
my $bedtools_path = "";
my $random_file_num=10;
my $bg="";
my $r_script = "/jet/home/dnaase/code/mytools/bedEnrichmentTool.R";

GetOptions( "anno_dir=s" => \@anno_dir,
	"bg=s" => \$bg,
			"random_same_chr" => \$random_same_chr,
			"random_file_num=i" => \$random_file_num,
			"bedtools_path=s" => \$bedtools_path,
			"r_script=s" => \$r_script,
			);

usage() if ( scalar(@ARGV) == 0 );

if ( scalar(@ARGV) != 3 ) {
    print "Wrong number of arguments\n";
    usage();
}

if ( scalar(@anno_dir) < 1 ) {
    print "Need to provide at least one directory for the enrichment analysis\n";
    usage();
}

my $outputPrefix = $ARGV[0];
my $input_bed = $ARGV[1];
my $genome = $ARGV[2];

##generate matched random files
my @random_files=();
for(my $i=1;$i<=$random_file_num;$i++){
	if($bg ne ""){next;}
	my $random_tmp=$outputPrefix.".random_$i.tmp.bed";
	if($random_same_chr ne ""){
		runcmd("bedtools shuffle -chrom -i $input_bed -g $genome > $random_tmp");
	}else{
		runcmd("bedtools shuffle -i $input_bed -g $genome > $random_tmp");
	}
	push(@random_files, $random_tmp);
	
}

##write file header
my $output_summary="$outputPrefix.summary.tmp.txt";
open(OUT,">$output_summary") or die "can't write to $output_summary\n";
#print OUT "#file_name\toverlapped\toverlapped_perc\trandom_overlapped_mean\trandom_overlapped_perc_mean\trandom_overlapped_sd\trandom_overlapped_perc_sd\tratio\todds_ratio\tpvalue_chisqure\tpvalue_fisher\n";
#output file name, overlapped, not overlapped, random overlapped mean, random not overlapped mean, random overlapped sd, random not overlapped sd, 

##go over each files in the directory, calculate overlap number/percentage;
foreach my $dir(@anno_dir){
	opendir(DH,"$dir") or die;
	foreach my $file(readdir(DH)){
		#if($file!~/\.bed$/ and $file!~/\.gff$/ and $file!~/\.gtf$/ and $file!~/\.bedgraph$/){next;}
		if($file=~/^\./){next;}	
		my $overlap = `bedtools intersect -wa -u -a $input_bed -b $dir/$file | wc -l`;
		chomp($overlap);
		my $not_overlap = `bedtools intersect -v -a $input_bed -b $dir/$file | wc -l`;
		chomp($not_overlap);
		
		my $overlap_r_mean;
                my $overlap_r_sd;
                my $not_overlap_r_mean;
                my $not_overlap_r_sd;
		if($bg ne ""){
			$overlap_r_mean=`bedtools intersect -wa -u -a $bg -b $dir/$file | wc -l`;
			chomp($overlap_r_mean);
                        $overlap_r_sd=0;
                        $not_overlap_r_mean=`bedtools intersect -v -a $bg -b $dir/$file | wc -l`;
                        chomp($not_overlap_r_mean);
			$not_overlap_r_sd=0;
		}else{
			my @overlap_r_summary;
                	my @not_overlap_r_summary;
			foreach my $random_file(@random_files){
				my $overlap_r = `bedtools intersect -wa -u -a $random_file -b $dir/$file | wc -l`;
				chomp($overlap_r);
				push(@overlap_r_summary, $overlap_r);
				my $not_overlap_r = `bedtools intersect -v -a $random_file -b $dir/$file | wc -l`;
				chomp($not_overlap_r);
				push(@not_overlap_r_summary, $not_overlap_r);
		
			}
			$overlap_r_mean=int(average(@overlap_r_summary));
                	$overlap_r_sd=stdev(@overlap_r_summary);
                	$not_overlap_r_mean=int(average(@not_overlap_r_summary));
                	$not_overlap_r_sd=stdev(@not_overlap_r_summary);

		}

		
		print OUT "$file\t$overlap\t$not_overlap\t$overlap_r_mean\t$overlap_r_sd\t$not_overlap_r_mean\t$not_overlap_r_sd\n";
		print STDERR "$file\t$overlap\t$not_overlap\t$overlap_r_mean\t$overlap_r_sd\t$not_overlap_r_mean\t$not_overlap_r_sd\n";
	}
	closedir(DH);
}
close(OUT);


##input the output file into R, calculate ratio, Odds ratio, p value (chi-square, fisher test), output barplot figure..
my $output_text="$outputPrefix.summary.txt";
my $output_pdf="$outputPrefix.summary.pdf";

runcmd("R --no-restore --no-save --args wd=./ input=$output_summary output=$output_text outpic=$output_pdf < $r_script");

runcmd("unlink $output_summary");
foreach my $random_file(@random_files){
	runcmd("unlink $random_file");
}


sub average{
        my(@data) = @_;
        if (scalar(@data)==0) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@data) {
                $total += $_;
        }
        my $average = $total / scalar(@data);
        return $average;
}
sub stdev{
        my(@data) = @_;
        if(scalar(@data) <= 1){
                return 0;
        }
        my $average = &average(@data);
        my $sqtotal = 0;
        foreach(@data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (scalar(@data)-1)) ** 0.5;
        return $std;
}

sub runcmd{
        my $cmd=shift @_;
        print STDERR "$cmd\n";
        system($cmd)==0 || die "$!\n";
}
