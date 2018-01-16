#! /usr/bin/perl

use strict;
use warnings;
use FileHandle;

#arg 0 - copy number calls produced by varscan "copyCaller" step
#arg 1 - segments of LOH, which give a reliable centering estimate (optional)

if(!$ARGV[0])
{
    die "USAGE: perl recenter_copynumber.pl cnsegments lohsegments";
}

my %stats = ();

my %chrom_cn = ();
my $neutral_cn_sum = my $neutral_cn_num = 0;

my $cn_file = $ARGV[0];
my $loh_file = $ARGV[1];

if(-e $cn_file){
    %chrom_cn = load_cn($cn_file);

    if(defined($loh_file) && -e $loh_file){
        parse_loh_file($loh_file);
    }
} else {
    die "$cn_file not found\n";
}

my $avg_neutral_cn = 0;
#no LOH regions, fall back to median segment value 
if($neutral_cn_num == 0 && $neutral_cn_sum == 0){
    print STDERR "no LOH regions to consider. Falling back to use the mean genome-wide CN value \n(may be dangerous in polyploid or heavily CN-altered tumors)\n;";
    $avg_neutral_cn = get_genome_mean_cn();    
} else {
    #use the LOH regions
    my $avg_neutral_cn = $neutral_cn_sum / $neutral_cn_num;    
}

open(OUTFILE, ">$ARGV[0].centerinfo") || die("couldn't open outfile");
print OUTFILE "$avg_neutral_cn";
close(OUTFILE);

my $recenter_baseline = $avg_neutral_cn;
print STDERR "recenter_baseline: $recenter_baseline\n";

my $cmd = "";

if($recenter_baseline < 0)
{
    $recenter_baseline = 0 - $recenter_baseline;
    $cmd = "java -cp ~dkoboldt/Software/VarScan net.sf.varscan.VarScan copyCaller $ARGV[0] --output-file $ARGV[0].recentered --recenter-down $recenter_baseline";
    print STDERR "submitting job: $cmd\n";
    system("bsub -q long -J recenter -R\"select[mem>4000] rusage[mem=4000]\" -oo $ARGV[0].recenter.log $cmd");
}
else
{
    $cmd = "java -cp ~dkoboldt/Software/VarScan net.sf.varscan.VarScan copyCaller $ARGV[0] --output-file $ARGV[0].recentered --recenter-up $recenter_baseline";
    print STDERR "submitting job: $cmd\n";
    system("bsub -q long -J recenter -R\"select[mem>4000] rusage[mem=4000]\" -oo $ARGV[0].recenter.log $cmd");
    #system($cmd);
}
exit(0);



##------------------------------------------------
##get median cn value across all segments/probes
##------------------------------------------------

sub get_genome_mean_cn{
    my $sum = 0;
    my $num = 0;
    foreach my $key (keys %chrom_cn) {
        my @a = split("\t",$chrom_cn{$key});
        $num += $a[0];
        $sum += $a[1]*$a[0];
    }
    return($sum/$num);
}


#############################################################
# parse_file - parses the file
#
#############################################################

sub parse_loh_file
{
    my $FileName = shift(@_);

    my $input = new FileHandle ($FileName);
    my $lineCounter = 0;


    my $neutral_regions = 0;
    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;
        my @lineContents = split(/\t/, $line);
        my $numContents = @lineContents;
        $line =~ s/\"//g;
        my ($chrom, $chr_start, $chr_stop, $loh_markers, $loh_mean) = split(/\s+/, $line);

        my $cn_mean = get_avg_cn($chrom, $chr_start, $chr_stop);


        if($loh_mean < 0.05)
        {
            $neutral_regions++;
            my $num_bases = $chr_stop - $chr_start + 1;
            $neutral_cn_sum += ($cn_mean * $num_bases);
            $neutral_cn_num += $num_bases;
        }
        #print STDERR join("--", $chrom, $chr_start, $chr_stop, $loh_markers, $loh_mean, $cn_mean) . "\n";
        #print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $loh_markers, $loh_mean, $cn_mean) . "\n";
    }

    close($input);
}


#############################################################
# parse_file - parses the file
#
#############################################################

sub get_avg_cn
{
    my ($region_chrom, $region_start, $region_stop) = @_;

    my $seg_mean_sum = my $seg_mean_num = 0;

    foreach my $key (keys %chrom_cn) {
        my ($chrom, $chr_start, $chr_stop) = split(/\t/, $key);
        if($chrom eq $region_chrom && $chr_stop >= $region_start && $chr_start <= $region_stop)
        {
            ## We have overlap ##
            my ($num_mark, $seg_mean) = split(/\t/, $chrom_cn{$key});
            $seg_mean_sum += ($num_mark * $seg_mean);
            $seg_mean_num += $num_mark;
        }
    }

    if($seg_mean_num)
    {
        my $avg_seg_mean = $seg_mean_sum / $seg_mean_num;
        return($avg_seg_mean);
    }
    else
    {
        return("NA");
    }
}


#############################################################
# parse_file - parses the file
#
#############################################################

sub load_cn
{
    my $FileName = shift(@_);
    my %cn;

    my $inFh = IO::File->new( $FileName ) || die "can't open file\n";
    while( my $line = $inFh->getline )
    {
        chomp($line);
        my @lineContents = split(/\t/, $line);
        my $numContents = @lineContents;

        #skip header
        next if $line =~ /chrom/;

        $line =~ s/\"//g;
        my ($chrom, $chr_start, $chr_stop, $markers, $normal_depth, $tumor_depth, $mean, $gc_content) = split(/\s+/, $line);
        my $key = join("\t", $chrom, $chr_start, $chr_stop);
        $cn{$key} = join("\t", $markers, $mean);
    }
    close($inFh);

return(%cn);
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub parse_seg_cn
{
    my $FileName = shift(@_);

    my $input = new FileHandle ($FileName);
    my $lineCounter = 0;

    my $cn_sum = my $cn_num = 0;
    while (<$input>)
    {
        chomp;
        my $line = $_;
        $lineCounter++;

        my ($chrom, $chr_start, $chr_stop, $loh_mean, $cn_mean) = split(/\t/, $line);

        if($loh_mean < 0.05)
        {
            my $num_bases = $chr_stop - $chr_start + 1;
            $cn_sum += ($cn_mean * $num_bases);
            $cn_num += $num_bases;
#			print join("\t", $chrom, $chr_start, $chr_stop, $loh_mean, $cn_mean) . "\n";
        }
    }

    close($input);

    if($cn_num)
    {
        ## Calculate average ##
        my $avg_cn = $cn_sum / $cn_num;
        return($avg_cn);
    }

    return("NA");


}
1;
