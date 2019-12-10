#!/usr/bin/env perl

use strict;
use autodie;
use warnings;
use Bio::Seq;
use Bio::Tools::GFF;

# Script to convert Prodigal's GFF3 output format to something that looks like
# the output of metagene_annnotator, which is used further on in Step_1_contigs
# Argument 0 : Prodigal GFF3 input
# Argument 1 : Metagene Annotator output

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

my $gffio = Bio::Tools::GFF->new(-file => $infile, -gff_version => 3);

my $feature;
my $segment;
my $i=0;
my $old_name="";
my $new_name;

open my $s1, '>', $outfile;

while($feature = $gffio->next_feature()) {
    $new_name = $feature->{_gsf_seq_id};
    if ( $old_name ne $new_name ){
        print $s1 "# ".$new_name."\n";
        print $s1 "# self: -\n";
    }
    $old_name = $new_name;
    my @id = split("_",$feature->{_gsf_tag_hash}->{ID}->[0]);
    my $strand = "+";
    if ($feature->{_location}->{_strand}=~/-/){
        $strand = "-";
    }
    my $partial = $feature->{_gsf_tag_hash}->{partial}->[0];
    $partial = $partial - 11;
    if ($strand eq "-") {
        if ($partial == -1 || $partial == -10) {
            $partial += 11;
        }
    }
    print $s1 "gene_".$id[1]."\t" . 
        $feature->{_location}->{_start}."\t" . 
        $feature->{_location}->{_end}."\t" . 
        $strand."\t" . 
        $feature->{_gsf_frame}."\t" . 
        sprintf("%02d",abs($partial))."\n";
}
$gffio->close();
close $s1;
