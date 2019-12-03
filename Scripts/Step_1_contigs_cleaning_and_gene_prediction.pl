#!/usr/bin/env perl

use strict;
use autodie;
use Bio::Seq;
use File::Spec::Functions;
use File::Which 'which';
use FindBin '$Bin';
use feature 'say';
use Parallel::ForkManager;

# Script to detect circular contigs, nett sequences, and predict genes with mga
# Argument 0 : Fasta file of contigs
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "--h") || ($ARGV[0] eq "-help" )|| ($ARGV[0] eq "--help") || (!defined($ARGV[3])))
{
	print "# Script to detect circular contigs, nett sequences, and predict genes with mga
# Argument 0 : Id du dataset
# Argument 1 : Working dir
# Argument 2 : Fasta file of contigs
# Argument 3 : Threshold on the number of genes 
# Argument 4 : Number of CPUs
# Argument 5 : Gene caller (mga or prodigal)\n";
	die "\n";
}

my $id                = $ARGV[0];
my $tmp_dir           = $ARGV[1];
my $fasta_contigs     = $ARGV[2];
my $th_nb_genes       = $ARGV[3];
my $n_cpus            = $ARGV[4];
my $gene_caller       = $ARGV[5];
my $path_to_prodigal  = which('prodigal');
my $path_to_mga       = which('mga_linux_ia64');
my $path_to_mga_bis   = which('mga');
if ($path_to_mga eq ""){
	if ($path_to_mga_bis eq ""){die("Couldn't find a path for mga, either mga_linux_ia64 (docker) or mga (conda)")}
	else{$path_to_mga=$path_to_mga_bis;}
}

my $in_file           = catfile($tmp_dir, $id . "_nett.fasta");
my $circu_file        = catfile($tmp_dir, $id . "_circu.list");
my $out_special_circu = catfile($tmp_dir, $id . "_contigs_circu_temp.fasta");

# Reading fasta file of the contigs
my $num_seqs = 0;
open my $fa, '<', $fasta_contigs;
my %seq_base;
my $id_seq="";
while(<$fa>){
	$_=~s/\r\n/\n/g; #Cas d'un fichier windows ##AJOUT
	chomp($_);
	if ($_=~/^>(\S*)/){$id_seq=$1; $num_seqs += 1;}
	else{$seq_base{$id_seq}.=$_;}
}
close $fa;

## DETECTION OF CIRCULAR CONTIG AND CLEANING OF THESE CIRCULAR (REMOVE THE MATCHING ENDS)
my $minimum_size=1500;
my %order_contig;
my %length;
my $n1=0;

open my $s1, '>', $in_file;
open my $s2, '>', $circu_file;
for my $id_contig (
    sort {length($seq_base{$b}) <=> length($seq_base{$a})} keys %seq_base){
	$order_contig{$id_contig}=$n1;
	$n1++;
	my $s=$seq_base{$id_contig};
	$length{$id_contig}=length($seq_base{$id_contig});
	# Test its circularity
	my $prefix=substr($seq_base{$id_contig},0,10);
	if ($s=~/(.+)($prefix.*?)$/){
# 		print "on a retrouvé prefix ($prefix) plus loin dans la séquence de $_\n";
		my $sequence=$1;
		my $suffixe=$2;
		my $test=substr($seq_base{$id_contig},0,length($suffixe));
# 		print "$suffixe\n$test\n";
		if ($suffixe eq $test){
# 			print " et il match bien $suffixe, donc c'est un contig circulaire\n";
			my $l=$length{$id_contig};
			$id_contig=$id_contig."-circular";
			$length{$id_contig}=$l;
			print $s2 "$id_contig\t$length{$id_contig}\n";
			$seq_base{$id_contig}=$sequence;
		}
	}
	# Update the length of the contig
	$length{$id_contig}=length($seq_base{$id_contig});
	print $s1 ">$id_contig\n$seq_base{$id_contig}\n";
}
close $s1;
close $s2;

# First, we split $in_file into pieces, one per $n_parts.
# But we have to account for the fact that a user might feed VirSorter
# a fasta file with a single sequence... or request just 1 CPU for some reason.
my $n_parts = $n_cpus;
my $filemax = $n_parts - 1;
if ( $num_seqs <= $n_parts ){$n_parts = $num_seqs; $filemax = $num_seqs - 1;}

my $cmd_split_fasta = "pyfasta split -n $n_parts $in_file";
print "Splitting $in_file into $n_parts pieces...";
`echo $cmd_split_fasta`;
my $out = `$cmd_split_fasta`;
print "\t$out";

if ( $n_parts == 1 ){
    my $cmd_mv_single = "mv $tmp_dir/$id"."_nett.split.fasta $tmp_dir/$id"."_nett.0.fasta";
    `echo $cmd_mv_single`;
    my $out = `$cmd_mv_single`;
}

my $pm = Parallel::ForkManager->new($n_parts);
foreach my $iter (0 .. $filemax) {
    $pm->start and next;
    my $filenum = $iter;
    if ($n_parts > 10) {
        $filenum = sprintf("%02d", $filenum);
    }

    my $process_file = catfile($tmp_dir, $id . "_nett.$filenum.fasta");

    if ($gene_caller eq "prodigal"){
        my $out_gff_part = catfile($tmp_dir, $id . "_nett.$filenum.gff");
        my $prodigal_cmd = "$path_to_prodigal -q -p meta -f gff -i $process_file -o $out_gff_part";
        `echo $prodigal_cmd`;
        $out = `$prodigal_cmd`;
    }
    else{
        my $out_mga_part = catfile($tmp_dir, $id . "_mga.$filenum.predict");
        my $mga_cmd = `$path_to_mga $process_file -m > $out_mga_part`;
        `echo $mga_cmd`;
        $out = `$mga_cmd`;
    }
    $pm->finish;
}
$pm->wait_all_children;

# If prodigal was run, combine the outputs and then convert to mga-formatted output.
my $gff_cmd = catfile($Bin, "gff3_to_mga.pl");
my $out_file= $tmp_dir."/".$id."_mga.predict";

if ($gene_caller eq "prodigal"){
    my $out_gff = catfile($tmp_dir, $id . "_nett.gff");
    say "\nGenerating combined $out_gff using $gene_caller ... ";
    my $cmd_combine_gff = "cat $tmp_dir/$id"."_nett.*.gff > $out_gff; "
        . "rm $tmp_dir/$id"."_nett.*.gff";
    print $cmd_combine_gff;
    $out = `$cmd_combine_gff`;
    say "\t$out";
    
    # Next we use a parser to convert Prodigal's gff3 output to the same output format 
    # as the mga gene caller to make it compatible with the remainder of this script.
    my $cmd_gff_to_mga = "$gff_cmd $out_gff $out_file";
    print $cmd_gff_to_mga;
    $out = `$cmd_gff_to_mga`;
    say "\t$out";
}
# If MGA was run, just combine the outputs.
else{
    say "\nGenerating combined $out_file using $gene_caller ... ";
    my $cmd_combine_mga = "cat $tmp_dir/$id"."_mga.*.predict > $out_file; "
        . "rm $tmp_dir/$id"."_mga.*.predict";
    $out = `$cmd_combine_mga`;
    say "\t$out";
}

# Special prediction for circular contigs if we have some
my $out_file_circu="";
my %circu;
if (-e $circu_file){
	open my $tsv, '<', $circu_file;
	while(<$tsv>){
		chomp($_);
		my @tab=split("\t",$_);
		my $id_c=$tab[0];
		$circu{$id_c}=1;
	}
	close $tsv;
	open my $s3, '>', $out_special_circu;
	my $long=10000; # we cp the 10000 first bases to the end of the contig
	my $n_circu=0;
	foreach(sort {$order_contig{$a} <=> $order_contig{$b} } keys %circu){
		my $id_c=$_;
		my $s=$seq_base{$id_c}.substr($seq_base{$id_c},0,$long);
		print $s3 ">$id_c\n$s\n";
		$n_circu++;
	}
	close $s3;
	$out_file_circu= $tmp_dir."/".$id."_special_circus_mga.predict";
	if ($n_circu>0){
        if ($gene_caller eq "mga") {
            my $mga=`$path_to_mga $out_special_circu -m > $out_file_circu`;
        }
        elsif ($gene_caller eq "prodigal") {
            my $out_file_circu_gff= $tmp_dir."/".$id."_special_circus_mga.gff";
            my $prodigal=`$path_to_prodigal -q -p meta -f gff -i $out_special_circu -o $out_file_circu_gff`;
            my $cmd_gff_to_mga = `$gff_cmd $out_file_circu_gff $out_file_circu`;
        }
	}
	else{
		`touch $out_file_circu`;
	}
}

# Mix 'n match of the two results of gene prediction
my %order_gene;
my $n2=0;
open my $fts, '<', $out_file;
my %predict;
my %type;
my $id_c="";
while(<$fts>){
	chomp($_);
	if($_=~/^# gc/){}
	elsif($_=~/^# self: (.*)/){$type{$id_c}=$1;}
	elsif ($_=~/^# (.*)/){
		$id_c=$1;
		$n2=0;
	}
	else{
		my @tab=split("\t",$_);
		$predict{$id_c}{$tab[0]}=$_;
		if (!defined($order_gene{$id_c}{$tab[0]})){$order_gene{$id_c}{$tab[0]}=$n2;$n2++;}
	}
}
close $fts;
if (-e $circu_file){
	open my $fts_c, '<', $out_file_circu;
	my $tag=0;
	while(<$fts_c>){
		chomp($_);
		if($_=~/^# gc/){}
		elsif($_=~/^# self: (.*)/){$type{$id_c}=$1;}
		elsif ($_=~/^# (.*)/){
			if($tag==1){
				my %to_start;
				# Some ORFs were modified, we clean up
				foreach(sort {$order_gene{$a} <=> $order_gene{$b} } keys %{$predict{$id_c}}){
					my @tab=split("\t",$predict{$id_c}{$_});
					if ($tab[5]!=11){
						# $tab[0] miss start and/or stop codon
						if(($tab[1]<3) || ($tab[2]>($length{$id_c}-3))){
							# And it spans the origin, so we can remove it
							if ($tab[1]<3){
								$to_start{$tab[0]}{"start"}=$tab[1];
								$to_start{$tab[0]}{"stop"}=$tab[2];
							}
							delete($predict{$id_c}{$tab[0]});
						}
						elsif(($tab[2]>997) && ($tab[2]<1001)){ # if we are around the zone of ~ 1000
							foreach(keys %to_start){
								my $total=($length{$id_c}-$tab[1]+1)+($to_start{$_}{"stop"}); 
								if ($total % 3 == 0){
									$tab[2]=$to_start{$_}{"stop"};
									$tab[5]=11;
									my $new_line=join("\t",@tab);
									$predict{$id_c}{$tab[0]}=$new_line;
								}
							}
						}
					}
				}
			}
			$id_c=$1;
			$tag=0;
		}
		else{
			my @tab=split("\t",$_);
			if (defined($predict{$id_c}{$tab[0]})){
				my @tab2=split("\t",$predict{$id_c}{$tab[0]});
				if (($tab2[1]==$tab[1]) && ($tab2[2]==$tab[2])){}# same prediction, we don't change anything
				else{
					if (($tab[1]<$length{$id_c}) && ($tab[2]>$length{$id_c})){
						# we span the origin, we replace the prediction
						$tag=1;
						my $stop=$tab[2]-$length{$id_c};
						$tab[2]=$stop;
						my $new_line=join("\t",@tab);
						$predict{$id_c}{$tab[0]}=$new_line;
					}
				}
			}
			else{
				# we predict a new gene, we keep only if at the start / end of the contig
				if (($tab[1]<$length{$id_c}) && ($tab[2]>$length{$id_c})){
					$tag=1;
					my $stop=$tab[2]-$length{$id_c};
					$tab[2]=$stop;
					my $new_line=join("\t",@tab);
					$predict{$id_c}{$tab[0]}=$new_line;
					$tag=1;
				}
			}
		}
	}
	if($tag==1){
		my %to_start;
		# we changed some things, we clean up
		foreach(sort {$order_gene{$a} <=> $order_gene{$b} } keys %{$predict{$id_c}}){
			my @tab=split("\t",$predict{$id_c}{$_});
			if ($tab[5]!=11){
				if(($tab[1]<3) || ($tab[2]>($length{$id_c}-3))){
					if ($tab[1]<3){
						$to_start{$tab[0]}{"start"}=$tab[1];
						$to_start{$tab[0]}{"stop"}=$tab[2];
					}
					delete($predict{$id_c}{$tab[0]});
				}
				elsif(($tab[2]>997) && ($tab[2]<1001)){
					foreach(keys %to_start){
						my $total=($length{$id_c}-$tab[1]+1)+($to_start{$_}{"stop"}); 
						if ($total % 3 == 0){
							$tab[2]=$to_start{$_}{"stop"};
							$tab[5]=11;
							my $new_line=join("\t",@tab);
							$predict{$id_c}{$tab[0]}=$new_line;
						}
					}
				}
			}
		}
	}
	close $fts_c;
}

## Generation of the final files
## One with all sequences nett and filtered (based on number of genes) - Fasta
## One of the associated gene prediction - MGA-like
## One of the predicted protein sequences - Fasta
my $final_file=$tmp_dir."/".$id."_nett_filtered.fasta";
my $out_final=$tmp_dir."/".$id."_mga_final.predict";
my $prot_file=$tmp_dir."/".$id."_prots.fasta";

open my $fa_s,  '>', $final_file;
open my $out_s, '>', $out_final;
open my $prot_s,'>', $prot_file;
my $n=0;
foreach(sort {$order_contig{$a} <=> $order_contig{$b} } keys %predict){
	$n++;
	if ($n % 10000 == 0){print "$n-ieme contig\n";}
	my $id=$_;
	my @tab_genes=sort {$order_gene{$id}{$a} <=> $order_gene{$id}{$b} } keys %{$predict{$id}};
	my $n_complete_genes=0;
	for (my $i=0;$i<=$#tab_genes;$i++){
		my @tab=split("\t",$predict{$id}{$tab_genes[$i]});
		if ($tab[5]!=11){}
		else{$n_complete_genes++;}
	}
	if ($n_complete_genes<$th_nb_genes){
# 		print "$id is excluded because too short ($n_complete_genes) \n";
	}
	else{
		## We check the first gene and modify it if needed
		my @tab_first=split("\t",$predict{$id}{$tab_genes[0]});
		my @tab_second=split("\t",$predict{$id}{$tab_genes[1]});
		$tab_first[0]=~/gene_(\d*)/;
		my $n_1=$1;
		$tab_second[0]=~/gene_(\d*)/;
		my $n_2=$1;
		if ($n_1>$n_2){
			print "We probably have a circular contig ($id) as the first gene $tab_first[0] is beyond the second gene $tab_second[0] ($n_1>$n_2), so we switch $tab_first[0] ";
			$tab_first[0]="gene_0";
			print "to $tab_first[0]\n";
			$predict{$id}{$tab_genes[0]}=join("\t",@tab_first);
		}
# 		if ($n_complete_genes<$th_nb_genes){print "$id is saved because of its circularity\n";}
# 		else{print "We keep $id = $#tab_genes +1 genes\n";}
		print $out_s ">$id\t$length{$id}\n";
		print $fa_s ">$id\n";
		my $seq_c=$seq_base{$id};
		print $fa_s "$seq_c\n";
		foreach(@tab_genes){
			my @tab=split("\t",$predict{$id}{$_});
			if ($tab[5]!=11){
				# soit on est au début de séquence, soit en toute fin (théoriquement)
				if ($tab[4]!=0){
					if ($tab[3] eq "-"){
						$tab[2]-=$tab[4];
					}
					elsif($tab[3] eq "+"){
						$tab[1]+=$tab[4];
					}
					else{
						print "%%%%%% pblm on a pas de sens pour $id : @tab\n";
					}
				}
				my $new_line=join("\t",@tab);
				$predict{$id}{$_}=$new_line;
			}
			print $out_s "$predict{$id}{$_}\n";
			@tab=split("\t",$predict{$id}{$_});
	                my $name  = $tab[0];
	                my $start = $tab[1];
	                my $stop  = $tab[2];
	                my $sens  = $tab[3];
	                my $frame = $tab[4];
	                my $frag  = "";
			# Regular case (not spanning the origin)
			if ($start<$stop){
				my $long=$stop-$start+1;
				$frag=substr($seq_c,$start-1,$long);
			}
			# Exceptional case, we span the origin
			else{
				my $l1=length($seq_c)-$start+1;
				$frag=substr($seq_c,$start-1,$l1);
				$frag.=substr($seq_c,0,$stop);
			}
			## WE GET THE PREDICTED PROTEIN SEQUENCE
			if ($frag eq ""){
				print "!!!! FRAG IS $frag\n";
			}
			my $seq_bio = Bio::Seq->new(-id => "dummy_id" , -seq =>$frag, -alphabet => 'dna', -verbose => -1);
			# Test to catch the Bio SeqUtils warning
			my @seqs;
			eval{
				@seqs = Bio::SeqUtils->translate_6frames($seq_bio, -verbose => -1);
			};
			if ( $@ ){
				print "We got the error $@\n";
			}
			#my @seqs = Bio::SeqUtils->translate_6frames($seq_bio, -verbose => -1);
			# End of test
			my $cadre=0;
			if ($sens eq "-"){$cadre=3;}
			my $prot=$seqs[$cadre];
			my $prot_sequence=$prot->seq;
			if ($prot_sequence=~/\*$/){
				# we remove the stop codon
				chop($prot_sequence);
			}
			my $id_out=$id."-".$name;
			if (($prot_sequence=~/X{50,}/) || ($prot_sequence=~/F{50,}/) || ($prot_sequence=~/A{50,}/) || ($prot_sequence=~/K{50,}/) || ($prot_sequence=~/P{50,}/)){
				print "we exclude $id_out because there is a pblm with the sequence -> too many succesive X, F, A, K or P\n";
			}
			else{
				print $prot_s ">$id_out\n$prot_sequence\n";
			}
		}
	}
}
close $fa_s;
close $out_s;
close $prot_s;

say "Cleaning up...";
my $rm_files = "rm $tmp_dir/$id"."_nett.*.fasta $tmp_dir/$id"."*.fasta.*";
$out = `$rm_files`;
say "\t$out";
