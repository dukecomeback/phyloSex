#!/usr/bin/perl -w

my $usage=<<EOF;
----------------------------------------------
I transfer mitochondrial genebank to fasta. Beside standard output，I also output a genome.fasta.

Usage: $0 mito.gb >mito.fa

			Du Kang 2025-9-1
----------------------------------------------
EOF

$n=30;
while ($_=shift @ARGV) {
	/^-n$/? $n=shift @ARGV:
	!/^-/? push @file, $_:
	die "Unrecognized parameter $_ !!!";
}
die $usage if (!@file and -t STDIN);

open IN, "cat @file |" or die $!;	# if @file was none, will be read from STDIN
while(<IN>){
	if (/^VERSION\s+(\S+)/) {
		$id=$1;
	}elsif (/(\S+)\s+.*?(\d+\.\.\d+)/) {
		$locus=$2;
		$type{$locus}=$1;
		push @locus, $locus unless $exit{$locus};
		$exit{$locus}=1;
	}elsif (/gene=\"(.*?)\"/) {
		$gene{$locus}=$1;
	}elsif (/product=\"(.*?)\"/) {
		$product{$locus}=$1;
	}elsif (/note=\"(.*?)\"/) {
		$note{$locus}=$1;
	}elsif (/^\s*\d+ (.*)$/) {
		$genome.=$1;
	}
}
close IN;

open OUT, ">>genome.fasta" or die $!;
$genome=~s/ //g;
print OUT ">$id\n$genome\n";
close OUT;

foreach $l (@locus){
	next if $type{$l} eq "source";
	($s,$e)= $l=~/(\d+)\.\.(\d+)/;
	$len=$e-$s+1;
	$seq=substr($genome, $s-1, $len);
	$tag= $gene{$l}? $gene{$l} : $product{$l}? $product{$l} : $type{$l};
	$tag= $note{$l} if $tag eq "s-rRNA";
	$tag= &standard_gene_name($tag);
	$rank{$tag}++;
	$Rtag= $rank{$tag}==1? "" : $rank{$tag};
	print ">$id\t$l\t$type{$l}\t$tag$Rtag\n$seq\n";
}

############################### sub ##########################
sub standard_gene_name {
    my ($name) = @_;
    $name =~ s/\(.*\)//;
    $name =~ s/\s+//g; 
    my $key = $name=~/tRNA/? $name : lc $name;        # 转小写做匹配

    my %map = (
        # --- COX genes ---
        'coi'   => 'COX1',  'cox1'  => 'COX1',  'coxi'   => 'COX1',
        'coii'  => 'COX2',  'cox2'  => 'COX2',  'coxii'  => 'COX2',
        'coiii' => 'COX3',  'cox3'  => 'COX3',  'coxiii' => 'COX3',
	'cytochromecoxidasesubuniti'  => 'COX1',
        'cytochromecoxidasesubunitii'  => 'COX2',
        'cytochromecoxidasesubunitiii' => 'COX3',

        # --- CYTB ---
        'cytb' => 'CYTB', 'cob' => 'CYTB', 'cytochromeb' => 'CYTB',

        # --- ND genes ---
        'nd1' => 'ND1', 'nd2' => 'ND2', 'nd3' => 'ND3', 'nd4' => 'ND4',
        'nd4l' => 'ND4L', 'nd5' => 'ND5', 'nd6' => 'ND6',
        'nadhdehydrogenasesubunit1'  => 'ND1',
        'nadhdehydrogenasesubunit2'  => 'ND2',
        'nadhdehydrogenasesubunit3'  => 'ND3',
        'nadhdehydrogenasesubunit4'  => 'ND4',
        'nadhdehydrogenasesubunit4l' => 'ND4L',
        'nadhdehydrogenasesubunit5'  => 'ND5',
        'nadhdehydrogenasesubunit6'  => 'ND6',

        # --- ATP synthase ---
        'atp6' => 'ATP6', 'atp8' => 'ATP8',
        'atpsynthasef0subunit6' => 'ATP6',
        'atpsynthasef0subunit8' => 'ATP8',
	'atpase6' => 'ATP6',
	'atpase8' => 'ATP8',

        # --- rRNA ---
	'12srrna' => '12S_rRNA',
	'16srrna' => '16S_rRNA',
        '12sribosomalrna' => '12S_rRNA',
        '16sribosomalrna' => '16S_rRNA',
        'l-rrna'   => '16S_rRNA', 'rrnal' => '16S_rRNA',
        'rrnas'    => '12S_rRNA', 'rrn12' => '12S_rRNA', 'rrn16' => '16S_rRNA',

        # --- Control region ---
        'd-loop' => 'D-loop',
        'n_region' => 'N_region',

        # --- tRNA ---
        'trnd' => 'tRNA-Asp', 'trne' => 'tRNA-Glu', 'trnf' => 'tRNA-Phe',
        'trng' => 'tRNA-Gly', 'trnh' => 'tRNA-His', 'trni' => 'tRNA-Ile',
        'trnk' => 'tRNA-Lys', 'trnl' => 'tRNA-Leu', 'trnm' => 'tRNA-Met', 
        'trnp' => 'tRNA-Pro', 'trnq' => 'tRNA-Gln', 'trnr' => 'tRNA-Arg', 
        'trns' => 'tRNA-Ser', 'trnt' => 'tRNA-Thr', 'trnv' => 'tRNA-Val',
        'trny' => 'tRNA-Tyr', 'trna' => 'tRNA-Ala', 'trnc' => 'tRNA-Cys',
	'trnn' => 'tRNA-Asn', 'trnw' => 'tRNA-Trp',
    );

    # 普通基因映射
    return $map{$key} // $name;
}
