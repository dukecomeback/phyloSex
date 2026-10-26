#!/bin/bash
function usage(){
echo "
I use barrnap to identify the ITS2 (internal transcribed spacer 2 that between 5.8S and 28S, could be used for phylogeny construction) sequence from a genome assembly (then go to http://its2.bioapps.biozentrum.uni-wuerzburg.de the old version for sophisticated annotation)

At first we plan to annotate 5.8S and 28S and then cut out the sequence in between. However, the 28S predictions are awful given its 4kb length. Good news is that the 154bp 5.8S are always delightful annotated.
After check the information in NCBI, we found the ITS2 are always ~670bp long in X. maculatus, we hence take the adventage of annotion of 5.8S, and cut its 1kb flanking sequence extending tarwards to the 28S direction. Hence we are sure that the ITS2 and part of 28S is inside the seqeunces we retrieved.

Usage: $0 assembly 

	Du Kang 2022-08-25
"
}

if [ ! -n "$1" ]; then
	usage
	exit
fi

assembly=`readlink -f $1`
basename=`basename $assembly`

barrnap -q -k euk -t 32 $assembly |perl -lape 's/Name=//;s/product=.*//' |perl -lane 'next unless /28s_|5_8s_/i; $m=int($F[3]+($F[4]-$F[3])/2); print "$F[0]$F[6]\t$m:$F[3]:$F[4]:$F[-1]"' |sort -k1,1V -k2,2V |dim_up |perl -lane '($scaf,$str)=$F[0]=~/(.*)([+-])/;if($str eq "+"){while(/:(\d+):5_8S_.*?\t\d+:\d+:\d+:28S_/g){$s=$1; $e=$s+1000; print "$scaf\trnammer\tITS2\t$s\t$e\t\.\t$str\tITS2candidate"}}else{while(/:28S_.*?\t\d+:(\d+):\d+:5_8S_/g){$e=$1; $s=$e-1000; print "$scaf\trnammer\tITS2\t$s\t$e\t\.\t$str\tITS2candidate"}}' |gff2fa.pl $assembly -fetch "" -id "" -extend 100 >$basename.its2_candidate.fa




: << EOF
ITSx -i its2_candi.fa -o out --preserve T --partial 10 --anchor 28
	# then check, maybe remove the bad result based on sequence length

fasta_spliter.pl out.ITS2.fasta
ls hehe.* |perl -lane 'print "RNAfold <$_ >$_.vienna"' |bsub -w
cat hehe*.vienna >out.ITS2.fasta.vienna
clean; rm -rf hehe.*

rna2pseudoAA.pl out.ITS2.fasta.vienna >out.ITS2.fasta.vienna.pseudoAA
echo "clustalw -INFILE=out.ITS2.fasta.vienna.pseudoAA -TYPE=PROTEIN -OUTFILE=out.ITS2.fasta.filter.vienna.pseudoAA.clustalW -OUTPUT=FASTA -ALIGN" |bsub

rna2pseudoAA.pl -bak out.ITS2.fasta.filter.vienna.pseudoAA.aln >out.ITS2.fasta.filter.vienna.pseudoAA.aln.bak2rna
EOF
