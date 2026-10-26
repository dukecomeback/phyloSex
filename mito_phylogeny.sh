#!/bin/bash
function usage(){
echo "
Build phylongey tree from mitochondral GeneBanks

Usage: $0 mito.gbs/

	Du Kang 2025-8-27
"
}

test $1 || { usage; exit; }

codon="/publicfs01/fs1-m8/home/m8s000630/micromamba/share/wise2/wisecfg/codon_mito.table"
	# replace this PATH in genewise.pl

cpu=32
while [ $1 ]; do
	case $1 in
		-t) shift; cpu=$1 ;;
		*) P+=($1) ;;
	esac
	shift
done
gb=`readlink -f ${P[0]}`
base=`basename $gb`

mkdir $base.mitoPhy; cd $base.mitoPhy

mkdir fasta
ls $gb/* |perl -lane '/.*\/(\S+)/; print "mitogb2fa.pl $_ >fasta/$1.fa"' |sh

########################### do protein coding genes
mkdir ref pcg.qry pcg.out

cat genome.fasta |perl -lane 'open OUT, ">>ref/$1" if />(\S+)/; print OUT $_'
cat fasta/*.fa |fatab |grep CDS |perl -lane 'print ">$F[3]___$F[0]\n$F[4]"' |cds2pep.pl -mt |fatab |grep -v "\*[A-Za-z]" |perl -lane '($id)=/(\S+?)___/; $flag=$F[1]=~/^M/? 1 : 0; $l=length $F[1]; print "$id\t$flag\t$l\t$F[1]"' |sort -k1,1V -k2,2nr -k3,3nr |cut -f1,4 |best |fatab bak |perl -lane 'open OUT, ">>pcg.qry/$1" if />(\S+)/; print OUT $_'

for qry in pcg.qry/*; do
for ref in ref/*; do
genewise.pl $qry $ref -mt |perl -lane 'if(/\tgene\t/){print ">$F[0]"}elsif(/coding sequence = \[(\S+)\]/){print $1}' >>pcg.out/$(basename $qry)
done
done

cd pcg.out/
ls |perl -lane 'print "macse -prog alignSequences -seq $_ -gc_def 2"' |waitn
perl -i -lape 's/\!/-/g unless /^>/' *_NT
ls *_NT |perl -lane 'print "Gblockc.pl $_ >$_.gb"' |waitn
cat *gb >../pcg.cds.mus.gb
cd ../

fasta_concatenater.pl pcg.cds.mus.gb >pcg.cds.mus.gb.con

export l=`falen -aln pcg.cds.mus.gb.con |cut -f2`
perl -e 'print "DNA, p1=1-$ENV{l}\\3,2-$ENV{l}\\3\n";print "DNA, p2=3-$ENV{l}\\3\n"' >partition.pcg
	# divide into two partion: one for the 1st and 2nd position of codons and one for the 3rd position.

echo "raxmlHPC-PTHREADS-SSE3 -s pcg.cds.mus.gb.con -m GTRGAMMA -p 12345 -q partition.pcg -f a -x 12345 -# 1000 -T $cpu -n pcg" |bsub -p
	# https://cme.h-its.org/exelixis/web/software/raxml/hands_on.html


######################## do pcg and rna together
mkdir rna.out
cat fasta/*.fa |fatab |grep RNA |perl -lane 'print ">$F[0]\t$F[3]\n$F[4]"' |perl -lane 'open OUT, ">>rna.out/$F[1]" if />/; print OUT $_'
cd rna.out
ls |perl -lane 'print "muscle -in $_ -out $_.mus ; Gblockc.pl -sin $_.mus >$_.mus.gb"' |waitn
cat *gb >../rna.mus.gb
cd ../

ls pcg.out/*gb rna.out/*gb |perl -lane '$l=`falen $_ -aln`; print "$_\t$l"' |perl -lane 'next unless /\S+/; $s=$flag+1; $e=$s+$F[-1]-1; print "DNA, p$.=$s\-$e\t$F[0]"; $flag=$e' >partition.pcg_rna
cat pcg.cds.mus.gb rna.mus.gb |fasta_concatenater.pl >pcg_rna.mus.gb.con

echo "raxmlHPC-PTHREADS-SSE3 -s pcg_rna.mus.gb.con -m GTRGAMMA -p 12345 -q partition.pcg_rna -f a -x 12345 -# 1000 -T $cpu -n pcg_rna" |bsub -p

