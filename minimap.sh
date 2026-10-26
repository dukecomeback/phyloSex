#!/bin/bash
function usage(){
echo "
minimap aligning (intra/cross species) followed by chain and net

Usage: $0 reference.fa query.fa (asm5/asm20/map-hifi) (-t 32) (-r) 

-t 		number of threads, default 32
-r		raw output of minimap without further chain & net or sort. This will ouput a .paf and the dotplot in .pdf
asm5		for sequence divergence below 5%, good for intra-species alignment
asm20		for divergence not more than 20%, good for cross-species alignment (default)
map-hifi	Align PacBio HiFi reads

* reference will be the first line of output.maf, which is also required by multiz
* https://github.com/hillerlab/GenomeAlignmentTools
* https://darencard.net/blog/2019-11-01-whole-genome-alignment-tutorial/
* https://github.com/lh3/minimap2/blob/master/cookbook.md
* last-dotplot --bed1= can be used to mark out contig joint point (Ns); Multiple \"-1 PATTERN\" can be used to zoom in.

	Du Kang 2021-8-3
"
}

test $1 || { usage; exit; }

cpu=32
asm="asm20"
while [ $1 ]; do
	case $1 in
		-r) raw=1 ;;
		-t) shift; cpu=$1 ;;
		asm5) asm="asm5" ;;
		asm20) asm="asm20" ;;
		map-hifi) asm="map-hifi" ;;
		*) P+=($1) ;;
	esac
	shift
done
file1=${P[0]}
file2=${P[1]}
p=${P[2]:-"default"}

t=`basename $file1`
q=`basename $file2`
o=$t"_"$q

target=`readlink -f $file1`
query=`readlink -f $file2`

tmp="tmp.minimap.sh.$$"

if [ -n "$raw" ]
	then
		echo "raw output without further chain&net !!!"
		# minimap2 $target $query -cx $asm --cs=long -t $cpu |paftools.js view -f maf - >$o.maf
		minimap2 $target $query -cx $asm --cs=long -t $cpu >$o.paf
		cat $o.paf |paftools.js view -f maf - >$o.maf
		# paf2dotplot.r $o.paf

	else
		mkdir $tmp
		cd $tmp
		
		minimap2 $target $query -cx $asm --cs=long -t $cpu |paftools.js view -f maf - |maf-convert psl >minimap.psl
			# align - paf to maf to psl
		
		axtChain -linearGap=loose -psl minimap.psl -faQ $target -faT $query minimap.chain & 
			# psl to chain
		
		faToTwoBit $target target.2bit &
		faToTwoBit $query query.2bit &
		falen $target >target.size &
		falen $query >query.size &
		wait
	
		export hoxd55q=`readlink -f ~/usr/local/GenomeAlignmentTools/example/HoxD55.q`	
		patchChain.perl minimap.chain target.2bit query.2bit target.size query.size -lastzParameters "--format=axt K=1500 L=2500 M=0 T=0 W=5 Q=$hoxd55q"
		cat jobList |parallel -j $cpu
		cat doPatchChain/*.psl >>minimap.psl
		sort -k10,10V -k12,12n  minimap.psl -o minimap.psl
		axtChain -linearGap=loose -psl minimap.psl -faQ $target -faT $query minimap.chain
			# patchChain.perl improve the chains by adding flanked loci aligning
		
		RepeatFiller.py -c minimap.chain -T2 target.2bit -Q2 query.2bit >tmp.chain
		mv tmp.chain minimap.chain
			# RepeatFiller.py improve the chains by adding TE alignments
		
		chainCleaner minimap.chain target.2bit query.2bit -tSizes=target.size -qSizes=query.size tmp.chain removedSuspects.bed -linearGap=loose
		mv tmp.chain minimap.chain
			# remove obscure loci alignments
		
		chainPreNet minimap.chain target.size query.size tmp.chain
		mv tmp.chain minimap.chain
			# remove chains that do not have a chance of being netted
		
		chainNet -rescore minimap.chain target.size query.size target.prenet query.prenet -tNibDir=target.2bit -qNibDir=query.2bit -linearGap=loose
			# Hiller's version of netting
		
		netSyntenic target.prenet target.net
			# add synteny information
		
		NetFilterNonNested.perl -doUCSCSynFilter -keepSynNetsWithScore 5000 -keepInvNetsWithScore 5000 target.net >target.filter.net
		mv target.filter.net target.net
			# filter net
		
		netToAxt target.net minimap.chain target.2bit query.2bit out.axt
		axtToMaf out.axt target.size query.size minimap.maf
			# convert net to maf

		maf-sort minimap.maf >../$o.maf
			# maf-sort conda安装last即可获得
			# flat '^a|^#' minimap.maf |sort -k1,1Vr -k4,4n |flat bak >../$o.maf
			# sort according to the reference
		
		cd ../
		last-dotplot $o.maf $o.maf.gif --rot1=v --rot2=h
			# plot
fi

/bin/rm -rf $tmp
