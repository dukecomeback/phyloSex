#!/bin/bash
function usage(){
echo "
Reciprocal BLAST Hits (RBH) in one button

Usage: $0 a.pep b.pep (gene.pos)

gene.pos	with gene position information provided, I only output those fit the synteny

	Du Kang 2022-5-10
"
}

test $1 || { usage; exit; }

a=`readlink -f $1`
b=`readlink -f $2`
ba=`basename $a`
bb=`basename $b`
test $3 && p=`readlink -f $3`

dir=$ba"_"$bb".rhb"
mkdir $dir; cd $dir

blast.pl $a $b |best >bla
blast.pl $b $a |best |perl -lane 'print "$F[1]\t$F[0]"' >>bla

cut -f1,2 bla |sort -V |uniq -d >rbh.list

test $p && {
	join.pl rbh.list $p |join.pl 1 $p |sort -k3,3V |synteny.pl -egg 3 -gap 5 >rbh.list.synteny
	mv rbh.list.synteny rbh.list
}

cd ../
mv $dir/rbh.list $ba"_"$bb".rbh"
rm -rf $dir
