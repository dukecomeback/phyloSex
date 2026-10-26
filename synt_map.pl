#!/usr/bin/perl -w

my $usage=<<EOF;
----------------------------------------------
I hold R's hand to draw the synteny map between chromosomes, improved from synima.pl

Usage: $0 <karyotype.txt> <links.txt> (options)

Options:
-l <label.txt>	set label for karyotypes
-c <color.txt>	set color for karyotypes
-m <marker.txt>	draw markers on karyotypes

>>>>> demo of input >>>>>
karyotype.txt (paint from bottom to top)
human h1:1:5 h2:1:6 h3:1:9
mouse m1:1:8 m2:1:4 m3:1:6

links.txt
h1 3 4 m2 1 2 (color)
h3 3 9 m1 6 2 (color)

marker.txt
gff format, the last column is flexble: "symbol=-;color=-;arrow=0;
<<<<<<<<<<<<<<<<<<<<<<<<<

* cat input.maf |grep '^s' |head -n2 |perl -lane 'print "\$F[1] \$F[1]:1:\$F[5]"' >karyotype.txt

			Du Kang 2022-4-14
----------------------------------------------
EOF

while ($_=shift @ARGV) {
	if (/^-l$/) {
		$labelfile=shift @ARGV;
	}elsif(/^-c$/){
		$colorfile=shift @ARGV;
	}elsif(/^-m$/){
		$markerfile=shift @ARGV;
	}elsif(!/^-/) {
		if (`head -n1 $_ |grep ':[0-9]'`){push @karyo, $_}else{push @links, $_}
	}
}
die $usage if (!@karyo and -t STDIN);

$rhead=<<EOF;
#!/usr/bin/Rscript --vanilla
library(ggplot2)
library(ggforce)
library(scales)
ggplot() +theme_classic() +coord_flip() +
theme(axis.title=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.line.x=element_line(size=0.2), axis.ticks.x=element_line(size=0.2), axis.text.x=element_text(size=7), axis.text.y=element_text(size=12))
EOF
chop $rhead;

if ($colorfile) {
	open COLOR, $colorfile;
	while (<COLOR>) { @F=split; $color{$F[0]}=$F[1] }
	close COLOR;
}

open KAR, "cat @karyo |" or die $!;
while(<KAR>){
	# human h1:1:5 h2:1:6 h3:1:9 ...
	# 1-based coordinates
	next unless /\S/;
	@F=split;
	if ($.==1 and @F>2){	# if multiple chromosomes, remove scales
		$rhead=~s/\.y//g;
		$rhead.="+\ntheme(axis.text.x=element_blank())"
	}
	$breaks.="$.,";
	$labels.="\"$F[0]\",";
	shift @F;
	$suml=0;
	foreach (@F){ /:(\d+):(\d+)/; $suml+=$2-$1+1}; 
	$gap=int($suml/2000*@F);
	$size= @F>=3? 1.5 : 0.5;

	$acum=0;
	foreach (@F){
		($chr,$s,$e)=/(\S+):(\d+):(\d+)/;
		$s{$chr}=$s;
		$e{$chr}=$e;
		$x{$chr}=$.;
		$ymin{$chr}=$s+$acum;
		$ymax{$chr}=$e+$acum;
		$acum=$ymax{$chr}+$gap;
		$color{$chr} //="black";

		$cx.="$x{$chr},";
		$cymin.="$ymin{$chr},";
		$cymax.="$ymax{$chr},";
		$ccolor.="\"$color{$chr}\",";
	}
}
close KAR;
chop ($cx,$cymin,$cymax,$ccolor);
$rest.= "+\ngeom_linerange(aes(x=c($cx),ymin=c($cymin),ymax=c($cymax)),size=$size,color=c($ccolor))";
undef $cx; undef$cymin; undef $cymax; undef $ccolor;

chop ($breaks,$labels);
$rest.="+\nscale_x_continuous(breaks=c($breaks),labels=c($labels))";

open LINK, "cat @links |" or die $!;
$i=0;
while (<LINK>) {
	# h1 3 4 m2 1 2
	$i++;
	@F=split;
	next unless exists $ymin{$F[0]} and exists $ymin{$F[3]};
	my $color= $F[6]? $F[6] : ($F[1]>$F[2] or $F[4]>$F[5])? "gray90" : "darkgoldenrod1";
	$flag= $color=~/^[\d\.]+$/? 1 : 0;
	my ($a,$b,$c,$d)=($F[1]-$s{$F[0]}+$ymin{$F[0]}, $F[2]-$s{$F[0]}+$ymin{$F[0]}, $F[4]-$s{$F[3]}+$ymin{$F[3]}, $F[5]-$s{$F[3]}+$ymin{$F[3]});
	if ($flag) {
		push @link_y, $a, $b, $c, $d;
		push @link_x, $x{$F[0]}, $x{$F[0]}, $x{$F[3]}, $x{$F[3]};
		push @link_group, ($i) x 4;
		push @link_col, ($color) x 4;
	}else{
		$cy{$color}.="$a,$b,$c,$d,";
		$cx{$color}.="$x{$F[0]},$x{$F[0]},$x{$F[3]},$x{$F[3]},";
		$cg{$color}.="$i,$i,$i,$i,";
	}
}
close LINK;

if ($flag) {
	$rlink.="+\ngeom_diagonal_wide(aes(x=c(" . join(',', @link_x) . "), y=c(" . join(',', @link_y) . "), group=c(" . join(',', @link_group) . "), fill=c(" . join(',', @link_col) . ")),alpha=0.75, color=NA)";
	$rlink.="+\nscale_fill_gradient(low=\"blue\", high=\"red\", guide=\"none\")";
}else{
	foreach $color (keys %cy){
		chop ($cy{$color},$cx{$color},$cg{$color});
		$rlink.="+\ngeom_diagonal_wide(aes(c($cx{$color}),c($cy{$color}),group=c($cg{$color})), alpha=0.75, fill=\"$color\")";
	}
	undef %cy; undef %cx; undef %cg;
}

if ($labelfile) {
	my ($cx,$cy,$clabel);
	open LABEL, $labelfile;
	while (<LABEL>) {
		@F=split;
		next unless $s{$F[0]};
		my $x= $x{$F[0]}==1? $x{$F[0]}-0.1 : $x{$F[0]}+0.1;
		my $y=$ymin{$F[0]}+(($ymax{$F[0]}-$ymin{$F[0]})/2);
		$cx.="$x,";
		$cy.="$y,";
		$clabel.="\"$F[1]\",";
	}
	close LABEL;
	chop ($cx,$cy,$clabel);
#	$rest.="+\ngeom_text(aes(c($cx),c($cy)),label=c($clabel),size=1)";
	$rest.="+\ngeom_label(aes(c($cx),c($cy)),label=c($clabel),label.r=unit(0.15,\"lines\"),label.size=NA,size=1.5)";
}

if ($markerfile) {
	my ($cymin,$cymax,$ccolor,$cx,$csymbol,$cty,$ctx);
	open MARKER, $markerfile;
	while (<MARKER>) {
		# symbol=-;color=-
		next if /^\s*$/ or /^#/;
		chomp;
		@F=split;
		next unless exists $s{$F[0]} and $F[3]>=$s{$F[0]} and $F[4]<=$e{$F[0]};

		my ($y1,$y2)=($F[3]-$s{$F[0]}+$ymin{$F[0]},$F[4]-$s{$F[0]}+$ymin{$F[0]});
		my ($ymin,$ymax)= $F[6] eq "+"? ($y1,$y2) : ($y2,$y1);
		$ymax+=50000;
		my $x=$x{$F[0]}+0.05;
		my $color=/color=(\S+?)(;|$)/? $1 : "red"; $color=~s/"//g;
		$cymin.="$ymin,";
		$cymax.="$ymax,";
		$cx.="$x,";
		$ccolor.="\"$color\",";

		if (/symbol=(\S+?)(;|$)/) {
			my $symbol=$1; $symbol=~s/"//g;
			my $ty=($ymin-$ymax)/2 + $ymin;
			my $tx=$x{$F[0]}+0.5;
			$csymbol.="\"$symbol\",";
			$cty.="$ty,";
			$ctx.="$tx,";			
		}
	}
	close MARKER;
	chop ($cymin,$cymax,$ccolor,$cx);
	$rest.="+\ngeom_linerange(aes(x=c($cx),ymin=c($cymin),ymax=c($cymax)),size=1.5,color=c($ccolor))";
	if ($csymbol) {chop ($csymbol,$cty,$ctx); $rest.="+\ngeom_text(aes(c($ctx),c($cty)),angle=45,label=c($csymbol))"}
}

$rest.="+\nscale_y_continuous(expand=c(0,100000), labels=label_number(suffix=\"M\", scale=1e-6),breaks=seq(0,60000000,2000000))";

open R, ">synt_map.r" or die $!;
print R $rhead.$rlink.$rest."\nggsave(filename=\"synt_map.png\",width=7.087,height=1.5,dpi=600)\n";
close R;
`Rscript synt_map.r`;

