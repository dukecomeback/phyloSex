#!/usr/bin/perl -w

my $usage=<<EOF;
----------------------------------------------
sliding windows with steps

Usage: $0 input(or stdin) (-w 15 -s 1 -m median)

-m	median/sum

			Du Kang 2026-5-8
----------------------------------------------
EOF

$w=15;
$s=1;
$method="median";
while ($_=shift @ARGV) {
	/^-w$/? $w=shift @ARGV:
	/^-s$/? $s=shift @ARGV:
	/^-m$/? $method=shift @ARGV:
	!/^-/? push @file, $_:
	die "Unrecognized parameter $_ !!!";
}
die $usage if (!@file and -t STDIN);

open IN, "cat @file |" or die $!;	# if @file was none, will be read from STDIN
while(<IN>){
	next unless /\S+/;
	@F=split;
	push @lines, $F[-1];
}
close IN;

$hw=int($w/2);
$start=$hw+1;
$end=@lines-$hw;

for (my $i=$start; $i<=$end; $i+=$s){
	next unless $lines[$i] =~ /\d/;
	my @values=();
	my $s=$i-$hw;
	my $e=$i+$hw; 
	foreach $k ($s..$e){
		$r=$k-1;
		push @values, $lines[$r] if $lines[$r] =~ /\d/;
	}

	if ($method eq "median") {
		@var=sort {$a <=> $b} @values;
		$n=@var;
		$result= @var? $var[int(0.50*$n)] : "null";
	}elsif($method eq "sum"){
		$result=0;
		$result+=$_ for @values;
	}	

	print "$i\t$result\n";
}
