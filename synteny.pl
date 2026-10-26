#!/usr/bin/perl -w

my $usage=<<EOF;
--------------------------------
I keep ortholog pairs that are in synteny, sort the input before using me.

Here is how I do it:
	1) walk through the last column
	2) initiating: when there are \$egg (default 3) genes in a row not melted, initiate the chain
	3) melting: when the chain is interupt by a gap of \$gap (default 15) genes, melt the chain  

Usage: $0 input(or stdin) (-egg 3) (-gap 15) >input.synt
	-egg	how many un-melted genes are needed for iniciating the chain  
	-gap 	how many interrupting genes are needed to melt the chain

Demo of the input: (Forbidden for the 2nd species to be a "polyploid" to the 1st speices, otherwise it leads to a shuffled scenario, and only one chain could be detected)
pla_ENSXMAT00000012342.2        gup_ENSPRET00000017102.1        1:1113941+1115892:10    LG7:8349497:8351489:199
pla_ENSXMAT00000021438.1        gup_ENSPRET00000017143.1        1:1137436:1163575:11    LG7:8373461:8401680:200
pla_ENSXMAT00000009103.2        gup_ENSPRET00000017275.1        1:1170788-1188607:12    LG7:8407280:8420174:201
pla_ENSXMAT00000026565.1        gup_ENSPRET00000017524.1        1:1192672:1218801:13    LG7:8432849:8449566:203

* as you can figure, large segment duplication and tandom duplications will be free from the purge
* watch out for polyploid
                                            		Du Kang 2020-8-8
--------------------------------
EOF

$egg=3;
$gap=15;
while (@ARGV) {
	$_=$ARGV[0];
	if (/^-egg$/) {
		shift @ARGV; $egg=shift @ARGV;
	}elsif(/^-gap$/){
		shift @ARGV; $gap=shift @ARGV;
	}elsif(!/^-/){
		push @file, shift @ARGV;
	}
}
die $usage if (!@file and -t STDIN);
$file="cat @file |";

# I quit to do sort cause I realize it should be dominated it by user himself
#$file= @file? shift @file : "synteny.input.tmp.$$";
#if ($file eq "synteny.input.tmp.$$") {	# if STDIN, put it in a tmp file first
#	unlink $file if -e $file;
#	open TMP, ">>$file" or die $!;
#	while (<STDIN>) {
#		print TMP;
#	}
#	close TMP;
#}

$log="synteny.log";
unlink $log if -e $log;
open LOG, ">>$log" or die $!;
	# output log to here

open STDIN, $file or die $! if @file;
@lines=<STDIN>;
$ini=0;
foreach $i (0..@lines-1){
	$l=$i+1;	# line number
	@col=split /\s+/, $lines[$i];
	($chr,$pos)=$col[-1]=~/(.*):\d+[:+-]\d+:(.*)/;

	if ($ini==1) {	# melting system should be put beofore initiating system so the melt line will go directly into the initiating system
		if ($chr eq $fix_chr and abs($pos-$pre_pos)>=$gap) {
			$ini=0;
			print LOG "line $l: gap size >=$gap, melt the chain: $pos\t$pre_pos\n\n";
		}elsif ($chr ne $fix_chr){
			print LOG "line $l: chain interupted, start the melting mechanism\n";

			$heat=1;
			$reconnect=0;
			foreach $n ($i+1..@lines-1){	# try to reconnect
				$l1=$n+1;
				@col1=split /\s+/, $lines[$n];
				($chr1,$pos1)=$col1[-1]=~/(.*):\d+[:+-]\d+:(.*)/;

				if($chr1 ne $fix_chr){
					$heat++;
				}elsif(abs($pos1-$pre_pos) >=$gap){
					$heat=$gap;
				}else{
					$reconnect=1;
				}

				if ($reconnect==1){
					print LOG "line $l: chain reconnected on line $l1, stop the melting mechanism: $pos1\t$pre_pos\n\n";
					last;
				}elsif ($heat==$gap){
					$ini=0;
					print LOG "line $l: gap $gap achived on line $l1, melt the chain\n\n";
					last;
				}
			}
		}
	}

	if ($ini==0) {	# initiating system, do not use else{} casue melt line need to go into initiating system immediately
		print LOG "line $l: start the initiating processing\n";

		$heat=0;	# heat for melting
		$force=1,	# force for initiating
		$fix_chr=$chr;
		$pre_pos1=$pos;
		foreach $n ($i+1..@lines-1){	# try to spark
			$l1=$n+1;
			@col1=split /\s+/, $lines[$n];
			($chr1,$pos1)=$col1[-1]=~/(.*):\d+[:+-]\d+:(.*)/;
			
			if($chr1 ne $fix_chr){
				$heat++;
			}elsif(abs($pos1-$pre_pos1) >=$gap){
				$heat=$gap;
			}else{
				$force++;
			}

			if ($force==$egg){
				$ini=1;
				print LOG "line $l: initiated\n\n";
				last;
			} elsif ($heat==$gap) {
				print LOG "line $l: melted the initiating processing at line $l1, initiating failed\n\n";
				last;
			}
			
			$pre_pos1=$pos1 if $chr1 eq $fix_chr;
		}

	}

	print $lines[$i] if $ini==1 and $chr eq $fix_chr;
	$pre_pos=$pos if $chr eq $fix_chr;
}

#unlink "synteny.input.tmp.$$" if -e "synteny.input.tmp.$$";

