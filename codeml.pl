#!/usr/bin/perl -w

use Parallel::ForkManager; 

my $usage=<<EOF;
--------------------------------
Usage: $0 -a align.fa -t unroot.tre -m choose_a_model_from_below (-fixw 1) (-wantmlc) (-mt)

Models:
	pairwise		pairwise model, with omega unfixed
	pairwise_null		pairwise model, with omega fixed to 1 (default)
	pairwise_test		run pairwise and pairwise_null, then do a likelihood ratio test

	branch_free		branch model, with omega estimated for each branch <free-ratio model>
				This model is very parameter rich and cannot be effectively compared to a suitable null model and so branch-specific rates must be considered approximations.

	branch_m0		branch model, with one estimated omega shared across all branches <one-ratio model>
	branch 			branch model, with omega estimated for each branch-group <normal two groups: background and foreground, two-ratio model> (tag #1 for the foreground group in tree file)
	branch_null		branch model, with omega fixed to 1 for foreground brach, estimated and shared for background branch
	branch_test		run branch_m0, branch_null and branch, then compare branch, by likelihood ratio test, to branch_m0 ana branch_null

	branchsite 		branch_site model, with omega varies across codon sites and branch groups (tag #1 for the foreground branch in tree file) 
	branchsite_null		branch_site model, with omega fixed to 1 for foreground branch, estimated and shared for backgroupd branch
	branchsite_test		run branchsite and branchsite_null, then do a likelihood ratio test

Others:
	-a 			I need the alignment in fasta format !!!
	-t			unrooted tree
	-m 			default model: branch_free
	-fixw			set the fix value for omega in the omega-fixed model, default 1
	-wantmlc		by default I won't keep the result file mlc, unless you give me this
	-mt			mitochondrial CDS

about the output:
	I always reveal the model type
	If there are more than one model type, then it's a test output, the p value[s] (2 for branch_test) comes in the end
	If you see numbers separated by comma, that's the site position (!: <0.01; .: <0.05)
	Normally the order for a model would be "omega lnL dS"; if more then one omega exist, then the order of omega would be "background #1 #2 ..."
	If there is a table, then read it as "branch t N S dN/dS dN dS N*dN S*dS"
	If there is a alignment with InK, positive selected sites are marked with "!" (**, lever 0.01) or "." (*, lever 0.05) in the InK. Be attention
	If a sequence name is output, I only keep the first three letters from the sequence name. So I suggest you to put a 3-letter species index before sequence name, like "hum_ENSG00000139618"

for more information:
	https://embnet.vital-it.ch/CoursEMBnet/PagesPHYL07/Exercises/day2/day2.html
	https://www.douban.com/note/526755885/
	https://www.omicsclass.com/article/37 (for site model)

                                         						 	  							Du Kang 2019-3-12
--------------------------------
EOF

################################ parameters
@ARGV or die $usage;

$model="branch_free";
$fixw=1;
$mt=0;
while ($_=shift @ARGV) {
	/^-a$/? ($align=shift @ARGV and $align=`readlink -f $align`):
	/^-t$/? ($tree=shift @ARGV and $tree=`readlink -f $tree`):
	/^-m$/? $model=shift @ARGV:
	/^-fixw$/? $fixw=shift @ARGV:
	/^-wantmlc$/? $keepmlc=1:
	/^-mt$/? $mt=1:
	die "Unrecognized parameter $_\n";
}

$align or die "Please feed me the alignment\n";

@options=qw/pairwise pairwise_null pairwise_test branch_m0 branch_free branch branch_null branch_test branchsite branchsite_null branchsite_test/;
grep { $_ eq $model } @options or die "\nSorry we don't have this model, or is there a typo?\n\n";

if (! $tree){
	$model =~ /pairwise/ or die "\nYou may forgot an unroot tree for our model $model\nor you may forgot to set the model to be pairwise\n\n";
}else{
	$notag=`cat $tree`!~/#1|\$1/;
	die "\nYou may forgot to tag the tree for our model $model\n\n" if $notag and ($model=~/branchsite/ or $model eq "branch_test" or $model eq "branch" or $model eq "branch_null");	
}

############################### if the model is branchsite or branchsite_test, readin and translate the cds sequence
if ($model eq "branchsite" or $model eq "branchsite_test") {
	open FA, $align or die $!;
	while (<FA>) {
		chomp;
		if (/>/) {
			s/>//;
			@_=split;
			$name=$_[0];
		}else{
			s/\s//g;
			$cds{$name}.=$_;
		}
	}
	foreach $key (keys %cds){
		$pep{$key}=&cds2aa($cds{$key});
	}
}

############################### run codeml and parse the output
$basename = $align=~/.*\/(.*)/? $1 : $align;
$outfile = $basename.".$model";
print "$outfile already exists, I'm adding output into this file\n" if -e $outfile; 
open OUT, ">>$outfile" or die $!;

if ($model eq "pairwise" or $model eq "pairwise_null"){
	($mlcfolder,$lnL) = &run_codeml($model);
	`cp $mlcfolder/mlc $outfile.mlc` if $keepmlc;
#	`grep 'dN\\/dS=' $mlcfolder/mlc` =~ /dN\/dS=.*?(\d.*?)\s.*=(.*)/;
#	$w = $1;
#	$dS = $2;
	$otherinfo=`grep 'dN\\/dS=' $mlcfolder/mlc`;
	$otherinfo=~s/(=|=\s+)/= /g;
	$otherinfo=~s/\s+=/=/g;
#	print OUT ">$basename\t$model\t$w\t$lnL\t$dS\t$otherinfo\n";
	print OUT ">$basename\t$model\t$lnL\t$otherinfo";

}elsif ($model eq "branch_m0"){
	($mlcfolder,$lnL) = &run_codeml($model);
	`cp $mlcfolder/mlc $outfile.mlc` if $keepmlc;
	@_ = split /\s+/, `grep "\(dN\\/dS\)" $mlcfolder/mlc`;
	$w = $_[-1];
	print OUT ">$basename\t$model\t$w\t$lnL\n";

}elsif ($model eq "branch" or $model eq "branch_null"){
	($mlcfolder,$lnL) = &run_codeml($model);
	`cp $mlcfolder/mlc $outfile.mlc` if $keepmlc;
	`grep "\(dN\\/dS\)" $mlcfolder/mlc` =~ /for branches:\s*(\d.*)/;
	$w = $1;
	print OUT ">$basename\t$model\t$w\t$lnL\n";

}elsif ($model eq "branch_free"){
	($mlcfolder,$lnL) = &run_codeml($model);
	`cp $mlcfolder/mlc $outfile.mlc` if $keepmlc;
	print OUT ">$basename\t$model\t$lnL\n";

	# print the "table" in the mlc
	open MLC, "$mlcfolder/mlc" or die $!;
	while (<MLC>) {
		@_=split;
		if (/^#/) {
			$_[0]=~/#(\d+):/;
			$tag=$1;
			$spe{$tag}=$_[-1];	
		}elsif (/\d\.\.\d/ && @_==9){
			$first=shift @_;
			$first=~/\d+\.\.(\d+)/;
			$tag=$1;
			if (exists $spe{$tag}) {
#				$spe{$tag}=~/^(...)/;
#				$head=$1;
				$head=$spe{$tag};
				printf OUT "%-5s", "#$head";
				print OUT "\t@_";
				print OUT "\n";
			}
		}
	}

}elsif ($model eq "branchsite_null"){
	($mlcfolder,$lnL) = &run_codeml($model);
	`cp $mlcfolder/mlc $outfile.mlc` if $keepmlc;
	print OUT ">$basename\t$model\t$lnL\n";

}elsif ($model eq "branchsite"){
	($mlcfolder,$lnL) = &run_codeml($model);
	`cp $mlcfolder/mlc $outfile.mlc` if $keepmlc;
	print OUT ">$basename\t$model\t$lnL";

	$mlc=`cat $mlcfolder/mlc`;
	$mlc=~/Bayes Empirical Bayes.*:\n(.*)\n\s+The grid/s;
	$positive=$1;
	if ($positive=~/\*/s) {
		@_=split /\s+/, `grep "ls =" $mlcfolder/mlc`;
		$aanum=$_[-1];

		$InK="0" x $aanum;	# this $InK string is going to show you positive sites with "."(*) or "!"(**)

		@lines= split /\n/, $positive;
		foreach $_ (@lines){
			next unless /\*/;
			@_=split;
			$mark = /\*\*/? "!" : ".";
			substr($InK, $_[0]-1, 1)=$mark;
			$reportpos = /\*\*/? "$_[0]!," : "$_[0].,";
		}
		print OUT "\t$reportpos";
		printf OUT "\n%-5s", "#InP";
		print OUT "\t$InK\n";

		foreach $key (keys %pep){
			$key=~/^(...)/;
			$head=$1;
			printf OUT "%-5s", "#$head";
			print OUT "\t$pep{$key}\n";
		}
	}else{
		print OUT "\n# No positive selected site is detected\n";
	}

}elsif ($model eq "pairwise_test") {
	my $pm = new Parallel::ForkManager(2);

	$pm -> run_on_finish( sub {     # reep son threads' $report
                my @recall = @{$_[-1]};
                $recall[0]=~/.*_(pairwise.*?)\.codemlstudio/;
                $mlcfolder{$1}=$recall[0];
                $lnL{$1}=$recall[1];
        });

        foreach $key (qw/pairwise_null pairwise/){    # launch threads
                $pm -> start and next;
                @report=&run_codeml($key);
                `cp $report[0]/mlc $outfile\_$key.mlc` if $keepmlc;
                $pm -> finish(0,\@report);
        }

        $pm -> wait_all_children;  #recycle threads

	foreach $key (qw/pairwise_null pairwise/){
		`grep 'dN\\/dS=' $mlcfolder{$key}/mlc` =~ /dN\/dS=.*?(\d.*?)\s.*=(.*)/;
        	$w{$key} = $1;
        	$dS{$key} = $2;
	}

	$p = &chi2($lnL{pairwise}, $lnL{pairwise_null});

	print OUT ">$basename\tpairwise_null\t$w{pairwise_null}\t$lnL{pairwise_null}\t$dS{pairwise_null}\tpairwise\t$w{pairwise}\t$lnL{pairwise}\t$dS{pairwise}\t$p\n";

}elsif ($model eq "branch_test") {
	my $pm = new Parallel::ForkManager(4);

	$pm -> run_on_finish( sub {	# reep son threads' $report
		my @recall = @{$_[-1]};
		$recall[0]=~/.*_(branch.*?)\.codemlstudio/;
		$mlcfolder{$1}=$recall[0];
		$lnL{$1}=$recall[1];
 	});

	foreach $key (qw/branch_m0 branch_null branch/){	# launch threads
		$pm -> start and next;
		@report=&run_codeml($key);
		`cp $report[0]/mlc $outfile\_$key.mlc` if $keepmlc;
		$pm -> finish(0,\@report);
	}

	$pm -> wait_all_children;  #recycle threads

	foreach $key (qw/branch_m0 branch_null branch/){
		`grep "\(dN\\/dS\)" $mlcfolder{$key}/mlc` =~ /.*?\s*(\d.*)/;
        	$w{$key} = $1;
	}
	
	$p1 = &chi2($lnL{branch_m0}, $lnL{branch});
	$p2 = &chi2($lnL{branch_null}, $lnL{branch});

	print OUT ">$basename\tbranch_m0\t$w{branch_m0}\t$lnL{branch_m0}\tbranch_null\t$w{branch_null}\t$lnL{branch_null}\tbranch\t$w{branch}\t$lnL{branch}\t$p1\t$p2\n";

}else {	# $model eq "branchsite_test"
	my $pm = new Parallel::ForkManager(2);

	$pm -> run_on_finish( sub {     # reep son-threads' $report
		$_[-1] or die "ERROR: $basename says, codeml failed\n";
                my @recall = @{$_[-1]};
                $recall[0]=~/.*_(branch.*?)\.codemlstudio/;
                $mlcfolder{$1}=$recall[0];
                $lnL{$1}=$recall[1];
        });

	foreach $key (qw/branchsite_null branchsite/){	# launch threads
		$pm -> start and next;
		@report=&run_codeml($key);
		`cp $report[0]/mlc $outfile\_$key.mlc` if $keepmlc;
                $pm -> finish(0,\@report);
	}

	$pm -> wait_all_children;  #recycle threads	
	
	$p = &chi2($lnL{branchsite_null}, $lnL{branchsite});

	print OUT ">$basename\tbranchsite_null\t$lnL{branchsite_null}\tbranchsite\t$lnL{branchsite}\t$p";

	$mlcfolder=$mlcfolder{branchsite};	# as follows is the same to when model is "branchsite"
	$mlc=`cat $mlcfolder/mlc`;
	$mlc=~/Bayes Empirical Bayes.*:\n(.*)\n\s+The grid/s;
	$positive=$1;
	if ($positive=~/\*/s) {
		@_=split /\s+/, `grep "ls =" $mlcfolder/mlc`;
		$aanum=$_[-1];

		$InK="0" x $aanum;	# this $InK string is going to show you positive sites with "."(*) or "!"(**)

		@lines= split /\n/, $positive;
		foreach $_ (@lines){
			next unless /\*/;
			@_=split;
			$mark = /\*\*/? "!" : ".";
			substr($InK, $_[0]-1, 1)=$mark;
			$reportpos .= /\*\*/? "$_[0]!," : "$_[0].,";
		}
		print OUT "\t$reportpos";
		printf OUT "\n%-5s", "#InP";
		print OUT "\t$InK\n";

		foreach $key (keys %pep){
			$key=~/^(...)/;
			$head=$1;
			printf OUT "%-5s", "#$head";
			print OUT "\t$pep{$key}\n";
		}
	}else{
		print OUT "\n# No positive selected site is detected\n";
	}
}

close OUT;

`/bin/rm -rf $basename*codemlstudio`;


######################################## subs #################################

sub run_codeml {
	my $model=shift;

	my $studio = "$basename\_$model.codemlstudio";
	-e $studio or mkdir "$studio";
	chdir "$studio";

	####### make ctl #######
	open CTL, ">>codeml.ctl" or die $!;
	printf CTL "%50s", "seqfile = "; print CTL "$align";
	printf CTL "%50s", "treefile = " and print CTL "$tree" if $tree;
	printf CTL "%50s", "outfile = "; print CTL "mlc\n\n";

	printf CTL "%50s", "noisy = "; print CTL "0\n";
	printf CTL "%50s", "verbose = "; print CTL "0\n";
	printf CTL "%50s", "runmode = "; print CTL $model =~ /pairwise/ ? "-2" : "0";     # 0: user tree;  1: semi-automatic;  2: automatic; 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
	print CTL "\n\n";

	printf CTL "%50s", "seqtype = "; print CTL "1\n";       # 1:codons; 2:AAs; 3:codons-->AAs
	printf CTL "%50s", "CodonFreq = "; print CTL "0\n\n";

	printf CTL "%50s", "clock = "; print CTL "0\n";
	printf CTL "%50s", "aaDist = "; print CTL "0\n";
	printf CTL "%50s", "aaRatefile = "; print CTL "/home/k_d239/usr/local/paml/dat/jones.dat\n\n";      # change it according to your path; only used for aa seqs with model=empirical(_F); dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

	printf CTL "%50s", "model = "; print CTL ($model eq "branch_m0" || $model =~ /pairwise/) ? "0" : $model eq "branch_free" ? "1" : "2";
	print CTL "\n\n";

	printf CTL "%50s", "NSsites = "; print CTL $model =~ /branchsite/? "2" : "0";
	print CTL "\n\n";

	printf CTL "%50s", "icode = "; print CTL "$mt\n";
	printf CTL "%50s", "Mgene = "; print CTL "0\n\n";

	printf CTL "%50s", "fix_kappa = "; print CTL "0\n";
	printf CTL "%50s", "kappa = "; print CTL "2\n";
	printf CTL "%50s", "fix_omega = "; print CTL $model =~ /null/ ? "1" : "0";
	print CTL "\n";
	printf CTL "%50s", "omega = "; print CTL "$fixw\n\n";

	printf CTL "%50s", "fix_alpha = "; print CTL "1\n";
	printf CTL "%50s", "alpha = "; print CTL "0.\n";
	printf CTL "%50s", "Malpha = "; print CTL "0\n";
	printf CTL "%50s", "ncatG = "; print CTL "8\n\n";

	printf CTL "%50s", "getSE = "; print CTL "0\n";
	printf CTL "%50s", "RateAncestor = "; print CTL "0\n\n";

	printf CTL "%50s", "Small_Diff = "; print CTL ".5e-6\n";
	printf CTL "%50s", "cleandata = "; print CTL "0\n";
	printf CTL "%50s", "method = "; print CTL "0\n";

	close CTL;
	####### made ctl ########

	`codeml` or die $!;

	@_ = split /\s+/, `grep lnL mlc` or die $!;
	my $lnL = $model=~/pairwise/? $_[-1] : $_[-2];
	$lnL =~ s/=//;

	chdir "..";
	my @return;
	push @return, $studio;
	push @return, $lnL;
	return @return;
	# return the name of the working folder and lnL value
}

sub chi2 {
	my $d = 2 * abs($_[0]-$_[1]);
	`chi2 1 $d` =~ /prob = (.*?) =/;
	return $1;
}

sub cds2aa {
        undef $out;
        my $in = uc(shift @_);
        my(%genetic_code) = (
        'TCA' => 'S',    # Serine
        'TCC' => 'S',    # Serine
        'TCG' => 'S',    # Serine
        'TCT' => 'S',    # Serine
        'TTC' => 'F',    # Phenylalanine
        'TTT' => 'F',    # Phenylalanine
        'TTA' => 'L',    # Leucine
        'TTG' => 'L',    # Leucine
        'TAC' => 'Y',    # Tyrosine
        'TAT' => 'Y',    # Tyrosine
        'TAA' => '*',    # Stop
        'TAG' => '*',    # Stop
        'TGC' => 'C',    # Cysteine
        'TGT' => 'C',    # Cysteine
        'TGA' => $mt? 'W' : '*',    # Stop
        'TGG' => 'W',    # Tryptophan
        'CTA' => 'L',    # Leucine
        'CTC' => 'L',    # Leucine
        'CTG' => 'L',    # Leucine
        'CTT' => 'L',    # Leucine
        'CCA' => 'P',    # Proline
        'CCC' => 'P',    # Proline
        'CCG' => 'P',    # Proline
        'CCT' => 'P',    # Proline
        'CAC' => 'H',    # Histidine
        'CAT' => 'H',    # Histidine
        'CAA' => 'Q',    # Glutamine
        'CAG' => 'Q',    # Glutamine
        'CGA' => 'R',    # Arginine
        'CGC' => 'R',    # Arginine
        'CGG' => 'R',    # Arginine
        'CGT' => 'R',    # Arginine
        'ATA' => $mt? 'M' : 'I',    # Isoleucine
        'ATC' => 'I',    # Isoleucine
        'ATT' => 'I',    # Isoleucine
        'ATG' => 'M',    # Methionine
        'ACA' => 'T',    # Threonine
        'ACC' => 'T',    # Threonine
        'ACG' => 'T',    # Threonine
        'ACT' => 'T',    # Threonine
        'AAC' => 'N',    # Asparagine
        'AAT' => 'N',    # Asparagine
        'AAA' => 'K',    # Lysine
        'AAG' => 'K',    # Lysine
        'AGC' => 'S',    # Serine
        'AGT' => 'S',    # Serine
        'AGA' => $mt? '*' : 'R',    # Arginine
        'AGG' => $mt? '*' : 'R',    # Arginine
        'GTA' => 'V',    # Valine
        'GTC' => 'V',    # Valine
        'GTG' => 'V',    # Valine
        'GTT' => 'V',    # Valine
        'GCA' => 'A',    # Alanine
        'GCC' => 'A',    # Alanine
        'GCG' => 'A',    # Alanine
        'GCT' => 'A',    # Alanine
        'GAC' => 'D',    # Aspartic Acid
        'GAT' => 'D',    # Aspartic Acid
        'GAA' => 'E',    # Glutamic Acid
        'GAG' => 'E',    # Glutamic Acid
        'GGA' => 'G',    # Glycine
        'GGC' => 'G',    # Glycine
        'GGG' => 'G',    # Glycine
        'GGT' => 'G',    # Glycine
        );
        for $codon ($in =~ /(...)/g) {
                $out .= exists $genetic_code{$codon}? $genetic_code{$codon} : 'x';
        }
        return $out;
}
