#!/usr/bin/perl
###############################################################
# Program Name: seven_filters.pl             Version: 1.10
#
# Description: A fast and easy way to filter RNA editing 
#  sites from our pipeline.
#
#
# Author: Shan Sabri          Date: 2016-09-01
#
#
# Revision History
#
# Version                 Date                  Who
#-------------------------------------------------------------
#  1.0                 2016-08-22             Shan Sabri  
###############################################################


use strict;
use warnings;

use lib '/home/ssabri/perl5/lib/perl5';
use lib '/u/home/s/ssabri/project-gxxiao/perl_modules/Statistics-R-0.33/lib';
use lib '/u/home/s/ssabri/project-gxxiao/perl_modules/Statistics-R-0.33/inc';

use Statistics::R;

$ENV{'PATH'} = "/home/gxxiao/apps/bin/";
my $paraf = "config_pre";
my ($para,$group) = parseParas($paraf);
my $p_ss = $$para{"P_SS"};
my $d_editout = $$para{"D_log_editing"};
my $d_output = $$para{"D_out_editing"};
my $d_work = $$para{"D_work"};


my @par = split(/\^/,$ARGV[0]);
my $gr = $par[0];
my $CHR = $par[1];
my $TYPE = $par[2];
my $gg = $par[3];

my @GR = split(/\./,$gr);

my $f_raw = "$d_editout/$GR[0].$CHR.$TYPE.readpos.0.dist2end.1.txt";
print $f_raw,"\n";
my $f_out = "$d_work/out/$GR[0]/editing.$CHR.$TYPE";

my @GG = split(/\./,$gg);

my %posRaw = ();
my %pos = ();

get_pos($f_raw);

filter_basic();


if (0){
	if ($$para{"Filter_snp"} && %pos){
		%posRaw = %pos;
		%pos = ();
		filter_snp();
		print "Filter_snp\n";
	}
}


if ($$para{"Filter_degree"} && %pos){
	%posRaw = %pos;
	%pos = ();
	filter_degree();
	print "Filter_degree\n";
}


if ($$para{"Filter_editing"} && %pos){
	%posRaw = %pos;
	%pos = ();
	filter_editing();
	print "Filter_editing\n";
}


if ($$para{"Filter_intron"} && %pos){
	%posRaw = %pos;
	%pos = ();
	filter_intron();
	print "Filter_intron\n";
}


if ($$para{"Filter_homopoly"} && %pos){
	%posRaw = %pos;
	%pos = ();
	filter_homopoly();
	print "Filter_homopoly\n";
}


if ($$para{"Filter_sr"} && %pos){
	%posRaw = %pos;
	%pos = ();
	filter_sr();
	print "Filter_sr\n";
}


if ($$para{"Filter_strand"} && %pos){
	%posRaw = %pos;
	%pos = ();
	filter_strand();
	print "Filter_strand\n";
}


open OUT, ">$f_out" || die "Can't open the file!";
foreach my $p (sort keys %pos){
	print OUT $p,"\n";
}
close OUT;
print "$GR[0]/editing.$CHR.$TYPE finished\n";


sub get_pos{
	my ($f) = @_;
	open IN, $f || die "Can't open the file!";

	while(<IN>) {
		chomp;
		$posRaw{$_}++;
	}
	close IN;
}


sub filter_basic{
	
	#my $outf = "res_basic_$gr";
	my $p_count_sum = $$para{"P_basic_count_sum"};
	my $p_count_editing = $$para{"P_basic_count_editing"};
	my $p_log_likelihood = $$para{"P_basic_log_likelihood"};

	#open OUT, ">>$outf" || die "Can't open the file!";

	foreach my $p (sort keys %posRaw){

		my @line=split(/\t/,$p);

		#print $line[4],"\t",$line[5],"\n";
		if ($p_ss == 0){
			next if ($line[4] < $p_count_editing);
			next if ($line[5] < $p_count_sum);
			next if ($line[$#line] < $p_log_likelihood);
			#print OUT "$_\n";
			$pos{$p}++;
		}else{
			next if ($line[3] < $p_count_editing);
			next if ($line[4] < $p_count_sum);
			next if ($line[$#line-1] < $p_log_likelihood);
			#print OUT "$_\n";
			$pos{$p}++;
		}
	}
	#close OUT;
	return(1);
}


sub filter_snp{

	my $f_snp = $$para{"F_snp"};
	#my $f_snp_str = $$para{"F_snp_str"};

	my %snp;
	open IN, $f_snp || die "Can't open the file!";
	while(<IN>) {
		chomp;
		my @line = split /\s+/;
		$snp{"$line[0]:$line[1]"} = $line[$#line];# last columan is the strand
	}
	close IN;
	

	foreach my $p (sort keys %posRaw){

		my @line=split(/\t/,$p);
		my $id="$line[0]:$line[1]";
		next if (exists($snp{$id}));
		#next if (exists($snp{$id}) && ($line[$#line] eq $snp{$id}));
		$pos{$p}++;
	}
	return(1);
}

sub filter_degree{

	my $p_degree = $$para{"P_degree"};

	foreach my $p (sort keys %posRaw){

		my @line=split(/\t/,$p);
		next if (($p_ss == 0) && ($line[6] >= $p_degree));
		next if (($p_ss != 0) && ($line[5] >= $p_degree));
		$pos{$p}++;
	}
	return(1);
}


sub filter_editing{

	my $p_count_editing = $$para{"P_count_editing"};
	my $p_ratio_editing = $$para{"P_ratio_editing"};

	foreach my $p (sort keys %posRaw){

		my @line=split(/\t/,$p);
		#print $p,"\t",$line[3],"\t",$p_count_editing,"\n";
		if ($p_ss == 0){
			next if (($line[4] < $p_count_editing) || ($line[6] < $p_ratio_editing));
		}else{
			next if (($line[3] < $p_count_editing) || ($line[5] < $p_ratio_editing));
		}
		$pos{$p}++;
	}
	return(1);
}


sub filter_intron{

	my $p_intron_site = $$para{"P_intron_site"};
	my $f_ano = $$para{"F_gene_ano"};

	my %spl = ();
	open IN, $f_ano ||  die "Can't open the file!\n";
	while(<IN>) {
		next if(/^#/);
		chomp;
		my @line=split (/\t/, $_);

		my @es=split (/\,/, $line[9]);
		my @ee=split (/\,/, $line[10]);

		for (my $i=1; $i<=$#es; $i++) {

			my $is = $ee[$i-1]+1; #intron start
			my $is_ex = $is+$p_intron_site-1; # intron start extension
			#my $reg = "$is-$is_ex";
			$spl{$line[2]}{"$is-$is_ex"}++;

			my $ie = $es[$i]; # intron end
			my $ie_ex = $ie-$p_intron_site+1; # intron end extension
			$spl{$line[2]}{"$ie_ex-$ie"}++;
		}
	}
	close IN;
	if(0){
	print scalar keys %{$spl{"chr1"}},"\t1\n";

	$f_ano = "/home/shansabri/runSixFilters/dependencies/ucsc.rgd.gene.rn4.060712";
	open IN, $f_ano ||  die "Can't open the file!\n";
	while(<IN>) {
		next if(/^#/);
		chomp;
		my @line=split (/\t/, $_);
		my @es=split (/\,/, $line[9]);
		my @ee=split (/\,/, $line[10]);
		for (my $i=1; $i<=$#es; $i++) {
			my $is = $ee[$i-1]+1; #intron start
			my $is_ex = $is+$p_intron_site-1; # intron start extension
			$spl{$line[2]}{"$is-$is_ex"}++;
			my $ie = $es[$i]; # intron end
			my $ie_ex = $ie-$p_intron_site+1; # intron end extension
			$spl{$line[2]}{"$ie_ex-$ie"}++;
		}
	}
	close IN;
	print scalar keys %{$spl{"chr1"}},"\t2\n";
	$f_ano = "/home/shansabri/runSixFilters/dependencies/ucsc.rgd.gene.rn4.060712";
	open IN, $f_ano ||  die "Can't open the file!\n";
	while(<IN>) {
		next if(/^#/);
		chomp;
		my @line=split (/\t/, $_);

		my @es=split (/\,/, $line[8]);
		my @ee=split (/\,/, $line[9]);

		for (my $i=1; $i<=$#es; $i++) {

			my $is = $ee[$i-1]+1; #intron start
			my $is_ex = $is+$p_intron_site-1; # intron start extension
			$spl{$line[1]}{"$is-$is_ex"}++;
			#print "$is-$is_ex","\n";
			my $ie = $es[$i]; # intron end
			my $ie_ex = $ie-$p_intron_site+1; # intron end extension
			$spl{$line[1]}{"$ie_ex-$ie"}++;
		}
	}
	close IN;
	print scalar keys %{$spl{"chr1"}},"\t3\n";
	}

	foreach my $p (sort keys %posRaw){

		my @line=split(/\t/,$p);

		my $flag = 0;
		foreach my $pos (sort keys %{$spl{$line[0]}}){
			my @loc=split (/-/, $pos);
			#if ($loc[0] <= $line[1]  and $line[1] <= $loc[1]) {
			if (($line[1] >= $loc[0]) && ($line[1] <= $loc[1])){
				$flag = 1;
				last;
			}
		}
		next if ($flag);
		$pos{$p}++;
	}
	return(1);
}

sub filter_homopoly{

	my $f_chr_seq = $$para{"F_chr_seq"};
	my $p_homo_length = $$para{"P_homo_length"};

	my $Cseq = ();
	my $flag = 0;
	open IN, $f_chr_seq || "Can't open the file!";
	while(<IN>) {
		chomp;
		if (/^>\Q$CHR\E$/) {
			#print $_,"\n";
			while(<IN>) {
				chomp;
				if (/^>/){
					$flag = 1;
					last;
				}
				$Cseq.=uc ($_);
			}
		}
		last if ($flag);
	}
	close IN;
	#print "length:",length($Cseq),"\t",$CHR,"-----\n";

	my @ty=split (/\.to\./, $TYPE);

	foreach my $p (sort keys %posRaw){
		my @line=split (/\t/, $p);
		my $chr=$line[0];
		my $pos=$line[1];
		my $start=$pos-4;
		my $estr=$line[$#line];
		$estr=1 if ($p_ss == 0);
		#my $estr=1;
		my $before=substr($Cseq, $pos-$p_homo_length, $p_homo_length-1);
		my $after=substr($Cseq, $pos, $p_homo_length-1);
		#print $pos-$p_homo_length,"\t",$p_homo_length-1,"\t",$pos,"\t",$p_homo_length-1,"\n";

		my $lseq1="";
		my $lseq2="";

		if ($estr == 1) {
			$lseq1=$before.$ty[0].$after;
			$lseq2=$before.$ty[1].$after;
		} elsif ($estr == -1 ){
			$lseq1=revs($before).($ty[0]).revs($after);
			$lseq2=revs($before).($ty[1]).revs($after);
		}

		my $hflag = 1;
		for (my $i=0; $i<=4; $i++) {
			my $seq=substr($lseq1, $i, $p_homo_length);
			my @seq=split (//, $seq);
			my %seen;
			my @unique=grep {!$seen{$_}++} @seq;
			if ((scalar @unique) == 1) {
				#$hflag="$chr\t$pos\t".($start+$i)."\t".($start+$i+4)."\t$seq";
				$hflag = 0;
				last;
			}
		}
		next if (!$hflag);

		for (my $i=0; $i<=4; $i++) {
			my $seq=substr($lseq2, $i, 5);
			my @seq=split (//, $seq);
			my %seen;
			my @unique=grep {!$seen{$_}++} @seq; 
			if ((scalar @unique) == 1) {
				#$hflag="$chr\t$pos\t".($start+$i)."\t".($start+$i+4)."\t$seq";
				$hflag = 0;
				last;

			}
		}
		next if (!$hflag);

		$pos{$p}++;
	}
	return(1);
}

sub revs {
	
	my $seq=$_[0];
	my $rev=$seq;
	$rev=~ tr/ACGTacgt/TGCAtgca/;
	return $rev;
}

sub filter_sr{

	my $f_sr_region = $$para{"Filter_sr_F_region"};

	my %repeat = ();
	open IN, $f_sr_region || die "Can't open the file!";
	while (<IN>) {
		chomp;
		my @line=split /\t/;
		my @pos=split (/\:/, $line[0]);
		if (($pos[0] eq $CHR) && (($line[$#line] eq "Simple_repeat") || ($line[$#line-1] eq "Simple_repeat"))) {
			#$repeat{$pos[1]}="$line[$#line-1];$line[$#line]";
			$repeat{$pos[1]}++;
		}
	}


	foreach my $p (sort keys %posRaw){

		my @line=split(/\t/,$p);
		my $flag = 0;
		foreach my $pos (sort keys %repeat){
			my @loc=split (/\-/, $pos);
			if (($loc[0] <= $line[1])  && ($line[1] <= $loc[1])){
				$flag = 1;
				last;
			}
		}
		next if ($flag);
		$pos{$p}++;
	}
	return(1);
}

sub filter_strand{

	my $p_strand_pvalue = $$para{"P_strand_pvalue"};
	my $d_out_editing = $$para{"D_out_editing"};

	my %posE = ();
	foreach my $p (sort keys %posRaw){

		my @line = split(/\t/,$p);
		$posE{$line[1]} = $p;

	}

	my %str = ();

	#foreach my $g (sort keys %{$group}){
	foreach my $g (@GG){

		my $f ="$d_out_editing/$g.$CHR.sam.all.mismatch.pos.read.info.txt.sorted.no.duplicates";
		#my $f ="$d_out_editing.$gr/$g.$CHR.sam.all.mismatch.pos.read.info.txt.sorted.no.duplicates";
		#print $f,"\n"; ###
		open IN , $f  || die "Can't open the file!";
		while(<IN>) {
			chomp;
			my @line=split (/\t/, $_);
			my $pos=$line[2];
			my $strand=$line[3];
			if (exists $posE{$pos}) { 
				if ($line[5] eq $line[6]) {
					if ($strand eq "1") {
						$str{$pos}{"P"}++;
					} elsif ($strand eq "-1") {
						#$str{$pos}.="N";
						$str{$pos}{"N"}++;
					}
				} elsif ($line[5] ne $line[6]) {
					if ($strand eq "1") {
						#$str{$pos}.="F";
						$str{$pos}{"F"}++;
					} elsif ($strand eq "-1") {
						#$str{$pos}.="R";
						$str{$pos}{"R"}++;
					}
				}
			}
		}
		close IN;
	}

	my @ps = ();
	my @rplus =();
	my @rminus =();
	my @mplus =();
	my @mminus =();
	foreach my $p (sort {$a <=> $b} keys %str){
		push(@ps,$p);
		if (exists($str{$p}{"P"})){
			push(@rplus,$str{$p}{"P"});
		}else{
			push(@rplus,0);
		}
		if (exists($str{$p}{"N"})){
			push(@rminus,$str{$p}{"N"});
		}else{
			push(@rminus,0);
		}
		if (exists($str{$p}{"F"})){
			push(@mplus,$str{$p}{"F"});
		}else{
			push(@mplus,0);
		}
		if (exists($str{$p}{"R"})){
			push(@mminus,$str{$p}{"R"});
		}else{
			push(@mminus,0);
		}
	}
	
	return(1) if ($#rplus == -1);
	#for(my $i=0;$i<= $#ps;$i++){
		#print $ps[$i],"\t",$rplus[$i],"\t",$rminus[$i],"\t",$mplus[$i],"\t",$mminus[$i],"\n";
		#print $posE{$ps[$i]},"\n";
	#}
	#print $#rplus,"\t-rp-\n";
	#print $#rminus,"\t-rm-\n";
	#print $#mplus,"\t-mp-\n";
	#print $#mminus,"\t-mm-\n";
	#print $#ps,"\t-ps-\n";

	my $R = Statistics::R->new();
	#my $R = Statistics::R->new("log_dir"=>$log_dir);

	$R->set('rp',\@rplus);
	$R->set('rm',\@rminus);
	$R->set('mp',\@mplus);
	$R->set('mm',\@mminus);
	#$R->run(qq 'cat(rp[1])');
	$R->run(qq 'dat = cbind(rp,rm,mp,mm)');
	$R->run(qq 'p = apply(dat,1,function(x){fisher.test(matrix(x,2,2))[["p.value"]]})'); 
	#$R->send(qq 'p = p.adjust(p)');
	my $pv = $R->get('p'); 
	#$R->stopR();
	#print $pv,"-0-\n";

	if ($#rplus == 0){
		#print $pv,"-1-\n";
		if ($pv >= $p_strand_pvalue){
			$pos{$posE{$ps[0]}}++;
		}
	}else{ 
	        #print "else-2-\n";
	        #print $pv,"-2-\n";
		for(my $i=0;$i<= $#ps;$i++){
			#print $ps[$i],"\t",$rplus[$i],"\t",$rminus[$i],"\t",$mplus[$i],"\t",$mminus[$i],"\n";
			next if ($$pv[$i] < $p_strand_pvalue);
			$pos{$posE{$ps[$i]}}++;
		}
	}
	return (1);
}

sub parseParas{
	my ($conf) = @_;
	my %paras = ();
	my %group = ();
	open IN, $conf || die "Can't open the file!";
	while (<IN>){ 
		if (/^[GS]P/){ 
			chomp; 
			my @line = split /\t/; 
			if ($line[1] eq "D_group"){
				for(my $i=2;$i<$#line;$i++){
					$group{$line[$i]}++;
				}
			}
			$paras{$line[1]} = $line[2]; 
		} 
	} 
	close IN; 
	return(\%paras,\%group); 
}