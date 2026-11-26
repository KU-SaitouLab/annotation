#!/usr/bin/perl -s
use 5.8.0;
use strict;
use warnings;
our ($f      );
our ($o      );
our ($chr    );
our ($type   );
our ($verbose);
my ($fh      );
my (@chr     );
my (@data_all);
my (%pos     );
my (%ext_ps  );
my (%ext_pt  );
my (@chr_all );
my (%pseudo  );
my (%partial );
my (%rna_id  );
$verbose=0 unless defined($verbose);
#$chr="1..22,X,Y,M" unless $chr;
$o = "-" unless $o;
@chr=split/,/,$chr;
if ($chr[0]=~/(\d+)\.\.(\d+)/) {$_=shift(@chr);for (my $i=$2;$i>=$1;$i--) {unshift(@chr,$i);}}
$_="chr$_" for @chr;
print "Reading file...\n" if $verbose;
if ($f=~/.gz$/) {
	open($fh, "gzip -dc $f|") or die "$!:$f";
} else {
	open($fh, "$f") or die "$!:$f";
}
	@data_all=<$fh>;
close $fh;
my @data_real=grep{$_!~/^[#\r\n]/ and $_} @data_all;
my @data_gene=grep{(split/\t/,$_)[2] =~/gene$/} @data_real;
if ($type=~/pseudo/i) {
	foreach (@data_gene) {
		my $attr=(split/\t/,$_)[8];
		my ($id1) = ($attr=~/ID=([^;]+)/i);
		$pseudo{$id1}++ if $attr=~/pseudo=true/i;
	}
}
if ($type=~/partial/i) {
	foreach (@data_gene) {
		my $attr=(split/\t/,$_)[8];
		my ($id1) = ($attr=~/ID=([^;]+)/i);
		$partial{$id1}++ if $attr=~/partial=true/i;
	}
}
foreach (@data_real) {
	my $attr=(split/\t/,$_)[8];
	my ($id1) = ($attr=~/ID=([^;]+)/i);
	my ($id2) = ($attr=~/Parent=([^;]+)/i);
	if ($id2) {
		$pos{$id1} = $id2;	
	} else {
		$pos{$id1} = $id1;
	}
	push(@chr_all,(split/\t/,$_)[0]);
	($rna_id{$id1}) = ($attr=~/Genbank:([^;]+)/) if $id1=~/rna/;
}
unless ($chr) {
	my %h;
	@chr=grep { ! $h{ $_ }++ } @chr_all;
}
# trenascript, exon などのparent_idにgeneのidを設定する。
my %name=();
my %type=();
my %gid=();
my %gname=();
my %ggene=();
foreach (@data_real) {
	my $attr=(split/\t/,$_)[8];
	my ($id1) = $attr=~/ID=([^;]+)/i;
	my ($name) = $attr=~/Name=([^;]+)/i;
	print "$id1 , $pos{$id1} , $pos{$pos{$id1}} \n" if $verbose and !defined($pos{$pos{$id1}});
	if ($id1 ne $pos{$id1}) {
		$pos{$id1} = $pos{$pos{$id1}} while $pos{$id1} ne $pos{$pos{$id1}};
	}
	push(@{$name{$pos{$id1}}},$name);
	push(@{$type{$pos{$id1}}},(split/\t/,$_)[2]);
	($gid{$id1}) = $attr=~/GeneID:(\d+)/ if (split/\t/,$_)[2] eq "gene";
	($gname{$id1}) = $attr=~/Name=([^;]+)/i if (split/\t/,$_)[2] eq "gene";
	($ggene{$id1}) = $attr=~/gene=([^;]+)/i if (split/\t/,$_)[2] eq "gene";
}
print "Reading file...done.\n" if $verbose;
if ($o=~/.gz$/) {
	open $fh,"|gzip > $o" or die "$!:$o";
} else {
	open $fh,">$o" or die "$!:$o";
}
foreach my $chr (@chr) {
	print "Analyzing chromosome $chr.\n" if $verbose;
	my @data_chr;
	if ($type=~/alldata/i) {
		@data_chr=grep{$_=~/^$chr\t/} @data_real;
	} else {
		@data_chr=grep{$_=~/^$chr\t/} @data_gene;
	}
	chomp(@data_chr);
	foreach my $gt (@data_chr) {
		print ("Parsing $gt\n") if $verbose > 1;
		my (@data_sub_ss,@data_sub_sa,@data_sub_aa,@data_sub_as,@dist)=();
		my ($st,$end,$end_dist);
		@_=split/\t/,$gt;
		my $attr=(split/\t/,$gt)[8];
		my ($id1) = ($attr=~/ID=([^;]+)/i);
		my ($pid1) = ($attr=~/Parent=([^;]+)/i);
		$pid1="-" unless $pid1;
		my %h=();
		my @names=@{$name{$pos{$id1}}};
		@names=grep{$_}@names;
		@names=grep{!$h{$_}++}@names;
		%h=();
		my @types=@{$type{$pos{$id1}}};
		@types=grep{$_}@types;
		@types=grep{!$h{$_}++}@types;
		$gt.=";Name=" . $gname{$id1}         if ($type=~/ifname/i and $gt !~ /Name=/);
		$gt.=";Name=" . $gname{$id1}         if ($type=~/ifgene/i and $gt !~ /gene=/);
		$gt.=";Top_id=" . $pos{$id1}         if $type=~/topid/i;
		$gt.=";Name_all=" . join(",",@names) if $type=~/name/i;
		$gt.=";pseudo=true"                  if ($type=~/pseudo/i and $pseudo{$pos{$id1}});
		$gt.=";partial=true"                  if ($type=~/partial/i and $partial{$pos{$id1}});
		$gt.=";Type_all=" . join(",",@types) if $type=~/type/i;
		$gt.=";GeneID=" . $gid{$pos{$id1}}   if $type=~/geneid/i;
		$gt.=";Genbank:" . $rna_id{$pid1}    if ($type=~/rna_id/i and $pid1=~/^rna/);
		print $fh "$gt\n";
	}
}
close $fh;