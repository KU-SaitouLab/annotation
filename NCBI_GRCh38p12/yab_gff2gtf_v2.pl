#!/usr/bin/perl -s
use 5.8.0;
use strict;
use warnings;
our ($f   );
our ($o   );
my  ($fh  );
my  ($fo  );
my  (%gid );
my  (%bty );
my  (%tid );
my  ($id  );
my  (%gflg);
my  (%tflg);
my  ($geneid );
my  ($transid );
my  ($gid );
my  ($tid );
my  ($pid );
my  ($gnm );
my  ($bty );

open $fh,$f or die "$!:$f\n";
open $fo,">$o" or die "$!:$o\n";
while (my $l=<$fh>) {
	chomp($l);
	my @l=split/\t/,$l;
	($id)      = ($l[8]=~/^ID=([^;]+)/);
	($geneid)  = ($l[8]=~/GeneID:([^;,]+)/);
	($transid) = ($l[8]=~/transcript_id=([^;]+)/);
	($gnm)     = ($l[8]=~/;gene=([^;]+)/);
	($bty)     = ($l[8]=~/gene_biotype=([^;]+)/);
	if ($l[8] !~ /Parent=/) {
		$gid=$id;
		$gflg{$gid}++;
		$bty{$gid}=$bty;
		print $fo join("\t",@l[0..7],join("; ","gene_id \"$gid\"","gene_name \"$gnm\"","gene_biotype \"$bty\"")) . ";\n";
	} else {
		($pid)     = ($l[8]=~/Parent=([^;]+)/);
		if ($gflg{$pid}) { # transcript/RNA
			$gid=$pid;
			$tid=$id;
			$tid{$tid}=$tid;
			$tflg{$tid}++;
			$gid{$id}=$pid;
		} else {
			$gid=$gid{$pid};
			$gid{$id}=$gid{$pid};
			$tid=$tid{$pid};
			$tid{$id}=$tid{$pid};
		}
		print $fo join("\t",@l[0..7],join("; ","gene_id \"$gid\"","transcript_id \"$tid\"","gene_name \"$gnm\"","gene_biotype \"$bty{$gid}\"")) . ";\n";
	}
}
close $fh;
close $fo;
