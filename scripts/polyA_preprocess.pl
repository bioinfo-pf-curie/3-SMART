#!/usr/bin/perl
use warnings;
use strict;

my $dir = "/bioinfo/guests/vbernard/Bureau/HiSeq/PolyA_201311/svn/analyse/report/"; ## fichiers initiaux
my $OUT = "/bioinfo/guests/vbernard/Bureau/HiSeq/PolyA_201311/svn/analyse/info/"; ## fichiers finaux filtrés
opendir(DIR, $dir) or die $!;
while (my $D = readdir(DIR)) {
    next if ($D !~ m/^A162T/);
    opendir(DIR2, "$dir$D") or die $!;
    system("mkdir $OUT$D");
    my $ligne = 0;
    while (my $f = readdir(DIR2)) {
	next if ($f !~ m/(.+).fastq/);
	print "$dir$D\/$f en cours...\n";
	open(IN,"$dir$D\/$f");
	open(OUT,">$OUT$D\/$f");
	open(OUT2,">$OUT$D\/$f\_polyAremoved.txt");
	print "$OUT$D\/$f\_polyAremoved.txt\n";
	my $NB = 0;
	my $NBg=0;
	my $leg;
	while(<IN>){
	    $ligne++;
	    if($ligne == 4){$ligne=0;}
	    my $L = $_;
	    chomp($L);
	    if($ligne == 2){
		$NB=0;
		$NBg=0;
		my $seq;
		#####poly A
		if($L =~/^(.*[CGT])(AAA+)$/){ ## 3A ou plus
		    $seq = $1;
		    $NB = length($2);
		}elsif($L =~/^(A+)$/){ ## que des A
		    $seq="";
		    $NB = length($1);
		}else{ $seq=$L;}
		###### G tail
		my $test;
		if($seq =~/^(GG+)([ACT].+)$/){ ## 2G ou plus
		    print OUT "$2\n";
		    $test = $2;
		    $NBg = length($1);
		}elsif($seq =~/^(G+)$/){ ## que des G (possible une fois qu'on a retiré les polyA)
		    print OUT "\n";
		    $test = "";
		    $NBg = length($1);
		}else{
		    $test = $seq;
		    print OUT "$seq\n";
		}
		print OUT2 "$leg\t$NB\t$NBg\n";
	    }else{
		if ($ligne != 0){
		    print OUT "$L\n";
		    if($L=~ /^(\@HWI.+)$/){$leg =$1;}
		}else{
		    my $tmp = substr $L, 0, 50-$NB+1;
		    my $tmp2="";
		    if($NB<47){$tmp2 = substr $tmp, $NBg, length($tmp)-$NBg;}
		    print OUT "$tmp2\n";
		}
	    }
	}
	close(IN);
	close(OUT);
	close(OUT2);
    }
    closedir(DIR2);
}
closedir(DIR);
