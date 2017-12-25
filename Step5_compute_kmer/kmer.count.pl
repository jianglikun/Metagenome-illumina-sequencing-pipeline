#!/usr/bin/perl -w
use strict;
use File::Basename;
my $usage="Usage: $0 reads1 reads2 ...\n\n";
unless(@ARGV){
    die $usage;
}
my @rsem_files=@ARGV;
unless(scalar @rsem_files >= 2){
    die $usage;
}
my $max      = 0;
my $length   = 0;
my $sum      = 0;

main:{
    foreach my $file(@rsem_files){
	    open(my $fh,"gunzip -c $file|") or die "Error:cannot open file $file\n";
		my $counter = 1;
		while(<$fh>){
		    chomp;
			my $line = $_;
			if ( $counter == 1 ) {
			    $line =~ m/^@(.*)\/[12]/;
			}
			if ( $counter == 2 ) {
			    if ( $max < length($line) ) {
				    $max = length($line);
				}
				$length+=length($line);
				$sum++;
			}
			$counter++;
			if ( $counter == 5 ) {
				$counter = 1;
			}
		}
	}
	my $avg = 0;
	unless ( $sum == 0 ) {
	    $avg = int( ( $length / $sum ) + 0.5 );
	}
	my $kmer;
	if ( $avg % 2 == 0 ) {
	    $kmer = $avg / 2;
		if ( $kmer % 2 == 0 ) {
		    $kmer = $kmer + 1;
		}else{
		    $kmer = $kmer + 2;
		}
	}else{
	    $kmer = ( $avg + 1 ) / 2;
		if ( $kmer % 2 == 1 ) {
		    $kmer = $kmer + 0;
		}else{
		    $kmer = $kmer + 1;
		}
	}
	print "$kmer\n";
	exit(0);
}
