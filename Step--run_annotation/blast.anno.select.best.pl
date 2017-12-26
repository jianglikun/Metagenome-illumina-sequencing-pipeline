#!/usr/bin/perl -w
use strict;
my($blast,$info,$out)=@ARGV;
open IN1,$blast||die;
open IN2,$info||die;
open OUT,">$out"||die;

my (%info,%hash);
while(<IN2>){
    chomp;
    my @arr = split(/\t/);
    $info{$arr[0]}=$arr[1];
}
print OUT "Gene_id\tDescription\n";
while(<IN1>){
    chomp;
    my @arr2 = split(/\t/);
    if(!exists $hash{$arr2[0]}){
	print OUT "$arr2[0]\t$info{$arr2[1]}\n";
	$hash{$arr2[0]}=1;
    }
}
close IN1;close IN2;
close OUT;
