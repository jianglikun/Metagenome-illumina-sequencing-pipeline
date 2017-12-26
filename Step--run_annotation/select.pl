my $in =shift;
my $out=shift;
open IN,$in||die;
open OUT,">$out"||die;

print OUT "Gene_id\tARDB_superclass\tARDB_drug\n";
while(<IN>){
    chomp;
    my @arr=split(/\t/) unless /Resistance Type/ ;
    my @id=split(/, /,$arr[-1]);
    for (my $i=0;$i < @id;$i++){
        print OUT "$id[$i]\t$arr[0]\t$arr[2]\n";
    }
}
close IN;close OUT;
