# Version 1.0, July 2 2009 
# Perl program for antibiotic resistance genes annotation
# This program requires the installation of BioPerl toolkit
# and BLAST programs
#
# Maintained by Bo Liu <boliu-at-cs.umd.edu>
#
# Copyright Bo Liu
#
# You may distribute this program under the same terms as perl itself
#
# History
# July 2, 2009 ardbAnno.pl written by Bo Liu

#!usr/bin/perl
use strict;
use warnings;
use ardbAnno;
use Bio::SeqIO;

my $evalueCutoff = 1e-10; #change this variable to try different evalue cutoffs


my $program = ""; #the blast program will be automatically chosen
my %genome = getGenome(); #get a list of genomes to be annotated
print "######################################\n";
foreach my $genome (keys %genome)
{
    my $seqfile = $genome{$genome};
    unless ($program)
    {
	my $seqi = new Bio::SeqIO (-file => $seqfile, -format => "fasta");
	my $seqobj = $seqi -> next_seq();
	my $alph = $seqobj -> alphabet;
	if ($alph eq "dna")
	{
	    $program = "blastx";
	}
	elsif ($alph eq "protein")
	{
	    $program = "blastp";
	}
	else
	{
	    print "The input sequences are not correct.\n";
	    exit;
	}
	    
    }
    my $outfile = $genome . "\.$program";
    my $db = "./blastdb/resisGenes.pfasta";
    print "BLASTing $seqfile against ARDB\n";
    if (-e $outfile)
    {
	next;
    }
    system("blastall -p $program -d $db -i $seqfile -e $evalueCutoff -b 3 -v 3 -o $outfile");
}
print "######################################\n\n";

#extract annotation from blast results
print "######################################\n";
foreach my $genome (keys %genome)
{
    my $blastfile = $genome . "\.$program";
    print "Annotating BLAST file $blastfile\n";
    if (-e "$genome\.anno")
    {
	next;
    }
    annoBlast($blastfile);
}
print "######################################\n\n";

#merge annotation files into excel sheet
print "Create Excel Table output.xls\n\n";
mergeAnno(keys %genome);

exit;
