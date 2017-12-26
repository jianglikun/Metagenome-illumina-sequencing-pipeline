#!/usr/bin/env perl

# This code derived from the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd 'abs_path';
use MIME::Base64 'decode_base64';
use Pod::Usage 'pod2usage';

#use Carp 'verbose';
#$SIG{__DIE__} = sub { Carp::confess( @_ ) };


my @samples = ();
my $sample_file;
my $temp_dir='.';
chomp(my $basedir = `pwd`);
my $cwd = $basedir;
my $INPUT_FORMAT = 'auto';
my $processors = 7;
my $calculateTaxonomy_previous_calc_coverage_stats_file = 0;
my $version = '1.1.1';
my $length_cutoff = 30;
my $identity_cutoff = 97;
my $quality_cutoff = 20;
my $prevalence_minimum = 2;
my $verbose = 0;
my $OUTPUT_FOLDER = 'RESULTS';

my $src_dir = '/mOTU/scripts/src';
my $bin_dir = '/mOTU/scripts/bin';
my $data_dir;

my %conf = ();
$conf{MOCAT_data_type} = 'solexaqa';
$conf{MOCAT_mapping_mode} = 'allbest';
$conf{screen_soap_cmd} = '-M 4';
$conf{screen_soap_seed_length} = 30;
$conf{screen_soap_max_mm} = 10;
$conf{filter_psort_buffer} = "2G";
$conf{MOCAT_rownames_dir} = './rownames/';

my $IS_PAIRED_END = -1;


my $systemType = `uname -s`;
chomp $systemType;

my $zcat = 'zcat';
if ( $systemType =~ m/Darwin/ ) { $zcat = 'gzcat'; }

# CORE
sub expand_data {
    my $basedir = shift;
    print("$0 will expand needed dependencies into ./motus_data....\n");
    open EXPAND, "|tar xzf - --directory $basedir"
        or die("Cannot open tar process");
    while (<DATA>) {
        print EXPAND decode_base64($_);
    }
    close(EXPAND);
}


sub mkdir_or_die {
    my $dir = shift;
    (-d $dir)
        or make_path($dir)
        or die("Could not create $dir: $!");
}

sub system_ {
    my $cmd = shift;
    if ($verbose) {
        print "Will execute \n\n$cmd\n\n";
    }
    (system($cmd) == 0) or die("Command '$cmd' failed: $!");
}

sub check_pairs {
    my $sample = shift;
    my @fqs = glob "$basedir/$sample/reads.processed.solexaqa/*.fq $basedir/$sample/reads.processed.solexaqa/*.fq.gz";
    foreach (@fqs) {
        if (/\.[12]\.fq(\.gz)$/) {
            return 1;
        }
    }
    return 0;
}

sub read_samples {
    if (!$sample_file) {
        my @inputs = ();
        if ($ARGV[0]) {
            (-f $ARGV[0])
                or usage(1);
            (!$ARGV[1]) or (-f $ARGV[1])
                or usage(1);
            for my $i (1..32768) {
                my $dir = "motus.processing.$i";
                if (-d $dir) { next }
                $sample_file = $dir;
                my $gz = '';
                if ($ARGV[0] =~ /\.gz$/) { $gz = ".gz"; }
                my $dot1 = '';
                my $dot2 = '';
                if ($ARGV[1]) {
                    $dot1 = ".1";
                    $dot2 = ".2";
                }
                mkdir_or_die($dir);
                @samples = ($dir);
                my $arg0abs = abs_path($ARGV[0]);
                system_("ln -s $arg0abs $dir/lane1$dot1.fq$gz");
                $IS_PAIRED_END = 0;
                if ($ARGV[1]) {
                    $IS_PAIRED_END = 1;
                    my $arg1abs = abs_path($ARGV[1]);
                    system_("ln -s $arg1abs $dir/lane1$dot2.fq$gz");
                }
                last;
            }
        } else {
            usage(1);
        }
    } else {
        (-s $sample_file)
            or die("ERROR & EXIT: Sample file missing or incorrect (was '$sample_file'). Specify sample file using --sample_file\n");
        open( SAMPLE, '<', $sample_file) or die("Cannot open sample file ($sample_file): $!");
        while (<SAMPLE>) {
            chomp $_;
            unless ( $_ =~ m/^$/ ) {
                if ( $_ =~ m/^\S+$/ ) {
                    $_ =~ s/ //g;
                    (-d $_)
                        or die("Sample $_ does not seem to be a directory.");
                    push( @samples, $_ );
                    my $pair_endiness = check_pairs($_);
                    #if ($IS_PAIRED_END == -1) { $IS_PAIRED_END = $pair_endiness; }
                    #elsif ($IS_PAIRED_END != $pair_endiness) {
                    #    die("All samples must be either single end or all samples must be paired-end.");
                    #}
                }
            }
        }
    }
}


sub detect_format {
    # format = detect_format(sample)
    #
    #
    # Returns one of 's' (sanger), 'x' (solexa), or 'i' (solill)
    my $lane = shift;
    if ($lane =~ /\.gz$/) {
        open(LANE, "gunzip -c $lane|")
            or die("Cannot open gunzip pipe to read '$lane': $!");
    } else {
        open(LANE, '<', $lane)
            or die("Cannot open '$lane': $!");
    }
    my $NR = 0;
    my $is_solexa = 0;
    my $is_solill = 0;
    while (<LANE>) {
        ++$NR;
        chomp;
        if ( ($NR % 4) == 0 ) {
            if (/[!-:]/) { return 's' }
            elsif (/[;-?]/) { $is_solexa = 1; }
            elsif (/[@-h]/) { $is_solill = 1; }
        }
        if ($NR == 1000) { last }
    }
    close(LANE);
    if ($is_solexa) { return 'x'; }
    if ($is_solill) { return 'i'; }

    die "\nERROR & EXIT: Unknown format of files in sample $lane. Please use the --fastq-format argument to $0.\n";
}


# READ TRIM FILTER

sub read_trim_filter {
	foreach my $sample (@samples) {
        &mkdir_or_die("$temp_dir/$sample/temp");
        &mkdir_or_die("$basedir/$sample/stats");
        &mkdir_or_die("$basedir/$sample/reads.processed.solexaqa");

        my $file_formats_array = '';
        my @files = glob("$sample/*.fq $sample/*.fq.gz");
        foreach my $lane (@files) {
            if ($INPUT_FORMAT eq 'sanger') { $file_formats_array .= 's'; }
            elsif ($INPUT_FORMAT eq 'auto') {
                $file_formats_array .= &detect_format($lane);
            } else { $file_formats_array .= 'i'; }
        }
        my $trim5prime = 'yes';

        my $is_paired_end_str = ($IS_PAIRED_END ? 'yes' : 'no');
		my $cmd = "$src_dir/MOCATReadTrimFilter_aux.pl " .
                    "-sample $sample " .
                    "-trim_5prime_end $trim5prime " .
                    "--solexaqa_or_fastx solexaqa " .
                    "-paired_end_data $is_paired_end_str " .
                    " -length_cutoff $length_cutoff " .
                    "-qual_cutoff $quality_cutoff " .
                    "-src_dir $src_dir " .
                    "-bin_dir $bin_dir " .
                    "-cwd $basedir " .
                    "-temp_dir $temp_dir " .
                    "-zcat $zcat " .
                    "-use_5prime_file no " .
                    "-file_formats_array $file_formats_array";
        system_($cmd);
	}
}

sub check_index {
    my $ref = shift;
    if (-e "$data_dir/$ref.index.ann") { return; }
    my $bin_dir = abs_path($bin_dir);
    system_("$bin_dir/2bwt-builder $data_dir/$ref");
}

sub screen {
    my $database = shift;
    my $reads = shift;

    check_index($database);
    my @screen = ($database);

    my $mapping_mode;
	if ( $conf{MOCAT_mapping_mode} eq 'unique' ) {
		$mapping_mode = "-r 0";
	}
	elsif ( $conf{MOCAT_mapping_mode} eq 'random' ) {
		$mapping_mode = "-r 1";
	}
	elsif ( $conf{MOCAT_mapping_mode} eq 'allbest' ) {
		$mapping_mode = "-r 2";
	}
	else {
		die "ERROR & EXIT: Unknown MOCAT_mapping_mode";
	}
	foreach my $sample (@samples) {
        mkdir_or_die("$temp_dir/$sample/temp");
        mkdir_or_die("$temp_dir/$sample/stats");
        my $inputfile = "$temp_dir/$sample/temp/SOAP_INPUT.gz";
        if ( $systemType =~ m/Darwin/ ) {
            system_("$zcat $cwd/$sample/$reads/*pair*gz $cwd/$sample/$reads/*single*gz > $inputfile");
        } else {
            system_("cat $cwd/$sample/$reads/*pair*gz $cwd/$sample/$reads/*single*gz > $inputfile");
        }

		foreach my $screen (@screen) {
			# Define variables
            my $db_on_db    = $screen;
			my $sf          = "$cwd/$sample/reads.screened.$db_on_db.$conf{MOCAT_data_type}";
			my $ef          = "$cwd/$sample/reads.extracted.$db_on_db.$conf{MOCAT_data_type}";
			my $mf          = "$cwd/$sample/reads.mapped.$db_on_db.$conf{MOCAT_data_type}";
			my $file_output = "$mf/$sample.mapped.$reads.on.$screen.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}";
			my $screen_index    = "$data_dir/$screen.index";
			my $ids_file    = "$mf/$sample.aligned.$reads.on.$database.ids";
			my $stats_file  = "$cwd/$sample/stats/$sample.reads.screened.$db_on_db.$conf{MOCAT_data_type}.stats";
			my $estats_file = "$cwd/$sample/stats/$sample.reads.extracted.$db_on_db.$conf{MOCAT_data_type}.stats";

			# Get lane IDs
			my @lanes = glob("$cwd/$sample/$reads/*pair.1.fq.gz");
			foreach my $i ( 0 .. ( scalar @lanes - 1 ) ) {
				chomp( $lanes[$i] );
                $lanes[$i] = File::Basename::basename($lanes[$i]);
				$lanes[$i] =~ s/.pair.1.fq.gz//;
			}

			# Get lanes in the sample folder, for fasta file screen
			chomp( my @fqs = `ls -1 $cwd/$sample/*.fq $cwd/$sample/*.fq.gz 2>/dev/null | grep -v 'trimmed\.filtered' | grep -v '\.single\.fq' | grep -v '\.pair\.'| grep -v "\\.2\\.fq"` );
			foreach my $i ( 0 .. ( scalar @fqs - 1 ) ) {
				chomp( $fqs[$i] );
				$fqs[$i] =~ s/$cwd\/$sample\///;
			}

            mkdir_or_die($mf);

            my $cmd =
                "$bin_dir/soap2.21 " .
                "-a $inputfile " .
                "-o $file_output.soap.tmp " .
                "-D $screen_index $mapping_mode " .
                "$conf{screen_soap_cmd} " .
                "-l " . $conf{screen_soap_seed_length} .
                " -v " . $conf{screen_soap_max_mm} .
                " -p $processors";
            system_($cmd);
            (-e "$file_output.soap.tmp") or die("SOAP failed");
            system_("mv $file_output.soap.tmp $file_output.soap");
            mkdir_or_die($sf);
            mkdir_or_die($ef);

            open(INPUT, "cat $file_output.soap | perl $src_dir/MOCATFilter_remove_in_padded.pl -format SOAP -db $data_dir/$database -step screen -ids $ids_file |")
                or die("Cannot open soap file '$file_output.soap': $!");
            my %inserts;
            while (<INPUT>) {
                chomp;
                my @F = split "\t";
                my $len = $F[5];
                my $mm = $F[-1] =~ tr/[a-zA-Z]/[a-zA-Z]/;
                my $as = 100-($mm/$len)*100;
                if ($as >= $identity_cutoff && $len >= $length_cutoff) {
                    $F[0] =~ m#(.+)/[12]#;
                    $inserts{$1}++;
                }
            }
            close(INPUT);
            open(OUTPUT, '>', $ids_file)
                or die("Cannot open '$ids_file' for writing: $!");
            foreach my $i (keys %inserts) {
                print OUTPUT "$i/1\n$i/2\n";
            }
            close(OUTPUT);

            $cmd =
			    "$src_dir/MOCATScreen_filter.pl " .
                "-zip gzip " .
                "-ziplevel 6 " .
			    " -in ";
			foreach my $lane (@lanes) {
				$cmd .= " $cwd/$sample/$reads/$lane ";
			}
            $cmd .= " -out ";
            foreach my $lane (@lanes) {
                    $cmd .= " $sf/$lane.screened.$screen ";
                }
			$cmd .= " -ex ";
            foreach my $lane (@lanes) {
				$cmd .= " $ef/$lane.extracted.$screen ";
            }
            $cmd .=
			    " -toremove $ids_file " .
                "-stats $stats_file " .
                "-estats $estats_file " .
                "-identity $identity_cutoff " .
                "-length $length_cutoff " .
                "-soapmaxmm $conf{screen_soap_max_mm}";

            system_($cmd);
			# If saving in SAM format
            system_("$src_dir/MOCATExternal_soap2sam.pl < $file_output.soap | gzip > $file_output.sam.gz ");

			# Zip soap file
            system_("gzip -f $file_output.soap");
            unlink("$ids_file");

		}
        unlink($inputfile);
	}    # End, each sample
}

sub filter {
    my $database = shift;
    my $reads = shift;

	foreach my $sample (@samples) {
        &mkdir_or_die("$temp_dir/$sample/temp");

		my $tosort_file = "$temp_dir/$sample/temp/tosort.tmp";
        my $stats = "$cwd/$sample/stats/$sample.$reads.stats";
        if ($reads eq 'reads.processed.solexaqa') {
            $stats = "$cwd/$sample/stats/$sample.readtrimfilter.solexaqa.stats";
        }
        (-e $stats)
            or die("Missing stats file '$stats'.");

        my $len_file = "$temp_dir/$sample/temp/lengths.tmp";
		system_("cat $data_dir/$database.len | sort -u > $len_file");
        my $output_folder = "$cwd/$sample/reads.filtered.$database.$conf{MOCAT_data_type}";
        my $output_file   = "$output_folder/$sample.filtered.$reads.on.$database.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff";

        if (!(-e "$data_dir/$database.len") ) {
            system_("$src_dir/MOCATFilter_falen.pl -infile $data_dir/$database -outfile $data_dir/$database.len");
        }

        my $input_folder = "$cwd/$sample/reads.mapped.$database.solexaqa";
        my $file = "$sample.mapped.$reads.on.$database.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.soap.gz";
        ( -e "$input_folder/$file" )
            or die "\nERROR & EXIT: Missing $input_folder/$file";
        (-e $len_file)
            or die("Length file '$len_file' missing.\n");
        my $input_file = "$input_folder/$file";
        mkdir_or_die("$temp_dir/$sample/temp/psort");
        mkdir_or_die("$output_folder");
        system_(
            "$zcat $input_file | " .
            "perl $src_dir/MOCATFilter_remove_in_padded.pl -format SOAP -db $data_dir/$database -step filter | " .
            "$src_dir/MOCATFilter_soap2sam.awk " .
            "-v MIN_LEN=$length_cutoff " .
            "-v MIN_AS=$identity_cutoff " .
            "> $tosort_file ");
        my $cmd =
            "LC_ALL=C " .
                "sort " .
                #"--parallel=$processors " .
                "--buffer-size=$conf{filter_psort_buffer} " .
                "--temporary-directory=$temp_dir/$sample/temp/psort  " .
                $tosort_file;
        if ($IS_PAIRED_END) {
            $cmd .= "| $src_dir/MOCATFilter_filterPE.pl ";
        }
        $cmd .=
			"| $bin_dir/msamtools " .
                "-S -m filter " .
                "-l $length_cutoff " .
                "-p $identity_cutoff " .
                "-z 0 " .
                "--besthit " .
                "-t $len_file " .
                "- | $src_dir/MOCATFilter_stats.pl " .
                "-format SAM " .
                "-length $length_cutoff " .
                "-identity $identity_cutoff " .
                "-stats $stats " .
                "| $bin_dir/msamtools " .
                    "-Sb -m merge " .
                    "-t $len_file " .
                    "- > $output_file.bam.tmp";
        system_($cmd);
        unlink($tosort_file);
        (-e "$output_file.bam.tmp") or die("Failed to create $output_file.bam.tmp");
        system_("mv $output_file.bam.tmp $output_file.bam");
        unlink($len_file);
        unlink($input_file);
	}    # End loop samples

}

sub profile {
    my $taxo_profiling_mode = shift;
    my $reads = shift;
    my $database = shift;

    my $rownames_dir = $conf{MOCAT_rownames_dir};
    &mkdir_or_die($rownames_dir);
    $rownames_dir = abs_path($rownames_dir);
    my $sample_file_basename = chomp($sample_file);

    my $taxo_profiling_map            = "$data_dir/$database.refmg.map";
    my $taxo_profiling_map_tot_len    = "$data_dir/$database.refmg.map.total_length";
    my $taxo_profiling_motu_map       = "$data_dir/$database.motu.map";
    my $taxo_profiling_functional_map = "$data_dir/$database.functional.map";
    my $databases = $database;
	my $profile_type;

	if ( $taxo_profiling_mode eq 'RefMG' ) {
		$profile_type = "taxonomic";
	}
	elsif ( $taxo_profiling_mode eq 'mOTU' ) {
		$taxo_profiling_map = $taxo_profiling_motu_map;
		$profile_type       = "motu";
	} else {
        die("BUG\n");
	}

    ( -e "$taxo_profiling_map" )
        or die("\nERROR & EXIT: Missing map file $taxo_profiling_map");

    my $counter = -1;
    foreach my $sample (@samples) {
        ++$counter;
        mkdir_or_die("$temp_dir/$sample/temp");
        my $output              = "$sample.filtered.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff";
        my $outputTax           = "$sample.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff";
        my $outputTaxRownames   = "$sample_file_basename.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff";
        my $folder_filtered     = "$cwd/$sample/reads.filtered.$databases.$conf{MOCAT_data_type}";
        my $fullname = "NOFILE";

        my $folder_tax      = "$profile_type.profiles.$databases.$conf{MOCAT_data_type}";
        my $input;
        if ( $conf{MOCAT_mapping_mode} eq "allbest" ) {
            $input = "$folder_filtered/$output.bam";
        }
        elsif ($conf{MOCAT_mapping_mode} eq "unique"
            || $conf{MOCAT_mapping_mode} eq "random" )
        {
            $input = "$folder_filtered/$output.soap.gz";
        } else {
            die("Unknown MOCAT_mapping_mode: '$conf{MOCAT_mapping_mode}'.");
        }
        (-e $input)
            or die("\nERROR & EXIT: Missing filtered mapping results file: $input");
        my $inserts;
        my $padded_stats_file;
        if ($calculateTaxonomy_previous_calc_coverage_stats_file) {
            $inserts = "$cwd/$sample/stats/$sample.extracted.screened.reads.processed.solexaqa.on.mOTU.v1.padded.after.PE.filter.and.within.padded.region.solexaqa";
        }
        else {
            $inserts = "$cwd/$sample/stats/$sample.readtrimfilter.$conf{MOCAT_data_type}";
        }
        (-e "$inserts.stats")
            or die("ERROR & EXIT: Missing stats file: $inserts.stats");
        my $covfile = "$cwd/$sample/stats/$sample.coverage.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.stats";
        
        $padded_stats_file = "$cwd/$sample/stats/$sample.extracted.screened.$reads.on.$databases.after.PE.filter.and.within.padded.region.$conf{MOCAT_data_type}.stats";
        if ( $reads eq 'reads.processed' ) {
        	$padded_stats_file = "$cwd/$sample/stats/$sample.extracted.$databases.after.PE.filter.and.within.padded.region.$conf{MOCAT_data_type}.stats";
        }
        
        
        my $mmfile = "$cwd/$sample/stats/$sample.coverage.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.multiplemapper.stats";
		my $folder_coverage = "coverage.$databases.$conf{MOCAT_data_type}";

        mkdir_or_die("$cwd/$sample/base.$folder_coverage");
        mkdir_or_die("$cwd/$sample/insert.$folder_coverage");
        mkdir_or_die("$cwd/$sample/$folder_tax/$outputTax");
        my $is_paired_end_str = ($IS_PAIRED_END ? 'yes' : 'no');
        my $cmd =
                "$src_dir/MOCATCalculateTaxonomy.pl " .
                "-zcat $zcat " .
                "-bin $bin_dir " . 
                "-i $input " .
                "-s $inserts " .
                "-sample $sample " .
                "-cwd $cwd " .
                "-rcr $reads " .
                "-dt solexaqa " .
                "-taxrownames $rownames_dir/$outputTaxRownames " .
                "-rownames $rownames_dir/$databases.rownames " .
                "-out $folder_coverage/$output " .
                "-taxout $cwd/$sample/$folder_tax/$outputTax " .
                "-match $conf{MOCAT_mapping_mode} " .
                "-datadir $data_dir " .
                "-pos $databases " .
                "-file_list no " .
                "-file $fullname " .
                "-falen $src_dir/MOCATFilter_falen.pl " .
                "-covfile $covfile " .
                "-mmfile $mmfile " .
                "-map $taxo_profiling_map " .
                "-len $taxo_profiling_map_tot_len " .
                "-mode $taxo_profiling_mode " .
                "-PE_filter $is_paired_end_str " .
                "-counter $counter " .
                "-padded_stats_file $padded_stats_file";
       system_($cmd); 
    }
}



sub paste1 {
    my $taxo_profiling_mode = shift;
    my $reads = shift;
    my $database = shift;
    my $sample_file_basename = chomp($sample_file);
    

    chomp( my $username = `whoami` );
    chomp( my $hostname = `hostname` );

    my $PUBLIC = 0;

    my $databases = $database;
    my $assembly_type = 'assembly';
    my $end;
    my $profile_type;
    my @levels;
    if ( $taxo_profiling_mode eq 'mOTU' ) {
        @levels       = ('mOTU');
        $profile_type = "motu";
    }
    if ( $taxo_profiling_mode eq 'RefMG' ) {
        @levels = ( 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'specI_clusters', 'taxaid' );
        $profile_type = "taxonomic";

    }
    if ( $taxo_profiling_mode eq 'identifier' ) {
        @levels       = ('identifier');
        $profile_type = "identifier";
    }
    if ( $taxo_profiling_mode eq 'functional' ) {
        @levels       = ('cog', 'ko', 'module', 'pathway');
        $profile_type = "functional";
    }

    &mkdir_or_die("$cwd/$profile_type.profiles/$databases/" .
            "$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff");
    &mkdir_or_die("$cwd/$profile_type.profiles/$databases/" .
            "$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff");

    my @BOI = ( 'insert' );
    my @NORM = ('mm.dist.among.unique.scaled') ;
    if ($PUBLIC == 0) {
        @BOI = ( 'base', 'insert' );
        @NORM = ( 'raw', 'norm', 'only.unique.raw', 'only.unique.norm', 'scaled', 'only.unique.scaled', 'mm.dist.among.unique.raw', 'mm.dist.among.unique.norm', 'mm.dist.among.unique.scaled' );
    }
    
    foreach my $boi ( @BOI ) {
        foreach my $norm ( @NORM ) {
            foreach my $i (@levels) {

            	my $to_paste = '';
            	my $rownames;
            	my $folder;
            	my $file;
            	my %hash;
            	my $counter1 = 0;
            	my $counter2 = 0;
            	foreach my $sample (@samples) {
            		$counter1++;
                    $folder = "$cwd/$sample/$profile_type.profiles.$databases.$conf{MOCAT_data_type}";

                    $file     = "$folder/" .
                            "$sample.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.$boi.$norm.$i";
                    $rownames = "$folder/$sample.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.$i.rownames";
            		( -e $file )
                        or die "\nERROR & EXIT: Missing $file";
            		unless ( $hash{$counter2} ) {
            			$hash{$counter2} = "";
            		}
            		$hash{$counter2} = $hash{$counter2} . " $file ";
            		if ( $counter1 == 100 ) {
            			$counter1 = 0;
            			$counter2++;
            		}
            		print ".";
            	}
            	my $norm2 = $norm;
            	if ( $norm eq 'count' ) {
            		$norm2 = 'raw';
            	}
            	( -e $rownames )
                    or die("\nERROR & EXIT: Missing file '$rownames'");

            	my $name     = "$cwd/$profile_type.profiles/$databases/" .
                                        "$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff/" .
                                        "$sample_file_basename.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.$boi.$norm.$i";
            	my $COG_name = "$cwd/$profile_type.profiles/$databases/" .
                                        "$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff/" .
                                        "COGs/$sample_file_basename.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.$boi.$norm.$i";
            	if ( $taxo_profiling_mode eq 'mOTU' ) {
            		system_("mkdir -p $cwd/$profile_type.profiles/$databases/$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff/COGs");
            	}
            	system_("cp $rownames $name.tmp1");
            	foreach my $f ( sort { $a <=> $b } keys %hash ) {
            		system_("paste $name.tmp1 $hash{$f} > $name.tmp2 && mv $name.tmp2 $name.tmp1");
            		print ".";
            	}
            	chomp( my $currentDate = `date` );
            	open C, ">$name.comments" or die "ERROR & EXIT: Cannot write to $!";
                print C "# MOCAT $profile_type abundance table ($boi.$norm) summarized at $i level\n";
                print C "# Created $currentDate by $username @ $hostname\n";
                print C "# $cwd/motus.pl " . join( " ", @ARGV ) . "\n";
            	print C "# Original filename: $name\n";
            	close C;
            	system_("cat $name.comments $name.tmp1 > $name && rm -f $name.comments $name.tmp1");

            	my @mg = ( 'COG0012', 'COG0049', 'COG0052', 'COG0048', 'COG0016', 'COG0018', 'COG0080', 'COG0088', 'COG0081', 'COG0087', 'COG0090', 'COG0085', 'COG0091', 'COG0092', 'COG0093', 'COG0094', 'COG0096', 'COG0097', 'COG0098', 'COG0099', 'COG0100', 'COG0102', 'COG0103', 'COG0124', 'COG0172', 'COG0184', 'COG0185', 'COG0186', 'COG0197', 'COG0200', 'COG0201', 'COG0202', 'COG0215', 'COG0256', 'COG0522', 'COG0495', 'COG0533', 'COG0525', 'COG0552', 'COG0541' );
            	if ( $taxo_profiling_mode eq 'mOTU' ) {
            		print localtime() . ": Grouping $i by COGs...";
            		foreach my $i (@mg) {
            			open IN,  "$name"         or die "ERROR & EXIT: Missing $name";
            			open OUT, ">$COG_name.$i" or die "ERROR & EXIT: Cannot write to $COG_name.$i\n";
                        print OUT "# MOCAT $version $profile_type abundance table ($boi.$norm) summarized for $i\n";
                        print OUT "# Created $currentDate by $username @ $hostname\n";
                        print OUT "# $cwd/motus.pl " . join( " ", @ARGV ) . "\n";
            			print OUT "# Original filename: $COG_name.$i\n";

            			while (<IN>) {
            				if ( $. == 5 || $_ =~ /^$i\./ ) {    ### WE HAVE 5 HEADER ROWS ###
            					print OUT $_;
            				}
            			}
            			close IN;
            			close OUT;
            			system_("gzip -f $COG_name.$i");
            		}
            		print " OK!\n";
            	}
            	system_("$src_dir/MOCATFraction.pl -in $name -out $name.fraction");
            	system_("gzip -f $name $name.fraction");
            }
        }
    }

    &mkdir_or_die($OUTPUT_FOLDER);
    my @input_files = ();
    my @output_files = ();
    
    
    if ($taxo_profiling_mode eq 'mOTU'){
        my $name = "$cwd/$profile_type.profiles/$databases/" .
                        "$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff/" .
                        "$sample_file_basename.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.insert.mm.dist.among.unique.scaled.mOTU.gz";
        my $pre = "$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff";
        my $file = "$sample_file_basename.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff";
        my $db = "$cwd/motu.profiles/$databases";
        &mkdir_or_die("$db/$pre");
        system_("perl $src_dir/MOCATPasteTaxonomyCoverageFiles_generate_mOTU_tables.pl " .
                "--wd=$cwd --map=$data_dir/$databases.motu.linkage.map --table='$name' --prefix='$db/$pre/$file' --prevalence=$prevalence_minimum");
        system_("gzip -f $db/$pre/$file*tab");
        # Make easy output
        @input_files = ("$db/$pre/$file.annotated.mOTU.clusters.fraction.tab.gz", "$db/$pre/$file.mOTU.clusters.fraction.tab.gz",
                        "$db/$pre/$file.annotated.mOTU.clusters.tab.gz",           "$db/$pre/$file.mOTU.clusters.tab.gz");
		@output_files = ('annotated.mOTU.abundances.gz', 'mOTU.abundances.gz', 'annotated.mOTU.counts.gz', 'mOTU.counts.gz')
    }
    
    if ($taxo_profiling_mode eq 'RefMG'){
        @input_files = ("$cwd/$profile_type.profiles/$databases/$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff/$sample_file_basename.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.insert.mm.dist.among.unique.scaled.species.fraction.gz",
                        "$cwd/$profile_type.profiles/$databases/$reads.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff/$sample_file_basename.$profile_type.profile.$reads.on.$databases.$conf{MOCAT_data_type}.$conf{MOCAT_mapping_mode}.l$length_cutoff.p$identity_cutoff.insert.mm.dist.among.unique.scaled.species.gz");
		@output_files = ('NCBI.species.abundances.gz', 'NCBI.species.counts.gz')
    }
    
    if (scalar @input_files > 0) {
        for my $i (0 .. scalar @input_files - 1) {
            system_("ln -sf $input_files[$i] $OUTPUT_FOLDER/$output_files[$i]");
        }
    }
}

sub check_or_expand_data {
    if (-d "$cwd/motus_data") { return }

    my $motus_dir;
    my $arg0dir = dirname abs_path($0);
    if (-d "$arg0dir/motus_data") {
        $motus_dir = "$arg0dir/motus_data";
    } elsif (make_path("$arg0dir/motus_data")) {
        $motus_dir = "$arg0dir/motus_data";
        expand_data($motus_dir);
    } else {
        $motus_dir = "$cwd/motus_data";
        expand_data($motus_dir);
    }
    $data_dir = "$motus_dir/data";
    $src_dir  = "$motus_dir/src";
    $bin_dir  = "$motus_dir/bin";
}

sub check_executable {
    my $bin = shift;
    ( -x "$bin_dir/$bin" ) or die "$bin was not found in $bin_dir\n";
}

sub pre_check_files_bins {
    check_executable("fastx_quality_stats");
    check_executable("fastq_trim_filter_v5_EMBL");
    check_executable("soap2.21");
    check_executable("2bwt-builder");
    check_executable("msamtools");
    for my $sample (@samples) {
        my @files = glob("$sample/*.fq $sample/*.fq.gz");
        if ($IS_PAIRED_END) {
            my %lanes = ();
            for my $fname (@files) {
                $fname =~ s/\.[12]\.fq(\.gz)?//;
                ++$lanes{$fname};
            }
            for my $v (values %lanes) {
                if ($v != 2) {
                    die "Input files are not paired up in sample $sample.\n";
                }
            }
        } else {
            if (! (scalar @files)) {
                die "No input files in sample $sample";
            }
        }
    }
}

sub check_args {
    (grep {$_ eq $INPUT_FORMAT} ('auto', 'sanger', 'illumina'))
        or die("fastq-format must be one of 'auto' [the default], 'sanger', or 'illumina' [for older Illumina machines, newer models use sanger format]. Got $INPUT_FORMAT.\n");
}

sub usage {
    my $hard_exit = shift;
    print <<EOF;

For single-sample mode:

If your data is single-end:
$0 <INPUTFILE>

If your data is paired-end:

$0 <INPUTFILE1> <INPUTFILE2>

For processing multiple samples at once, you must set up the following:

1.  Create a directory for each sample with the read files inside (you can even
    have more than one lane per sample).
2.  Create a text file <SAMPLE_FILE> with the directory names (one per line)
3.  Pass the name of that text file to the ``$0`` script::

$0 --sample-file <SAMPLE_FILE>

This is the same organization is used by the MOCAT pipeline.

The most important results are saved in the RESULTS folder:

    RESULTS/annotated.mOTU.abundances.gz
    RESULTS/mOTU.abundances.gz
    RESULTS/NCBI.species.abundances.gz

See webpage for details: http://www.bork.embl.de/software/mOTU/

Citation: Sunagawa et al.. Nature Methods 10, 1196-1199 (2013)
		http://dx.doi.org/10.1038/nmeth.2693

Options

    --verbose           Be more verbose while running the analysis

    --processors=N      This should be an integer and defines the number of
                        processors that the script will use.

    --length-cutoff=L   The minimum size per read (after quality-based
                        trimming), default: 45.

    --identity-cutoff=I Minimum percentage identity in alignment (default: 97)

    --quality-cutoff=Q  Basepair quality cutoff (default: 20)

    --fastq-format      The format of the input files. Must be one of 'auto'
                        (the default), 'sanger', or 'illumina'. Note that new
                        Illumina machines actually use the 'sanger' format. The
                        auto-detection should generally work well.

    --output-directory  Where to place the final results file (by default it
                        uses a directory named ``RESULTS``).

Website: http://www.bork.embl.de/software/mOTU

EOF

    if ($hard_exit) {
        exit(1);
    }
    exit;
}

sub copy_map_tables {
    system_("cp $data_dir/mOTU.v1.map.txt $OUTPUT_FOLDER");
    system_("cp $data_dir/mOTU-LG.v1.annotations.txt $OUTPUT_FOLDER");
}

sub main {
    check_args();
    read_samples();
    check_or_expand_data();
    pre_check_files_bins();
    #read_trim_filter();
    screen('mOTU.v1.padded', 'reads.processed.solexaqa');
    filter('mOTU.v1.padded', 'reads.processed.solexaqa');
    profile('mOTU', 'reads.processed.solexaqa', 'mOTU.v1.padded');
    paste1('mOTU', 'reads.processed.solexaqa', 'mOTU.v1.padded');

    $calculateTaxonomy_previous_calc_coverage_stats_file = 1;
    screen('RefMG.v1.padded', 'reads.extracted.mOTU.v1.padded.solexaqa');
    filter('RefMG.v1.padded', 'reads.extracted.mOTU.v1.padded.solexaqa');
    profile('RefMG', 'reads.extracted.mOTU.v1.padded.solexaqa', 'RefMG.v1.padded');
    paste1('RefMG', 'reads.extracted.mOTU.v1.padded.solexaqa', 'RefMG.v1.padded');

    copy_map_tables();
}

GetOptions(
    'sample-file=s'             => \$sample_file,
    'dbdir=s'                   => \$data_dir,
    'processors=i'              => \$processors,
    'length-cutoff=i'           => \$length_cutoff,
    'identity-cutoff=i'         => \$identity_cutoff,
    'quality-cutoff=i'          => \$quality_cutoff,
    'prevalence-min=i'          => \$prevalence_minimum,
    'fastq-format=s'            => \$INPUT_FORMAT,
    'output-directory=s'        => \$OUTPUT_FOLDER,
    'verbose'                   => \$verbose,
) or usage(1);


&main();
1;

__DATA__
H4sIAAE3z1IAA+w8bVvbxrL9in/FVjjBDrZsixBaDARKHC7n8lZwenIbc1xhyyBiS44kE0jC+e13
ZvZFu5J5S1POvW39PBhpd3Z2dnZ2dmZ21nHUq333B3/q8FlaWoT/i0594QW9Ly49p//i813DcV4s
vFhYWEC4Rv3F4tJ3bPGPJgw/kzhxI8a+64Xe8Cy8Ge4kjN4/Bj2P/Ilh/nf3Nzfah57bb0f+6LU/
TLyo604u7fHw2/Rx+/w7S88XFjPzv7i4WP+O1b9N97d//uLzP/t9bRJHtRM/qHnBBRt70bAwiT0W
J5HfS5r0/NGNAj84jfnblpeE42R5eScMTpuFwixrn/kx64V9j8H/sRslLByw5MxjJFfMDdzhVYxV
/tgb+oEHTTYFdKlXhufxVeSfniWstfvTDoD32U9b2xXm1BtOFb4Wsl1E3tBzY6/PJkHfi9jW3hu2
dbDDLhbswuiKldiT2L3w+t0PE3fY7U2ScDCoyLKBH8VJt3fVG3oVVhx6wWlypmCKRovi2PUjaOIF
/W7fTVwogdF3F8ewRjxARBiASSNZBIBQAuup2/cjeAKWiqc4HHqX7ge3G0bdgRsnl1DW+4jQiTca
CyDE2B2E0chN4q4bRe4VFP4KDEQE7mis9ddVHbJyE8dcbO+3R+4l0z6rrK6q+DinVsWT0U2t/CD2
oiTOVeEIPnSJkgHpCqiycoXdi8UuzqfVZLPsYhHnDUVi6CZenOAEswsviv0w4Dhp3N1eOAkQoewO
RG1/nABQXCrMzBnTtRrPCbrWWCczkwCrTaWCFLD6LANkZpoFNEFmJQCgMyKgQ2elA6Az4qFDZ2cS
wYXsaAQrcClWACXkahqUFDnElZU5gue4cuII8CCPGYwKK4oqDkbIqgnGByPFGMDycowNCGyKiEOD
Tz03yXdNDUj8cTAk/1kgPhiqKpRJE20PWDQJUFex/aO3FdY7c4NTjwFb3OiKYfeF3lk4GpcYCl18
FQPh7auxB9L22yRwRx6rxr/hmvIHoEcMgH+zUe2VG30EPcnK7HNhZuo6mFJoM6sLxFjNwjWSCCLN
Bh9YGMG3ffqJiIoLs0DPujscvsY3pGYYs2qDtERNDLH2DJplCxCDs1brexe1YDIcsi/sNPLGrHrB
SPZGXt/mZHj9Ob3SjoFHQ88efDCLUeLtud9ooX8tPWbZySdHpxBQF0ACiP8+8LjObJvFPXcIe3Da
YZU1OJPFbBVlzbuif4zzM4P0/QMIi8dDP2G1Tq1WMaEAxmwFk/OPd9XGsZiHV7CJJB4Lh30+BezE
A7pw34MNDJhTgJUZJMzKm0VsmR1t7B7stKRe7gQwubeCbx5uH7TZL63Do+39PbZwZ4PNndbG3psD
eIpGrDpI9wnFVyyoPbOzkwwT+qzJRu8BlFXH5lTAyJK4yTFGZlUERMT2OAp7Xgwbq51XETchvbsl
DbbAlxKzvmI80Fy2vmVgGtTXjfDObu6FAiTrFFY47jF+coXSlPhx4vdilHrP7Z2R5J+zkhJ2knOu
kLhqrHDlhJorglWwigiF1iwJ0S9+8se464LmhD5nuLY6RyVV68ASLJKKUlCnkwCfqj2LXWeBYXlm
oE8I2pHgt4npVqvNfn6zsbPd/h921N5oH0HZCIy01TxrGB/PqhwXH8+qGDITCh7IWmatt63NN+0W
J8iYAaj+ouyqGuHtCk53SQhYNbxJtornNm38BGdjkeQxCeiMnPvH6hX7vJW3IEULjMwDBuYBJ7KY
s2w/W6K/4rl1LWRlQRkVJc10PCe9OXs/FERvrjUpzgO+ABjIr9i8stbbKmwVzW8r75rdLkq0AaBq
v31Uao1kbDTmfWBzVx7YFUTVjN6PIERCRyHt6F6epXJWtLZZlmbov76NnCC8hRpd7rhe+Rb0zIRj
L2BgL82trc1VYA2g+JPzZnNmChphm+5N0ILvozAK8X3LJPJOUjyHPw03l9qZ3jAEv/GtPvKcIwND
t4LQmj70xleyGeRdyWsQBgzNmyr2xrf8M9hiJEEZWz9DD0rcCficZCNCR+dEUFrybxbXYK/i6rR2
QyUq5hvqRIXoqTfqoyZWOmeKVVkdTXErWdXNqS7cbXVmVj+Ya6c6zPjBuO8t1tWCRN32sF2wpoZn
iTHdpucODvc3W0dH7Gh7b2unVW3tvWLmLgDMEGKENl84SdAcxdLfsGwSDHFyS8qGJIBVtigmbqbv
e8wS9lrr8HD/kD0F3NttsMpoQCz2QO0nIchHws5AwpgbXDEaJYhLjHYyvTHkPuPct9lROPKSM6z7
GIXB6Ut24EVn7pj7uMBHhoy8cIcTClgkYcjO/NOzlzSM6zs4MmW299+0D94AydqyK8JA39XRzJ2J
z/xBQiOnxaY5OMK5X5XQXD9oEDDzOkQjDyHDCgLCyUPIOIGAWBAQeghh1QgozKs2Ek4jVLzNE/Wy
XgYwVvVwxjyRDybXjFjCIgiywokW8y9LV6lU8P9abGVhnDAhyhndULhZK8B+IdQCOjKcB7T1YcyC
NeXLyup0v6aZQouHeeZItZd6LFR1zP3Of5Xs+XLHboCiqDU1bdSV6ojrx9w2zJFg3UP222+uhSzw
gfMjA7+YVU9ug0z502ACfsAIJL8TyF1walfXvLHzoMZm7wLFIynPrtKe99afBxvbh1x7ppLxZ9ei
fyu+r1R8e34PGDz03GAyxmA9V3m3hjDe7O5uHG7/2hKOHkbpP0Z+4k1xxeWbjTOOEy7mOy/tNsGT
YN4vQJBxpppsdHEP2AcHO7gjf0RtY6+XyDC1JhOcuzwiPm2WeM1UOWgWyN5GRqLJrSzu38VESwYt
hRyt4mZEonDn2pXa75usYYpz4tDdi1MYNkhUSVAluFRTJJZBkuv2ojzDeD/yIjkIbPwEtkU1CrDk
oayL/ePSBXXgXXgBeMTYCtmLLWrMUZ4elRsoZlJg+j9PbgUsCm8I3sm0eqfJV4wAyNMQ9vspCYLs
eYpd3kRK4wZS6neQ0hCk8DWKosMXatxJfgJ9D/923ctOsnFx2kn+G5p0km0ubJ1AsBs8MzUJ/HmE
LZBm+H5PbTQhpUXJHTfs7Y545T5KV7v1ih22Nl6x9uH2Lnu9vdNuHRIa79JPcISF2Wr6Kcwevfnp
ECzb7b3WkVkTT06m+t7IGW67kLzKpdjldjAPNGjHA1jTkDWmQ2nhHkKr8PXBNizClVvcXgT9eIbd
l1YAeo1PIO0Dcj/FU00ZjG4KtUw6GCtw+wG7UZ5UQPHTpyAW8sNhGhzmnLb8L1+AGFEMNp91cxX4
lLfVgjsqqlkZABh9893B5Adv5ByrbUL47DBeJcZ6C7SDpW9M+qUTcNWC+wJXLnlG0jYDGqEXRhEo
VZvtgPpivTOv997rixDgsuRTpnboBlh3bjO0zlh8Fk6GfXaC6EQska3whmudTrKC4MzFI2MsSkTX
aBRRNZpIaxZfcpGXTKIgEzuAXVIKYRpNy4jfVMlLhY6YVo1hQu4bGpQMJbnc3uu2W7sHSjbvjSSj
SGh2snOz6Qao5Kmje2Pm9uO1WDewcCWBa/cl0CqACIKMfhOaODl04F+M3I+klrvBZFTRXnGeycUp
sXoFtwFosSI4u6Yta1nEd5qZWbk2pV1ctEl3lz/r652p4yerk1jQbVf6ZnyNqGXFvXSyH8DiYmnF
olbR96FKViwdN5VxR23WVk3HAzSIbJOtEuYfV9dylgRxaDQbzFI9gv7h3Rm1K2ktN21Nzs6nrZsG
33QVInjblO+SJKmADzf+CeKzNtWOPLehQ5sbbFyC7ic+90HF5YezCWjIbKeBVcmKlZhlgwVYSnj4
2ABPM9Uo2qwYGoXbpF+pTG6zZIUH+TDF8lUIb1Yy0rbc9bnF+Hvwf/vF/ZjL949aaNOW1l9sJaml
lD8Pedi6MqQBX2jueW6KUQJ0n3Arvyx3nScbFfZkE/624K9N8yYW2EM27T9gT0wFX8i8YZ8a4l7c
+KyknF2nQsotweLm9NoFXrs1vfY5r21Pr+WLZjyJz0ombyspeNkUcjF5JCrCiUxjqzoKZQ+AZ9OF
2SmO/KC7QTriKUoLRcYiDM2VWAenz8BZNppviuabNzXfvLX5lmi+dVPzrVubt0Xz9k3N21OaG0wV
GYRzdfTW57Ce4tXqSHbESgbvuKSQalKeycbn4uiarQl+FriDIUpXBHdl6aYOu2mWClhVuqXDbpml
AlaVtnXYtlm6opgk3Roc2YVLjvMII7Y4eB5HReZkeIOQ0ipSW1kOxXvPG98HxbXmZ4ojS3zHRYq5
tAFGQ85DPygxi8H6MzRNWTlZOjgP/Xf68xRwKc/jYxnJqRmnpEipM8WXMWyOVHqEgizxQHFWjADZ
eleJYjwZCeWMYQL5iKXhoBuD3om8WJYCJxRA0vcuxP7cNMUOc77ee1cxe/KZ+r+WQ5Fd0V5HVZ+L
/jUfluqcTvolYM0g/SH9ZIjHLktap6yaDpeVnz0TzFUjlFRk0dRoW9e4SWlsTUzbSBkCmutDlJQU
snK6zZW0bueJ/WmzZ3hGpE0D4J4CIJI05KQLx1jsiOINCJicgIjdkHFsHiJVlBzzE400GMCRYcIC
z57Q8M9x2LmmOgdZZVb1Z7awoCUgaQgu8wgoxmkimN7Wz7X1sR2pr0zbzNDm5w3Oy7MxdTCGnPxP
J+c/wkfd/zjqRZ4XiJO9b3b1gz633/9oOM7zhcz9j+eLS42/7388xufPdf+jUBA2lB9wHQwL2u97
QQJOAl3EcMdgR4zQj0zCQ28UXniv+U0OPxhPEvEcTpL0xbtMIreXeH3x/us2+M6YJTj0LrwhmP5n
bozdkP8D+zoh2nEDgF3niMSLQsTfi2SxC6Qt9SJPQ9R1DnXvQt3iUCXq8oYqaaUnVLIkPZySJemh
Y908pdJg0oMrjNvLU4dX4KAwa6q6WMZjgsP29t6WzTN+zbsbCezNyG11F8BgPwD4wWr8uXKtLkho
fIRaYKSqplqNs1DtXWYam7zGGwTIX/3+Rcp9bK9VU23LqIbJNq+boAzwchKCVV8MSkkFXcAQsgZN
xbULKXw4XCGUih9KStVFF40cJc50z0F5mdaKwUbgetb1k6cUca1Ts5/xLLMiCuznYveaW6nXhdTP
utc8D0O3D4sPsbDwfUUm0ICCsG0x91LMMIRPQm0I34HrRw7LCLICNUoVqBBEeXXgPHt1QElLenfg
XoM5UMTTEQNgXi4qXO+K58dyTDN8ocheRBjBABVAqWim6QN5KENCEcoo4IB5k9aiCxoNcfIDEqDu
cYh3qnfkyZCIZ6fKjRJ6FNF20bdkoFmjsox0f9bVIB+IpcFQU7SqeGODeDLIDLVJjQ0Y2V6FuX06
UzKGZ4ay75hI3hIzFbyEoYajKAo+JCGri9yXmSlSOcOdgoz8ca8wR5riLDqhvNyYgYeT7OBVII6D
p1mwnY2jNp6etphQwssGfZ6tivXByBF+bc8/tV7vH2pdyg1A704/Gp4hsRL7xvyqsYYFq1Xay6qp
DX4XnRuv261D9jA6HzLpatY16cfTFlOeUcTvN4JUv2irsLpm4C+hrJo9gK/XxxtK2BA2hsj3Yjzq
NPQ9jo8vIZWCSO4OZQOmm8Qpv/zBqn2dhC/ZlU+90srPT88/I58ogaWkQUv+yvNA6OsL2UesqnZC
7NgYrIVdwIROOTIQ0U6RVqQ1ub/AICmwN+nMxS2Dvwp680qLD/2OcRvQxtBbbyvs5qEb7ax7Dd1s
8vDRm7I6lQG5GHGa0EDSJCLGPEosVoSSsDSNRaQ4iFjVesl+Vu7UMHxMeaqyITc4Gteq0Qy3ltVa
JJL0IJwG0dAgriUxvFaiy0mxXJwglCIhQqwViUPraaooiPatt9Ob5/jhMIMWgzp5jksZetyMKxHS
csqPmSJP1TOq+cAFb2ak/T8/HQrVGYU05IvUbbIwz2JOV+sOwlp3U9a6g7SWQVvrJuLoq6iHZ6Zw
ejGVvIzGo/b4l9qzMpKdFxBVhR91MC1WZ2nKAltmXhRBNcKiSuDYyKddZsXvrXLa27XZb07N5LsG
Obt/zwKfj/dZpvbOszzTHDxYYyrBt2jmBmay9NIUvUx+3vXdGXqPmor3uEl3er6mtZa6jlbzK/Lx
kl3Yx58w6e/xd5mWd7S/ccBwwY38eOQmvTM8QC3y7D2Vujclb09aPEkxxZt6jfCsnM9cXl9hmpS0
bhWTVionrZyg/OXFo/Xo8tHiAtJSs92aIiKt3yMj9/SlAQgz+Si/BW3afgiennbb/Q4Er7f3to/+
q/VKRBFk2uh/OkT6p/6o+H/rEqYhwEwCEAMndkff7gzgjt//WnScpezvPz2v//37T4/yuSn+n8b8
KVwfJGByLLPhGUbjf+G/3LMMSr9hNwr3PC84SvqA7KkUL54W3uSnlaDL0LRB2wgD+MW4gQFyRx5H
Mx4Ai/1Pwn8XplXcqK7hzRzwQOaezZGfLoso49qhZ7TnZ/FAImDqBrPAxmNxlw3cKnjTxjE4ZfXL
Rr38UiBbOAZdL+1aKvkRUC6rWhHQu3Q4EieHxMkhcTQkjkSixoeoqkgVjhONyVk2GfeR6BF+9cIw
6vsBPCo2ODobaDNbB+revbDtH/CHT0opCLGFswhos1YtRYNzXFHUVAR7+OG1YMsX4PylU9e6TIdJ
lDK1+xpNfpDDmDJlklbHoHXKFKa0NhStDUFrVSfWmUZs41ZinTyx4hBdyqoUzCcgyLFIAmN4xI7v
JWtsVViH6qic7Pg3sXsKtrmmTNm76viYrbjDwMbSNdjkykTf+sbh1i/c+gEJrtI2ub1XTsW+y6+h
Qsd9bwCOVR+cGejs8/iagGZBJCKPDcNwzJusxw1FJL06+ivwozsE/wtZ2O1NoohyKDvruObg2yE4
GR/g0YG49q7TWFqqdhaWlo6/vOuA0q526o06PTcceH5eP67VTmEGAnBPaFRqoeMdkVKxq7qraEMq
g+dZl9kv64Iw9r1ghXiH6ZFXKAgBvRNdjD3lmiM/pibWcquDsm541pvsoqzSDW8AIhwakKKNGKlA
ZEJgKk53d8pyY+XE0pEH5ZZQTZOlXYremiztthiT9clX1ld1mZFxusqjFDD678hMc7KkNuaXIKWT
z4UskQmGpVonnq+JdD6qFlkW1QanI4FJ/5Ei6DJ2VPtXJ65hyfdFvK1JjdaL6UqbZQP/kqHxW3Vg
y7lkJ5PTGGGw03VoYttO5bltF2eT4zRpn65kcuz9+WJNYKIrbHgTpIApIly2gJ38Mr1WhL/MwGNZ
RZJsbOsHfuK7QzYYuqes9NEfDvEKClfOffptu6gscTSO5VYlX7+gq/GFNVZWkLbnJNJz7txL9gIU
21I5A+pw3ZVyXyP/qfyBIdnkRzGGRjqGBo2rJPechDbCNZa+N47LoFZ5MhCvr2AitQmAKjehi0Gk
Z/xTN5IdLB5rYaFY7Gg2GPeWVEqhgnUEfUvHTV6wIAp+OM4MGzU0H3ny7gVnUXVOjH3kjn+W0M8F
AkQEfH3JFupAa10CZvdJavMC28CuI2hYIhT49IOaq1kGywBXFOZSrqP8W3u7y/4y9PTjsaV0Ml2N
n5tThkhCg+eJVuuXTfHTAA06ogShFMqKI72swDrBBTsoWU/shX6l2MB8XKdclriKXSG3pXcbm1vt
43J1jXL/ar7IjF3HiFwcRglgkz+R44poKvUMxbxLWsyY+vCprNanVbG4xsFh2CAj6NcXr8rVoosT
WPzEK106uLmibLjUBRdtDCkhkMyuShwyoKQZoLF299Xyr8sAaulqov4Xybz6v/FR/p+40DlwYdK+
bfrXXf6f4yzk8r/w39/+3yN8pvl/sxNuvIKB+xoMBvfQAwOzv0OLOUZLtlr1A/r503d4g8el25XH
UBpOEl48DHsuRagxCkS/RREOmKg9/v+cX8azQAbizBivEYuUj4GvLhaLH+jL5kBZgmexJbJ6+HsF
aiTfoIpqxDtWASKLaalH8E5JQLRVwAvpeQxAltj2HmYNVNIf/Ov3FLFfGOYmyMNH/V4JMpJ+22Yg
DwvxjezeNJ6a6WBFoP0qnBzX/pu2QPa/7P17fx1Fku6Pzr/Uq1hj3BszGKG7LGh8OklM4s/gyw8H
PXtOw09b2AI0bUuMJUMzbJ/Xfirj+0RWraskc2mmBzdtS2tVZWVlRkY88TxRWbfjdq/W2ljTpFzo
6D/v95b2XvzS26uqhpTMsFNPD5+PTrDrTxA5/8h93R5vOTp3nD9AITVqrokaq29HVb72jqUDtxul
PTyR6bUd16K/lXj1Y1u5wrjvemZTx7bi7QW9uDnZGMDbXO9bLbojFDWX9JgTsGASuD0ups3k2qWv
j0/joaopelz9HqmACwZxtIvBBUPxOwb4Nf7Mxv9GWBx+/7OFuwvi/8bG3mz874/b/T3+/xp/RvG/
n/HJ21/91iLu65PPHqVy593J5I9np89qWvfsWd/i7T6krb0ztta6g/W9u/cPPrlz//0+GdRv6dH7
++uT7oM75e79yY8Tl8kffPSohurPzz0e/Ugs3Z38MU5/c1JZrEpNnz+/d696vvsf9WHlGS75ZEhs
/Ps3VZd649grESbHdeu3k/7ft956c9JJ7n88uG4/6ebkuHfYIdvXDhz3t/63G9eS5dKng/0oOk0z
eeaVE33e1B1WPmRjff3tG8+evXN9981/6X9Gje+/+aNutnW9z5X/U5Jt//W/TCrvORFs2Ksp87W3
rr0JkdH/vH7tPVK3+GBj91qT/a73gaV+3oeGW/3/9/uJ3Nnpf7hx/6O3NyrD9C+fn6/3/3kW2/9/
q/8oPeqT5ms1Y//PSKH73571bb7s6i5lP/6e4/02/szlf/7Pwzs/Zwp4gf/f2ljfnPX/6xu/+/9f
5c/Vnv/5jcWG6/+3OscK6799fvTdARUO8mf14z989T3/9t8OecAfXV9B2WhyxeT6gbjiPzkTHb/W
3dHfG20UsOApes8Y+LE9wu4Fkb5t5zs3Kof85jvtEQLvI08GxMOVT46fs8tcPMvL8+5fHRw/mYw2
nNhsW4FEM/WJhfGttz0u/lDDpLd34032ZqnVajRZ91c81BZvXx71iP2s/vXlD5MvT8+/mXz1vU9K
32gNHqfRUt/E+GGEv7bna/vjY8dBjvzx+l9fvvVWHdSjvx2f1d33+rbqh1NV9Ty4EY2ot2+O86Wp
yx39MJk7ljJAxaivvv+xHvWyVQJe8Xz6OG5gKOMbiiaXtTkMwtX7U0f6Ep157bV6lUkoTq/5eU1/
CrsYWZIKm+oNvH5Yp9Wlg28Oz5jv8298cVWD8EfJ+/7y28vK7/bmvxbqF0negvY3R+3X3qy8gN/T
8ivMrGAZ+HuqZJw259WGeLEZqkjtIhO8lAHNT/fLq505nuuaG7/WBntJ/0Zlcn9Pe/y9POsf4s9c
/u+bMP2q/P/21uYc/7/9+/v/fp0//2DPfy/AgtdGCNEj1Ji5v3v/0Z1P7dH089AfpEd3HgXpHh/W
nU3nPryX/ncr5Z95nPq90Yvt9LnKbsePc4ut9lLYoyf81h4Y0Jf+3qbhHsYn1vuZe9aaa747ethZ
m1uMnjt+t70Nbv7B4/GZeoq+Ps+sPg5nxic8ga0+vxsPNLdP2uPXanf0/LUrGQtQ+Xij7ngQ6UIY
Ph7UBTuFPEr32otVlqL067yfQ2Yyxt7sAjeC74HfX3t9jIoCwNdPBwA/jeDrd+w2GrvJ+HfaVWv0
ssbgmXTAfpxcq3g4d/L+VFNvtUc43h6VzM8PxYP08JcYi1gJVxmRvRUjcmvBiIy33rvaUFyQOtV6
s/HvtUCxDZLchD+U1LA2lURHY9fSGonvhibcf/BUkxxMLa+YWtFD0+r6bTyMWqg/vj/Y+ICam1eT
P1gKpl8bX2+qrVafP/3IXI9Hx/c+3MSyW5jenUA7Da585tGPd2WPwyc3eNzp+z6PuP7Pb9YnEQZv
8GgDve529QbjzzfHn7fnv6q78K+mnMWjTWjTf+5/+svuFzySxc+tqqf+ujf96634VQ9XndLd95uH
bp/4S338Z99xYTjymn4c3m/3VFOw9m1zxG0/9/ZYySNU0qnRvHEt+xbIbQj1QBqX8vG8oVPejEfE
hodRHv1yT6Ncm7kONvP5ORbT/9ubcT8PtdAIG3onYmtNB8/rSG9/0f8rs+OT3S/4d0//3vpi5hmV
R78/L/Jqfxr+z1Eab4d/Oz05ffbDr/T8R4/w9rZn3/++s7299Tv+/zX+/IPh/7P+u8fnN3jQrHd5
1/8vjwcPn/cBbvQ5j8k+fXpwfvi3g7q5XrwYZgzox9F09Gf6+xFOmv2+99F3T3ItR43NVGYThvkv
x9+fn54fPnVYc3bw+PS7GjgWfK8XEeGInyw5f2H7la47qHTd9E3OnO9HjV8K3r73B7EPZl4lPvr+
yZf90B49P/z66ODrPmEY8Mfw/fjz2fP/wPV16RiBg5Me0I6/nxqf4dvr3x9V4zpim+E6x+OutS9p
/aA9CztskVkfRar/bvrumDrvGdtbfzRV9XVj8qdvD5886VvzuOvvrIh3x899cXPyh4d3Dr/66qjm
SAq3ddfXZ/Fu+c/u3/1/Hn12r//p3j3+fXhH8txNZ+meHX7b//Di5Pg/X1Sn3V/lD89ePD2v20/c
O/z2Wz9s+oM/13cp9Z8+unPnfmxZdvD89PtaoXfWN9F/1c/BTceFbR+0ujT6L/q/P/DtPetP9Pfm
5E/V6u3w+ddH5xtTv236Tmj1IH6vXandvc7feq99vK988oc+O+m9zBMHJdzd0fOzo/iVLdb+8O3p
2XHNKfofNTH+b3+tPxzHAP6p/XR9ymxii7cYXc33eD444ObohZr1zeeA6np+ncmD8x++rXV5Ghu/
oX693RxtIBm91NZufpbfpDYseDLqklvI02N/cEWd8AJcP1hLRp9PT1Z88ODk6Q93v+JFPuu1N9VH
soXcn3wXEjEHXtFZLXXyj/in63OJyWk/FP44Rg0N9Unk+jRr/1EPWHkXyeTR6fPnP6x11+YeaH74
/LQf5foE4ru+JdjXzw+fTb7jKcfJd7tfLdrR7nj+nfWtrrPZWaU+Vhw3a4LLXnU/etN9PaY3y/lG
/ZhqsP0Bzx8/X3IAttwf8uR8Qb/8kLGZv9Gvdd97b747I+cQG/QtvKIWSn+Mr5TZo9jfztdQf0i/
dpY0M7Wqav/7Xvaravroof91wVUqLNbX+DgosbbydNjcddthfkRdlQv7znqts3L63Xw7zMrUUq7z
o1U8dSy3Ob3KmYH5wxcca0yENkWpWxIuanraX/THf3k8d1s6Ppxz3ejw8eH89LITYnXkPrffLmmm
OnyoxSUHaPj8xeELr1G/qUe08Dc3mUNgrA09WzCdNPSszcBcOG4bMc594yTljXj0MlxsIzXu1V8r
WVEzcR2kDQxen3yQ7tUDa6LcxTZAvan7xnG9c/ry6Oxce8YNOxhf60+aKr4efVUb8u/6Dr19FHv5
HcxWTy9+FYcfK6pldJ539JPTHr0+1BrzvR0HZNde6V2Dipdet5Lhawfp/ocH12aiXg1A8c4orbKZ
N22vcL+Rh9ZfY+H0Pr3+qgpw4uXJUT8Z56eTZ4d/rS+LPDs/fTZZ82e+/Ki2b2Pb1MTHy7/yo65N
ybCLu9J7y3MuPGrZK9C9oXdVsr2iiTv/+45C+kTPTHDqJJ6V4Ne183i/0vAi+suddVEH0pMnk2fH
T548rZWTT188O+kHxs8RT7dZ9y7749Dc5ILd2fy0sSX989Ban1ptOlU1HuWLWxNzNW5u2J5sc+H+
ZHNvXBuq6zfjFWjrX3x+vlH3zeHFKaPNz4IgdFq6bf60eYmpfP5s6XyNv1JNg9/dwwePbrYRnhqT
G7NjMk3o5fb8YtB5QxON0lvV4bqm64+xNM9qXv786Ks+Pzp5HCaEXXulT131bdmsaLfPJfQqt7Nx
I9MvH+lvWzM3vaPX0vfm8fiB0+PR4fEbRH6MZuq2tX9pr8O52YSVm7wFCMlgJqN8a9Bf6rOK8aod
HmqM/dmqH/ct5Z3r/Or46PkbrbKpb2fqdSarG3w53+hXL05cDjl8OjT6o2dDP/Yh++s3Xk41/7L+
42GrvoS+z0rrYqnPGOsJkTjzr6fLTrzozL5jL54eLTz7gjO/PTz/5vvDHxaduuTMKw/eeCencZXK
k8kNjz8MoDx6wL13rj8ZHPvqQDg8GSymfL6R9piRi0nj7wd35ytbb/pZ0ItLLfLHc4t8rqVhrf+U
xT5/h/NrfsHind2Pb/DCN5698/9+/uitftF+/kR/X39n2cgvuHytE6tT3wfR6wefn9z527fOhDhx
6H6h792X/tyg9OzJ/Qf3Hz1M+U7+OH366PPPzz+8W+5a/HttYWgYnMrP4FWu6FYu4VeuujbGu0+u
ci0/wbf8BOfyE7zLq7qXdi50WNzvzHFTN90k7PqrmGXfycH5g7O6UT02WLfK0F7nddIndz+cfNOv
ppuTL3uw0pvW8A7Sr/zxuy+PHh9Wfrye/PZGSPfHZyefn8NJHJ780N/ZXJd9nJf3eEZ1nz63DfSS
81edO4z04pPH576ysQ7vHPcMo7/waJfBZsInpydHMt4Z4/706Kt7JexajFatAp288dfe2z05fdZ7
3ze+/eaHHtnWnx73s3ZWf+g9TL/O+h++Onx2/PSH+lM/iS/8u7Pe0xwftR/vHjx+2qcPR8/9kz7d
Pjx+8sbkIn87BYL627o4FQgP/UxnDm7XI8G99NBxYv/1MtDc3qk4DZL7My8GyXrcaWqBawSXrfKR
7W18MQMCGPEVZ86UvIw8kk/RqjNnKmKGU5nUS526PXOqzOAyHd6ZORW7udRVd2dODUu7xAjvLTh1
ZJkL3eG4SmY4UxZ8mVuNpT3ljmatYmNpn9eHKzenMmMYm6svPuu9p21ja8VQLzh7xjy2r3b2rIXs
XKnnM0aye7Vrz9nJ3lXGfJmp3JozlQUnX91ahgy6dzw8W1LD5ORcdQOKfOHXPrlDLVD/6SUxsSKn
8HD/2wgBh8/rW72Sz5sNXbOObVYCnn7j6uie+wt7hjl7wvvzMvI7MxX89bM3Z6ryWqB79sA+e+Ni
SurnDDujcMrl3/z1g9GfGmSsHVgM+sKTrIVXGE5tOIizFxy6ZAUshTKz8+02vnjK5oH3Kvbo669P
TodZiPlbO/LPf6ZpnDj+7VFMjyn7v4UOK0IS1uvtr67o13k4h9cU1ESMXOv50dnRef3ZobCDXg45
O62ia63M8F00nx+93Q7tjb5+eKZDz9jR/Hn9qgLho/4cVYOcPX5+/O35AgvTGLy6oQWi7DPS19+J
FGiwvhXmN7DDF2VMbUf8RdltHKsnjY7Zs6S/rfsPykTAuvLqSLk1f6jOToc97lPjp71puDQ3+b5P
Meokf+vs9pq6yU772N8FHa3SwjSeH712Qede68+9tnqxTV14/NYA3oK6oqGbSxppJ2rR6tzhYD93
apm+nF+Kq1fJv94pZX6N/bWfiZ8Tn3uDv7q5+nfxAqyqz3oWtHShx3x5AdKLZzQdb5mtZ1cX1p8+
ftcDLyrf0bEvF56ia82ftrvytPmOxXl7i85zK6gHjFaeBubsfBiYm9dujo9VS9Ovzars/A0/b2hp
ft17T1et+7/X0h9GZsYBLOvxKgcw8gHLTvc3hPVjtuDyYzcw6wkWNXdzeVOzzkCn++EzfiBOeX16
7WJa746m4LYmOzZSGq73shv80EKfsjC+zxN2Y8Q0+ja2Nepen3zqcXFABhEVx0G1uwhHjC9TT3iD
3VpfT999PfnwgwjMDrKX19bNkpXvjF+k8ha7Qa3whl48NXnypS707mx7C3b1nzq/FtL1Z9d3V/0w
1caC3kotL/3I8TJyhwzeAZVKMWT/3ETvmcIZlcxP3rxQXX7WtHAvdAwxfGF7NyuN91W8ayV2ez0+
eXz6vD7cgpnNFFeyadd06Vl8OK4ZdCJUH0TtodT4QftYcR8I4qNi/2U3QS+HhzKWD97qcLbyIvUS
Xkhbb/NQj2LM/Dw8lREbJhy+t+CBrnc+P697kNXvlhRT9of6KfODP5exOphfMYrmNsbp47YwyHEo
2vniYoj/CWRxzPpXz2tFwjBFgSbmjGNMO13Y2ThtupW5B/WXrBMdPqyXC27pbtxLnX1C3NLb+D/n
h8dP65K56Lr/p54uG5hq5BIAbW58p1qmnzVPOevdjHvd+fEfP6D0Sy2GKcsf23t7udqGB9rD9qTN
5c7ajLOmB85PHi2w0Sv53rj9xius+nF5xrJz9WK5b1+cr42g9fD82MbnJ/RZ494KLd6b3jfwsub3
5PTo7OSNc8xwbSjOOW5wajTXDbt7ff6hP+M2emzvzFGS/7SmVzSOcJKfUEd65oC1oTay1ZIMmGjm
rHg6do0LLjr3orqhWlznRcjvqPjzHb/GO//Cy07/Ze3r/1r4NS/LrF9/fjL5v5Ovnx/5CyCZJj01
N798f4bLjS/2f7R/5Jy3m/7gncn274vyt7IoK8h8ePj87IgKr27uyeEP2kPUK2z3I47vs+/+8Jv9
zb/f28v56WmPaL87Pvp+XHw4+b9z+ERlp++sOGdySSK5NVGfPJrc8Lq9oZ1glYW4Vj0lfambfZQe
PPS71UKa6vHcbS445pJ3VYtxL76jVQXu/QzHgnl3fH54zRXn/uXuRxP7+O6jSf+ffXyHlxV/cvf+
nZv+68NPH+Q7jx7Vfz9IH3zy75NP0323vgcfTe7duffg03+nDm92mwR+zq3aVE/Z5JynUN161JW2
gt+ZutKZ92UrkVId2tzjNn26+VLgfFnJ0aUbGaqTD/oM8OD89OD0q6/eUMvLXst+0XYGPC+vX4aH
29+b3uyAXw6nv9ST4FO7iy3cF+GDS+6L4Lzx10fnkTrFrk5X2Clh0a4Ay3dJ+Al7JFx2WwC4cL+H
SO+eH333uLr4fpyG3var+9rb/cScnU411n98+PT7wx8qk/7V+dvPTs/OPZn46RsvTKb//Erjvrd0
3C/aieE3MODda69P/q2qE2zp56PwCls+rDM/3bCP98yzkeO3Hy94MpI6OW3aPXuq70Ze9614eGci
/3WD794ET/zh7n278+mjO9kmsc/b6nKZ6OYyhzhLhs64M/HDX1bgg1gprfLHyfX6kOCPb/iDdm+8
dC/38r3RKZsrT9kcThnTuv2FbuvdTn0Do83+aPUP/tjW2dHj80VM7rDV2+LejWlOLjfXn7r/2zQb
2q7oX7XXTY+4ypdD/4ajZ29+1PFu1IMb41P+6Pdf95ac/qw+WFpHZPh0emCWBJ32RuXXXqsD0ixH
g/H+gtGcGk45gBjTWgk6NTQ1OPXH9A313ufuh/VWbraz3mxNzrbZW1j/97ihtgjbSPcHTE+Db7ES
/X9D2NbJSk1u1WeZTb6M6awtDV0Z+OXhp5eXuHc3o9/IzXtfltz7yNx/zlvf/A3d+uaKW9+8+q03
AaD9rd0LW6TM4/rm9sxTXd/H/3X04KvNDz/uF8YNXx7spHrJNRRKne/rPfZc9aBhJ8xhPEaP0b9V
tfyNujMKndGoV72kVTFHCz9yharebH6x6Ew/9fVJbOvhv14UVBiJU0emvRU8fH7UG4Ku2F7rvvyx
/B/rqT57Q1dGsYMbPT/8WyCJhQfUytz6WPvIKV8Qyi5pugN2nrAFEt1dFHIW2fUwdu1B+MEqp6y7
zteiY156rcb6zen/Rmvr5SB6DVBhOrRNu/upkDY7euuL1+xUOJ2NHgvW8V/9FYx1RGZW8LLZWrI0
XxvN//uT1X5/aGmI+eMAQJyb7/+gGo6Hsr7G9vjkxVGPPfu5r5Vv4zGdvo9LzWOF71NGXm/qvUuc
tzG3Un046jvKR3vYRcfJx6/1Jv8uE/DuZFWXpF68Nlj746VYbXowh67rTTdTaFezN2w4sSgYz8z/
8n5urhy6FSduLThx4bjNKsZ1Q425Eax7Wsz2qx/AboHdzsjfbpLT+2vMDshUhBmMdfXQv3Z9waYd
F7S8YiQvHRL7G/L41bYMnw+QNXcZ4uPYoByZvNfndh89+HRy/8G/+fY9W2+v77y9uTn5tzuTzx7d
mTzK6ZP06eRf7/x73Ufz9W9Ov588Ozz5oU/7vvLHlIb9yX1r8ru+jfT/p5796MHHDz775MO+YfOW
KtWk1j5Ojz7uwqfMRW3v1oVhewr+LY3bftTCwL1of565QDz5l9kmInJfLrhvLQjuq9r8hcP+/G47
i4L+RTf9jwYHZrzJgiAy629+hwILocBcFrQICUxn/0uAwIJx/Ik4YC5kLEIBF1r+ivY2rtTeOPZF
y6uD/iJO57VLx/y5JHRxyF8UWF9plJY2uBAIXGmYrhbbp258ZWhfMsBLIvuqdq80Yq8U8wl+wTZe
j1/XXwkCbMLwXja0f/ighfRLh/PNS4XzzUuF882fHs43f4FwvqjN31Q4X3bTv4fz38P50nC+eZlw
vvkbD+dLLX9Fe6vC+Xx7rxbON39SON/8mcP5qlF6pXB+uWG6ejjfvHQ4XzDAK8L5snavNGKrwvnq
eL45Hc83Wzz3nT8n7aVRVTO9l+7en3zy4MFDLx1yH1KF3ZBO52JYyNejTWvhzpdG2vkmmhdyJbEf
iydeGbY0XvqYRwn4SGv+Y5z99u26aqu0NlKLb4++3fgiisN11bMXX/r7KQYvOk3MNzt7P2TrY7DL
1LuttLvK8RPF52DqqsL9+bnaqntQjerch7U56uvMnQxPbI36OXXE2zMFIc0aVLR8uHQUFrXcfn5r
fOrbU6fOX2T+ZWPdyPMsgeFS1324lh3z1vujqoHWuVh5y856f9YjcGDTSbV65Bsu6sMv2YWfB+lc
BGWWiYZjBz+fGKx619vs9C5wy3PTu+CYS0zvgrMWj+3mZaZ3VWu/QBd+9eldGr8XT28LjV5k+mGb
5mFHm41airi5NrmrbXQnL07+enL6/cnE382x1p6woslFlS4bkbIs3p18iQn0fvbBJx+O07Cl57WU
KNWHRGqpv7/+4sifT/7Kn8SrDl1PARwdTVpSCvqZfH/0htz5JZOpX3VaVysrS/P5cWH9zL7z8fKa
C14OM/syme5le6FK99riAtWHz0+/fHr0DPzN3pC1PjV2behef0gFfd0oqT5v8/T08V8v2IG5Hk6N
ajsjqlO71/s/vj1TfeBuQSXVojqqy2w6tLSGauXEX7l+6urVUxfVTi2onHq1uqnLVU0trZkau5nL
Vku1y/4cpVKXrIO6fBXUz7OkL66D+TmqnxZyfVepevp73OyrVTu90q1u/gZu9erVTRcyuP7UxjxP
e2H16pLa1foIyFQV1Hjd1tVxJYa3m2pqYUHVpRbjFYupRqVUwT9SV7WgrOqVi6qmyNiL0MNFROyr
FFNdQLYup1pXA5nLrI2VJOvqZXPZ8qlXLJ4KamSOXl1Nri6nVl+FWF1Fq64gVWcrjl+hUOqqZVIz
Y7aMSl1FpL5addRFtVHL63uaB79yWdMyRvQqJU2vWND06uVM0yYxTWz+MhVIP2f90cqwNVVvtEBq
XFlr1F1U8nM55HPVcp+fo9jnFyj1+SWC0c9S4vMPEKSWqoCvpAH+jwhQS8t3XqV452cJTz9T0c4v
VbLzCgU7l4pmq4t1fuZSnV+qUOfywe+nluj8wgU6U7ezJCaO63GmqnHmA+RFVTKXS4ivWiHzc9TH
/ALVMX//CLhE9f49Av5PjYALK15epd7l14+Ayys4fqkql1eocbl0BFxe3/IzV7f8UrUtV4uAP6Wq
5ReuaZm6nRURcHMcAaN+ZRB+qrGwU+b8zkJjRnP2TZ9//OO1Ow8+uvZeN6U33W+bDrI3VN1qxe/m
6MnUaz2O+0BZdch6mhxtvPhPO1CNXhP37ux7n1dfs13PF/bsZRY2f7ys7Q9raPuWXV7WJv/2/Pi8
bRa4trbWLdjn0QPdJIL5oui36CxxygTByq8OVOv04drUtO1dNz11y5ue3eZu+sjRlm91p9+pt2x2
XT/PXbx9+0AX9Dt6b/bTZiu+y4r2e5w6ZTLaema85WATn5aethAsvTP1MZz75MuvnsvwfHfJpz9M
HutNiEdP1ibp6dkpWwB9+d3x6Yuzyfff1F2Hzt94olfefvXiaS3gWvObqG3NdXlGSR/6ztGrezp9
cu3y3LgCbMbD8/bCRt/rFgz+cGp88PasmfCueWG8m3XjJY8VvWmuRnsrhObHo3dNYnd178B4J3g1
aTbUqnG2/vbs8Pmxdof0eCWH8+LZ+L3odeBdI1loYzPnxcgtO08W4me1d5WPLwgMmvt6mBEv/PMq
icnbG6vfr7Mc9b7xdp//MMELOrdSveDU9xevOO3EH80v/DP/ZoX3qF2cA7ACDcvuYvouh1F6a7Hx
C+97BH9H6DmkuUWDfXEz3u8Lur1oDOc7HrP/1uyyuXrnr9IUNzCcX2+0TgqztGBM3pmz9dja6HX7
+M6jO5MH9+88muT0Sf7sk2R3Ppx8fOfTO5P0af3ik3+vD4Z8OLl735/3/NYXLq/JO5t8dPeTO4/W
6hc6/v6df5v825275WObau+jzz6tx0w+fPBvtZlkfvCj9OehYbZYpsXxzUlkqbc3e3PD5oXzc/Le
So/zQepv+sO7f67PIV9qvKauHeN94Ss6795/5FeZu8yqnk9darj7Cy9Ww9JshJ9f7N7MywvevF6/
4nWwUe0zbJrof7Rz4ngXStpf04vGtSHsWqCBNYckerXtFd9Fe9WrRHfvP/j03s1X7G5dbJd/Ee+r
drdeJXp7tz56defDm6/SWwIw/f3lestVor91CS0d3bpCWjv+y0JDmPy00b3oKuOuYgxX6+rYCH6t
rq6yghWNjC3gl+tqswCHt+5C4r33Y6DrXxy0b6bL+OZOqyyA9+jjO6n3RkttKk5pfYsPfrohLW96
um9164qFRrS0gbXKgfykRXlB0+ON6kej/+Dk6Q93v3rkTaxXZibK8xYO/9Lxn57JVxvomTZmLjc1
pNOHvtrYLWqDd0u0qMc+vIzvaKdVBwVLPq/uY9F3yz6XN1/Y3KKvfOVMz80oa2Ko6pbW/iDLCAkM
g7j828m1tzeWnaNvXr7Kmx0YrwVvTziosGNIm3WF6dcpaLBn05iAJ2MOYFG6U5HT4rdYNK5pmLb7
yZt9hVt0S3nFO8TKliKwZbc3goYX3ODUzQ357EzaE9nncMBcYsEBs/T92enzc+VEoRFNK1V6Ymzq
yClB6a8D2V9lEh3vHfrTzIE/qjVxoEO9olJPvtWWtH7x/i5Dd4tvt9q3ww6pw7cbXwxPUjnpjuuD
heqzuylNaJr7jQty3Pvjq08dN2Xa7Zipd16NnEm0OPU9VjN/ns4ZHb90ut96f+jtKlrAD+SgTtLF
vAN6bcaVXP9rPHUy70z67w6GAR5+3vhC50Sd7OvDn8mnR98eHZ47vwO8mLw++6dbzKtszlFdA7e1
jKjZnKpEHa+KqWaXtnxwMntOa3h8ztxYz/QoMsjNmc+HdG9TMWFB94ZovqBB1vLgixYfsnA057Pe
zRYZFt3won6MbmB5T6YPWjJPi1LjzYtS15bX3/BXfR+fTA6/7E3wsD4d6ozFm4sz/oU3vzDn37xo
a/bI+Zf3YBkbsOSel/ABmwNuuGxQm4vBQiLz0XVtcm2MGO4uODDiVBz68vcA8rMFkBjudsq/LLLE
YZbmp2nmnPGKa2dVV3xVT9yN3lTBjzUwxc8Eg+nfamiITyL4jX+Pn9X59lX82p39cFZfJ3ft6cnk
7a/O5rLL5bnrKOVY3YYnCZdoKJKJq/RohlP5KX1a0lRLcV6fPKw7ses9Oagl56f+nrrnPUqtT8tL
tTxb5cA+PD47f3785QuXY0KGhll8fhZ84GqtZ1XJzPi7Dx7c7ZuImuY+DaTESz9stjqhqfZQO0fl
JTPae2119pG66j16z+pPn2iZ3piqhqnP7n19/N2RnTa548aoeof34T34SoUos8URy3oQIn1T9sc3
4vNzidsYNTL7iOCowGJ0eG14XGTDsXWs/XWKDPZUoYaGRqUosxXy3p485/TxGyuO3xodz8hGz+Yf
oJ0pqFjQn6GU4zK9mT36Kn2ZfXoypr74eYsMgBZphrFedM774yckmzHO2uFrYYXvX82wlluW2vu5
zWH9iuaw8cVPnvPVl5w9eu6CCx5P/HqYn58yXNP3OzT68i/bC0dkfER/u6C+uPl/mZ11an/mp74N
xgJfPH5nan18d8XL0vkzlDGO3wY709nxrjsXbzu/6q53LhyXjalx2bjCuMwJrdNmN/wwszvR372/
o+aWFjrOLpqx3S82vMXf/2Ob3eJ73rlgTP4+Jvd36+tic1seDMd9/4050c1/aGt+dae09d/Mif4M
/f01neg/ttm9mmP6+5jc362vV3Sirb661gI+eHjnfp+CJ3u0+l3p0yn4UGn77JlzfBcWBlHPM9tM
1CGPXonpfWlirJpfKsBOvw4zDtebT8/j1afRr77lyTXnAT4/f3vDK50nzw+/H35xZaP9BgH1+XnN
kQ6fPq2HUh09fOQk18xn4q1Gn2qXqZnz9elsE/p4aSuHT/7jxVktSV3cXPt6Sbvt++kLzA1TPwhi
OWKE9GsbI/2+aJRCZZwep7lP1YWpz0djtejz+Yamx2vRV1MjtvKApa3Pjlo7hlLhk1bBP1389vl5
MNTff3PcO7vvj57XavzjWk99dPIkCrwl3tb1+Pzo7edHZ0fnE+iqyp3VkeK3uZLdeW8smstz+jce
n35deay/nta/+5NePD2qP9VXNX9/+MMbk2HnT5Zn76RW8WZwyH84h3748sXx0yfV6dy4fhy7oHxc
7+/7fpaO2Ibt8enJd9VG/P3slQOcDIFlKKKuR9Y7OT6vdzJ8XpvUBQnWdUvx131xn5x+X69z/vwH
XvF+3nfaH4s4Ofqeweqn4KvT9jbww2dHvgXcTT+8fnT6/Phr31DMu1VfyM0rJc9O69fPj272N3he
r+F72PTT8uToq+OTI07WIB8+f374Q+09j308fvqin5fa3qjn1etvNnFzdSRdVZ6ukRyNX9N06gZ3
fSf62Ro+CoF0CiTXOx2RfQMIou51lqycf6SwHtYbUN/o46M3gi76O+/xMn+Vn/sRxcArTVRan6yt
TXYu0/m/XP+rx/8/LRruH308ddR7q/voDU1GTfkBKxqaCfivzS9bn9JrPpfX4uVwr4fJOvs8st/4
jjvQd7Eq67dL34a84At1+72Lbb4+4OgPLOET/qv+PG/nnfDOz+6A5lft1EMS0xHz+nEfRidrS5TM
JV+9M/dEw+jo+eL6f7lKE6p2mL2FZw/ss+V3MF1UdHB8clCPP6hn/5Qvl3Vm7AWnEo1xBvIKoz3f
meFKS/s7d0jrdXWHePuRcLNxc7J5c7J1c7J9s3oDaax/qr145BrPcr9SjwMaPBoR8LPHdYPb+Y/e
rXoHFqhRU1S7/xKeaanznZFnFvgvjljtG1/OX2SGAZ+7xOz3l7nA9TZSf7n+Hyvo/f5LTruuKfDD
J8vI+Ti8vQX09fCRqG5jc2st9gn2sD7jw43RhzemPl6kmfdjMt/E5qJ2txa3u3VxuzemWn5r+G37
iyWHbU0dtjNz2OoDV/ZmEfZZ7kWXVXdEdxcXKc76v+lbXPQYz79ctakrO9P/Jo506YD/zB506prD
gp5aT8PHMytq6ovFNSVTq2o4fnNx+1vL2t+6TPs3Zq7w1vj3uRU21fjUoQtW2QUHX9C3Obe1bCUt
QWtRltOtxGWeLM4hsUbjwON8UGv4r91WoUhdUjzH0efga9ePx0+0jP4Ey3P9n0cNndxc2JA/GLOs
pYUNvVjc0OnJ0x/WyPPne7e4ob5LFzQ017tFDd2dHyPV16wYpYUNzY+RGloxSgsbmh8jNbRilBY3
NDdGCxqa6d3CwT5bPGt62Ojyt/bB2YsLZ22u0YW3Nt8j3dqKPi1uaK5HC8ZoptGphlpLW8u69OzZ
2pPjs/O1w2enJ19fZf62ls3fwhYvMZF3t86u0uKq244Znb9rn9HL3PMSG5m/6+UtXmqJfzB/18tb
XG1+zSl/cDb93M7w+YvFX9xdcsLdZSd8sOTjkyWfL2vmxZIT7i75eNnhy+5rtv3hm60lZ2wtu8TW
2ZKmPljS1AfLmvpgQVOXeSoOg7H0v+MZuN5spusxrbee3jiGR+FWPpZGb85PVz2aNtviiJJpHfFc
W3TiShw/f2I886VMbubsBcj0qk0skiPnmph+HOL0+egpr4XHT4OmS5xQK+UrFvZxdyh8ePJk+O0K
J64+dPqwqdz1EgnWB8sedXtt5Fh+Guc03drJ1U4cnbl1cVe3ft6+br1yZ+8uf8DutZED/mlZ6FyH
755c9dTRuVsrH6J4bewsX7l/Wz/XLY9H8UV71HBks4s+XPjZyfyHdxcceHfqwEWra2lKfbkltvxJ
l3GgVQ+uuCoW3KMvlaXXvLohX8p0Ft3AVexufla2Vl/6t2Amry98vu4AdHfQNORFj3S9N3f6sGZW
NjAcNt/ETA9e9cJLr1dPndyvUsvw6g8OPBveIOmw4/yHb6soumy7oBG7+s6YP31v6pyRF3l/hkN6
Z5ojEsW+fCRGf96fokrfGROdc83MD824mSlO6J1pLueiHg2Ti2wwJmzfmhqgN/umb0wRsG9NUc5v
XtzrsSXdmOGu3poZWV1uiod6a4aRm6LLFxBG149n9ztaNuErtji6+ALx1OO7q41j5d5GV7gNWcC7
y+1n4d1oDq56U8PVVpnZknu74jVnb7GZy7uTGdVggXX2DrnvVT2u/llpqX7sijFqdvqKgzXT70vY
+tD5C21+rvPLl9l0Frg0eZnGEaMHv2ci1PSX8Xz3/EHDtCxoYKSHTWUgLxccPp7gtyerz10CLKZy
sJlAWuv8VmGO9lDxopu5y82M1MALbufuzKS/PasmXuaWnAlZPkMni6fo4GTBccOCWNTI8qFecPx4
ca2ep42F83RywUTNfb/sxvzI8XJZ1NCqQV9wxvTyu2jWFt6gE2fLZ+3sUgvr7GKBdx7Pn100G5dv
56KLXzjxKy81BfVXW8PZRev27HLC3Tz+P7t4cq/S1mU6cQl7uuCS02nFCjN7MasNz38/QhQLvlwx
vQuOH6OT1baxudApvLjADC64HT9/DFkWfL1y5BecMQ2BLpq6hbcFNb1ikk4umqWx796aG/YLnPei
E8bAcfVEbS2eqIvc90X3RAtjOLno+5VjveiUaYB60WwtvDeUjRWe+6I1dfbiwvKZyOgWON0XF83G
1du7ZGcutIPLXHmGGlnt1C9azmcvZiZ0sVtcevd3NZorLeBV2rxCpy5hg5frwQw5t9xAty4FLbaW
JVYjFzZ19ArbuFwDo8utNrW55l4llGxdYHxbFyGKrVUJ3KJ7vLt1wVTTyGXjx/KeXGRSCzr7anFr
a3Xcmv16aeqxdbIsM99aaCsnK+d7J4bxEiFraR9Wm9RcJ18lRm5dFCPnD1ia42ydrKIIFt3wXY3i
8jnfudgclzY81ZeLjGtBd18tMm+tzqlmv17m+c6WM0Vby0pM59jMBRY2V5f8Cob7E68+uqtLtXPh
MrhUK2/6G64viTMue/HLtvLmwrV3UUK5dWFGyYCuXHTLaiQvnq67Wwsq5V9pff7kPkzd4SXbusSK
v2RLMp3LA6rLd+LyLY1MaL5OdPzn9V8xv399Ul+pVh8UfH70ny/qY5zDS6G+P5o8OT1543zy/WE/
keffHD2jAqc+iPT98eOjtdrAk6OnR+dHPLIx04WmIi54tGO0N2F7vsP789rlHr+ox73SExj1xJcL
rzQayuXXmj/oElerP12qXkprZlyoVC8SK+rlLKG+aDenCzjXto/TBTRM2xHsIh6gbQa2miBs1700
MfhiYQOXddhLSPvFm25dSOYuGbUFXM/0uK2mGqaHbhVduHj0VtOEc1np4hG8hDNflX20+50PI/pq
+3K4cejdsoZ2Fja0dPle0NSVQ9ryRHV2DGa/WDQCc/nHbKdnv1h09wvQ4CWauTL6c7/DBqLNNb3X
Pvpg9OPJ6OcX459HX9wd/Tj++MX453FDZ+Ofx0edjX9+8d7Qy7tbo2+2xhfZGje2Nf55fMV6VDc8
YczTc+OdRc++uqDudDKuj56uR40dCriWbxjS6cL+jMvlXiT07eGTJ3X3hvF+Iwu2G+Gw0d7Fl914
ZO7EmS1IooP+ZNin9S2Pn59/4OnS+b3Dv31+nr77+vPzf3129Pzz87tCvuf3jk8mf9Dzyuc/8Du7
9H5+/uhBejh5dvi3ybPjM17EB1Se78hB3dqgipjvT964n95YfgiFEYtf0rfklB0/ZfZFbgNMWHLa
7hfxFpal379fS6uu1bm9sLW9C1rbu1Jrty5o7da4tfGs/sfp8ckNf7zq5uRPi08eYdopA77YfKfs
ltdSzL/UMSrixm+inJ+fYX/vNmEvFy2FqVdoXnYZTJ20aglMp4HnM0nf+ewrNNsXY6LhfKoOcu4c
fTx6f0W/kJ7/oBdY+EhOdWkmNT2fzUTP++Hs/5rqwUJm43zxG0D7zxe/S+NqBnHv9Lv6k++gjD2c
vfiyT0eqv6yl6nWWfuQtkHWLk7pP7DfHX9X3i9ePnhydnY8+4rQbfuRNvn2z991M9o3eUdWvJ18d
Htc9b/7yhh/3xuTt25M3/Ng3vnh3cv2fr72p15beP/3+Zk2Dvnr64uwb37am7pVy98Hkyxe1OvCs
T5VqbjPa12GtP6vHj4cnk/sfPZqcHZ3Xxwq9jaO/fXv0/Pjo5LHvuHLoCZXusrbztO7/4sZ3Wndr
+f74jG7Wl6cOL+nVCwKfHGOjN344On+zb7k+d1PfH1q3efFNXyb5Ybo/YYOc+iLQerf1Vbrvvvvo
h5PHN7V5T9+l07qAzmsKONTu9L55cnjWt3fYPnzP9445Pq+N1wsf1gV88uSw7nx9dHj+ol72ZPLw
6PnTNe3AfePaWX+lOpLd1FTeeMX3273Ki/Hq1X/CxduLyl7hLXc/8dKjF4+90kvrLrp8Xb8Xvknu
qq+f+wkXXTrUF5zzEy65YogvPKtedsVu7LM9WvHEcL304q/rNS5safzE3qKm9P2l2lryPOyiZucP
vfIVLur5gmMvuoaWwvK+Tx1wydZW9HP6iEu2d7lRXnrwK1zl4jt4hbH2OZp+unPRLLYjrmwfF7e9
8OhLjs/K1mePeYUxv0z7P6X/Fz0QvPiKK876Kde92MBWnfZTrnyZYV594qXs8oqjfdE5r37Ni5zm
Txnn5Q1cvBQvHuPVlMaHpydHQH/lf3o5YIWstbB80csCr5ESLsjh7t5/dOdT83deT/4tPZrcu/vo
0d37ZW3y9sbko09TtrsP7k/uP7DRq7RvTj74zCYP/G3af06ffFbf2/3g00/vZFtrR/VnPZp8nP58
p2/5/t1HH9/5cO2a3sSl91PMP8v1x8hcF3Tz/mf3Pugv9+Cjyb308KG/srv2+9Hk7sI9bNYmj06f
HfVQvR+yHoV//7wO9RpdqC8cmnz02X2/tUedJ1C+y2EkTsftkZ4DXvigrQ0Pjp/w4UZ86Lvn1SOv
VeqPGzuo1KZedMr37/Phe/7KdO2sqAb6qWwnLn56t134/zd59s7/e+PzJ2+9+fna2r+8U0m2enpt
3d8fMtp/uW2+fOrbAnLhNvIzmxfV3YLopz+UoL0dubbapwuvdpErbCVdrwNVNVzxZV2Gfcp0olF7
2ebr4fOjYcra7ERq+6pDqov5kI7ux//MjFx/9LKdlma/mrrX7rUfRxca32nbVdoXwIdsLFoT0vpm
em74jT4bPTqpyerR+96+dt15GSPTNnictuZpYw67nDLbjUVmuzGY7RTB6ye9wjz/r/8VV+9Na9ra
xjnCf7QN7aZf28ZW24c/Xv+P9t6cl9NDN7cd6ejtRljssM3p7BWnzFqbg4bYOHfdYSPtmEsfoGai
Z8f/dfTgq80PPz48+0Zz8Qe9ze1PB5qHegxM3lRf+h63fn8Tk+LHvjX/4qdvfqzHtxUYdlUP974c
/e34vF6h+6cL/pw9f/zOvQc52Ufir9a+fXrROVf9s97/2dvb6f/d2Vzf2vXfd/a2/d/+z9bu3t4/
bWxu7m7tbm1t1eM2tnZ2d/5psv5zd2TRnxdn5/3ATv7p8enR029Olx/35enzv/4a/fmV/7z+z++8
OHv+zpfHJ+8cnXw3+fbo+dOuU8nD47qa+3+/PXx+Xve5ruya28rksF/YP5zVr46/PeohylF/StbR
Nx6/2f/87Q/Pq3o2uXPvg098/+UPyt2bk831jc23+7+2Zi/x/Ojp0WF9b+WLkydHzyfl/meT8vCT
yXdba96bw78enVXGcHil5fnhl389mhzWbZDrPif1Ct8+P33you6MHFTsGc8h9weBxmpTn50dfn30
7mTW5CdvH5/06MIB0dv90fVFg/7LX95+fPoUBnJj8nbTrvqfz46+nbzx+fkbX1QW9IHdqXvX6sja
m3bo4fOj0V7h3xwd1hvsv+Ww4Zt2gn93o/3qO1R/f/hDvejdydmLb3vnePRmbGpd2UZvou/F0d98
N+ozNsVuDT8+ffasPqVd5+msspLPfSO374/Pv5m8vsY8PDs6PDm7OTn67qhunn364utvNNP9OFc6
9Zt+ELenW7o53O9xbfa4j1gb/XS9cMnu+fHjPiJ7jczh85P+emf8Vo7OT789f/fdT3ps9p6/bvvG
9eOTmy4G3pz86ezFszdRO6bfFM9rNZ/OfFKn4P3YcLJv+cG3Pu03utfeOD559+yN3u3ennxe2+8/
6S9QP/JP6sX6j6LF9nl8UL+MDrQv44P6ZX9lb99de/2y/6B7rVK56Cp37/dJwB/7K1+wX0/v/p/U
qqtrcYsvTs6PpPN8/0213ht/vHv/9pt9MHj8zemzb3t3H3JWj2Jef0dhQue99VaNUf5C5Wjp9jBu
CmoVi9YJrLjp26d9mLhWO98nLdcP3lN8DMZumIS1tQhBnPt2b4xvEiOjP9frF3+5fvyFA4D76Y22
Dez1s7ot6zG7uMZRfq2Zl1O8jFdy3n/v1cfRz+uXb2zBdcGZ3z8/rsVZ9cD3uktNgGOpuo38gYax
AarxlEwfcf3A8dpE0LJd54+z0zM9lctbOVpw9AUzO5pY38x8mN1+NgNozVxw+LVNnLcpDHZVY/F7
X2Yny69dn82u+mv9p11+PA4XnTy66DuTMMc3Z5tr0M67OW7Qc4bXcYkLJ3j8W31G/EaLQW8C1F+b
ae/snXrl6+/U/IO0nxewD0c5xB+viNEbcBu++8nxf8B/x097Qzp4fvTs9Lujuq0revbPAQcvwH87
W3u7M/hvd3d973f892v8WYT/rhTBf1tYseKJqGK+6Unns8OKK558ebNmhodPbg7JZf/J2fnRt/0/
x0/OJoIdvNvF3T9JWQ8x+HWGnbtg0bxbdf1PrZJqsHbT8ISONWihfvZfPPny3bMfb74UdKn9rmjj
POAGWKP2umKcJ2cCIcCcJ2fdm+IH//lGvLpEbb8ZlE//s2fhj9K9NyIxb589SA/dHTvxoS9q6csD
m3x456O79+98eG0S8X72Au/5WX/tYeCzoxvXPuJkV9qfVR38/PT5D5NTHwGf42f9wpt8eTSpbFNv
L313apz2Gqcb3/doM5pdEzcKYVFv3ePG2ePnR0cnYhGACh8+8pDfD8O1lo/TnWBFjp4/n+FEvj72
6f366Pzbr59/6yTH9DRPY4baAtRLbXjy9v5Nb0K1EJ+c9iCkvnFmSvn0HORGP5NDVx8+eOSwpn61
9vj09PkTBynq7I1rU6BmdFDtXoCSvo3bRKLAJfPx/53Pz956571WSz96b72Hw/UvKq/xly/alHqZ
xNmi49jtp5Z4/OnHBV+/vClgt/FF/CTq0wfmod4M3g/q3fsNVvlv3APLzlGM380NH7mIrQdeKj6P
brzu6nrQKo1TY9Cfn517DY4+aK+6918Op798clzfyFSLOLpgxmT81dR8pcSqqH/qETfGCOfNHsdt
rI8PGS+FT2qn5VienB5RHnJ2dPSs7tP5pZeFVOvnitVbvDlMWcVYMcrAkBtrb735+Ts3/rKx+cWb
1yGE3U4afa2j/RscYXzjDNZws/7Z5nsjbnTcAi+6Hg1k+4aNrTSiOke/cMC+DohRngCQWkNvxXRU
ZDhDIY/HvXqjIML+LmOwt3QMbi0YA77ZufLN+2RXM/yyB/NPano7WmH06iUWMfhcDn2z9xrNVMLc
4kaqE62W9lU9dEIdqRMiMrF22RPiW/1VS7P/6I8Nwv+oi70c27evotPvv/Sd1Pzrt2//5frJF+9N
HfHi29kD3toYHdJfps8IN4cP9C4fyr76+VSVa3uB10mtEiPYTr487YfQq8T80HrEce+7+lOOphbq
5J9vjBZ8Tbi83x4Nh1m5TWffnMwuYkcA1TGNP1wRiyYzf1SD/+GjtsvCO2xue4kDN+cOfNnNn9K7
UVKHljTMH61EaXzy65c7++VFsVf5yYeLirJt5kV5dabO3tXAukOcKslccdLkr0c95mQ2FpxzIRhr
eiRg7GdLn/7b/2n538N+iRzZ4d9OT06f/ZBVzFXV5TOvW33eryteSQIjeJW8cHX+t7G+vbU5nf9t
bm39zv//On/m87tx5teDzLv/+79H2nf9+yekbI6d+fHZ6fkLfvq2D4rHf1Md/YunTw+0PzqZXf30
4ad3/pw+uXM/3zm41zuLe5/dm3hsms7bvn/SaN+J0q7vn9Rc7Nnht6Nv/AvvSP3OV8w0XVx7Vr+i
X/Gdf8VH+vK7w6e1HrolivO9vOlu+o16U29zU0F7j+7zZksMW1hPn5Y/92FvW48NLLz96/Wgv2xX
RHNj4Qj5Ywfda8NLCx62Lk+eHZ8cP3vxbHL2jdds90D3696VV9Kzj9onk/86en46ufF1D1AWNPzm
kMHVPdCfn34/aP5nC7TT7/oE66C/Eb6vNPN3701Gumj/aWuwSh21xYPjrw76vvZI6ceWIkya0upq
LXGvfv7Hhbdfr6hLrK/5Bd3q/HHX6PefDnw71drKUKFQ84W1BadHf595bx9/0+NUDNtZ7RuTdL9P
xW9O3vjjGzdl6P3wx+i3vNFzTD1kXDGuZ6BvcPwbUb1fjf4Pw0ED2+3XIC0bMd5/Oj/969HJ2WSU
W56/o/SpGvOoECOO7b+9Plzgxzis5pyft0NEcN7wq8LB/CkY8Yg3fzvlkWkf12sfnj47PD6pRVYP
v/nh6Ytn9afcA7yz+sOD571vqD98dPjs+OkP9adydPLCv3v07dHj46Mzjq8P8PmB7WHstXq1tcd8
40dNf+Av52Ma7vVf9N9//eLkv46/nbz9mCH4v9cWTYZ/5VMwHvk/UUd95jf8Bx+Z+rTQMAv1EkxC
NRnpCpOTo7+dVzv5yRMTGX/0Qy8S1m+1WqEd+lq9KFlaXOv9+srjvjtfPT09fX7j+sGbvfGOzsDP
YP839PGbkci1e11iEHAGGEUdhDpgM4Xkk1s3J/u0x7r97ujpgT+18v4gp8wajzSnOtyc4NKJPqlc
/XdHtT9nM54lhjAqM0bLihtyFun56YuqRC409xC7WPocWgv2qp4xmtNxbjF23XMJyaivP9KcqlPm
Yb9fshphEDuLBj/64A6ISrXebanca8HxL9+bKbg5jrfjNk6kNuCa0iA56sGBOkgysiYCXh/Nx4/X
uKU1HXVtmqGK+1h5xnvDZgJ9R1z0aZsHvDj7hoqf1desfFY7F6Vy+bhPVQV5yVTNYGRZHmH0s7aG
mLEwDczMUiQ9+lHfvlTMewmBiAc6yA8+u28H9+6k+8PDrA4e1kZrYq1HH9fGQXo6TIzeExP+SQzo
zCXiXXOLv5x6EjK8yKq2nESdelG1L4ypbS9GIz5abM8O/ybG4rVLjKOf0k/Ci94jUMZ4kfGMF1u1
osa5DK0sSvTrwuEIrZ3RRV/Op/txG7FkdOoURNNnt3Xwm7NcQLW5Kbfhhy2GKv/rf03+ecatjFb9
kllicK799x5pt5ke81V8PQ/+bqipN5cxMguWgTe3gMKZXbMVh9ZDFxzp117o7Bac8XIRf7Oyh+sr
eaPrckSzbmyFJUS96WtDYB59768tHjmmKCG/0DeFSv7TnNTU5Rb4qenvL+OqZlu8jLeKjX9YDcsb
u8KimjXj8/a0fjOy9+YO+q4etMS+Vqy58+oihhXo0uB5TfCWLMHr9aEJ3kvy3TvXz1cun7kJ0smz
VnoJO59rat7U5wx54WxO2/LUIdWcX/4P4QYruHvnF77GBfUf9U/j/7bXN/5pfWN9b3fvnyY7v3C/
/M//cP7P598T3e821iDMf/Zr+PzvLJr/7d2d9c2dnc1x/c9erf/e2tn4nf/9Nf7c3ljfvnVra3vN
8vbGwdbm/uZkY31jsrGzu9GVXHJOOVnq/+Ri/f/6X/sfcknWf9P/aql+5L/2x/b/t3pCTqnjkPp3
Pa/Uk/r/9Yf7j7W5/vT+X6sH9W1ZSbXt/u9iXW3VL1lb7f/mWz8r1Xat/7//5efkepH6Sb2UdX2T
/fVy7V09LRX/shQ/srZSO9U3Wu/Lu9B32MwPyF3fgXqSH1zvpjZQv0v1JupX3qN6F9lvwI+qZ/f/
7+pHpfbR6k32N+ZjUq9Q77IeUn/0W6rXYFBro/1HXe1YvUw916/CAPVX96GIbmXzhuolvX2/nU6D
Ww8pDJv3wAelcKjp/nOdkeS3Vnxu6sl10nwiNRhMqY9k4q78PP8kMc50v5+qOi7ml0/1cG7cL9jm
xQ/1XvkU5oQhdT529cx6kz5Efgf1StxOHZ46YDTEMCc633Ex04T4bGW6kmRq3kmTWfqH9Zja/85N
zxhTP5B26me13xpiNaXZNP++73bWMnADlO17b+tBfgUfH//PFwt35NfSwnAbwJLcCrCE5JcodKIe
7h+kOKbec8LE+pXn4+GD4Oe7JbixMmQ+HxyTfcY7NygWji8D8xWry9ZPfDUmzMqnvq7kenv9qjLs
v6S4H78gBl20Xnzm/DTZidtx7Xbtkn/tHTRf3lhOwng0N94gxq+1Uc2TeapL0Ax7zHIh9ZpaBjLZ
2hU/vB6ZOp2iwdKyktHVgaoH9p3wVWR+sk+Td7/DefhQFB9dt4qC76s/JHc+xb2At1I/r02V/uQi
q/c+Zhyhu4PqT91paHn6kJkPPv+VTiuoHusz5O7MO+JLobR7qcPO6qu9cxOt5lnHJGVft5oo7Kmt
AZ/75M7IP6yj4FNl7tmyG6/7UbdTbCzLb2b8kh+E68D/dH5iwqVnC8tMRT7YnYGbUGbhytL9Qp3f
abhkRkchwF0DXqI6j+72xs6t7f39zZ29t8/OT0+fbvS5/3fHZ3X3rceHX311+vTJ1vqtzVsHGy5d
7+/s7W0oqm3vdu4/fN1Vy6CDuE5cU0pyUgVbiwhXD+68R5iYBzcWca4en6Wl8WUISjgYN9+O9ZRL
eGXGWRHWcPhMqdEDWTpRTW7Y3T998O+SOxzMkECU3Fq9y0aA6twwFT7Ct9TvzB2uGVceRTV64v92
BMfaxRbVkq91ppTlS0cZx5QV1frA5Fd1cy+MpPnVEiZaG6FbiqKseb8vd1jut33czHDOhAEABdHS
54jPTe4dt+EzVIdRF663l4YmuFvNmNuA3w7opZiM0/uIwymKAHJ+xtWJ+B6nGZUOj+E35i6ZqFCw
Zx/uHAGGz33t+qIsnRtSjG/KsVQSk6U4o4HwC/CDh73Ox9ZikHH9Dp+0hvmpZBlh0ldup10YEEEL
F+3og5DqKwDniT/wlexLp3g8jWsWoUXAmuM8azfvw1yENoyl1eEKcp4NaYnx83YV0kwrxO09ezDG
qjGpFCFtdEHsgIWrOMufamE+zx7mCvdcTB4+Aoy8P9/6+lXw6by9GonNIqT5vMksWYDevSyMUYRn
UjUSD4F+poWFEQwAv9UJeFBgeqursER3O2xdWNEv31YxmCEBMcJzyVP7SHbECkIV/+M/Jlxwnzb9
dP/QJ7q/Msf4p1mLWKdpoQQE8VvzQeO6/b8eyROzUTBRh1apDYB5O0KSrPkiYNPiKS7ayBhYXIkb
SiDLzOdurrLhTj/qrv1HdzHyfj5C3rrJRWDofq0OJ9qsqQir4pfcdPxvny71Ta0Wh7j+n/sxx2Gg
/L5Tt3f2d7c2Ntc+vvfw0zsfre/s7h/0ifIthaqN/Q4HmbkHMGRJhE/lBho2MIF8gaOZ3sA8Dnl3
CoANK/VPI2aYEYAUqkjABI/DS0yFKm9NgSkrCQQFKsmqnpdQlcBT3h+ux6LIkYARKpQ/slJ8XYwT
sBIJmJGA5VECpkRPaNWRT0nALYUqoI2gOvFFyQIJmd8hWFAJmF9FoTKBab3rajAJwODOQU/FYVOO
IW2xBRzgYwECwiMVeU7FLyHcLAiScgmgVy3ME3ANmHKaLIMSHOg0SQQ5IX5lrKYWgaTunZSiFk86
a5xzqJJJBXzUA95mZVDcgrJBppD+dEq+MMDI6VrKFImBwmtpN+W30Cn8JQBw0lARiWkgUr9Agpoq
q9FGTj6p80WpgRVlZUBtMoMSoYW10in0FKVVYT9ujzIY3CfAnqH3nzOJNqE/kr4kSJYj4Wc6mIKp
hKtrGd8QnjwdcIzHd+QQUA8l0pFqVF1ulshxHJRMp4ldUMal3ChD1RQi7Cg8kT+QwClrY50webIc
UFanjCsP4Qn0oCwAd05+QZhQl2po7HLMetY1SwyMbI4VJCMqaj+51XfkKYIxOLBSAmQWYSPxPlkL
JsmPdfQfhDf9dQqYLbvHvHxxJ4ymk7sCiWmZu1PQ3wrrU+cl3VCXU6zHMroDFqGoFWA6Vgh7oYS9
Uy5bYHNY0yYPLm/j3wVk5cI+KsFMhMXLHFJL9UsQD4U2YTA8Eeyzqnsfr6/vbs7lUttbOxsHG5ue
S23urO9sEaC2b235IoaO4nImkgk0lcjsFFABm+RX1fF24W8cK1tStM1CfgXPj01kNWiiJBwykjAQ
V0uDdFmWl3Nz0EFHJLpQT06kTuEoiBtQdXCHWlkRLZScOMkH8jDhcB/IFEaWSLRIVnwoSg7atFSv
59HA+5RM6RBJHEsvCY/5YhaOKKDA0iXdfAamAQyyCKEBxGo9yl+KfOuE6OVfLYi0Il6Pc1NSQJDx
g2Bzl8h5AH4y+MRqJCvHh0V2r0iDfXcpRyLiJ4tScv4zAYjdUFO0Tv/9H+t0Rxk6yDmjEtiBuEpE
CXROIujpTurkHX3MCoklDIv4DLI7pooO0ApglYOAokpNk6ab1a4emPgzkJ3nvoE3w/S8aW6eaGOK
oz4NSouKkG9XlFWS0WgJJ/lUQQglAwK98hUF5CbiJsU05vDyRdQpixDKoyhtwbbFVQMMLSxVuCCJ
8nK2I7XRgbIhEfNFxcXdtBRpWUSQlMKELLqieNdp6FqOWSJBT2Lk8teHjw/PH3/9WEf63JOvdroT
vw7pdUunAnX7SJLKEeUSxIF1IsborTgxMig4HofLLbjirItiX0dogxeFPIBAJbUQPyr2CjBrgG+S
z5SDri6QESn8nTOYpN/MeYbLECHs61n+T0QxUCcJG+BlSywpjsO9OGlvcF5hZzZyR+EbfL0qgdaK
JXMVr2QNisrJK4YDy+mUadTxl5W0d6djhA0jza051u2Nnf3N7b29vQ0ovfkw1AedvT1Relt767e2
fw9D/9BhqPwmwlBpYSj9HoZ+/TCUF4WhPBuG0q8YhuROYpB/tjBky8KQKLN/6DBkv50w1P9Z3721
Vr3IB5892vrso52Drf31jXVFm83KyuV2ncxyYbWKw8cyfQpYMiwa59s7LDBFl6NwQnQaqStSziwr
5yeLwRRnYRossXLEAdwoLIYFK2cqESjKwiGJYW1yrHGijeYZYoaohIDkeXMmb04SkBSAJNqUPCJt
CbOWXUAibqYQQSLz19Q1Sg9KN1gZN21IClalVggxOAUrh2CdhiUjtSa7aBZ0r0UA9ZAsoR3+0c2z
CUgSMkJAwp9rjSS4Ecg4RVbiShFbKiGiw9NDyJPq6zqplTLIl8j2i6UQ/LtYg4WyCHdNqVW04BSo
rghetlg41a5EOYAmPkyFBR5MKIUz3ECSCmKAYVmIhSOAHx8kDTm1nDXaYKDiJwsPSQ/zOYT0jFoJ
ZGUcJk6afnUqC6L5XIKkQxA0C/4ZbhhNleNLcaEe5l4f5SwNiQmUHKB4XHIrZqmfdCIxTWOr60f3
DPiXQmtoJIrff2cgM0EhfALBkAEPks79syK8VlfqBpYuTWtIKTQkKnWizoKpQIlFTJGeJquNaCMc
YjkcttwNITUj1KtSRWZvYuDkyLRIHPqBRAULqpOU1hbuJClwIPaZxlZVATF0gICqrOIMQomWhVhU
V6F1JDkzGH9Js7mqT/K9Jmgs5YXD3fAc+ChY4+TxoFVNlnmZBoZFTTYQOCNmWVpSUWSJLEGhMYU5
Y4bd7f393d2tvbVH5d7B7t6mUpM+bnR4WYxG4Y/IUJQTBAhVPNfoeYToFBnDnvEgcMWFmWGUpyUc
XAFU4IjgziEDLpRwHHjwiXwfHwkBiosGETk3W0wBC/ApiFmmJRwn7z0Ii6OPWZa/Dy0f3++RhnPk
JrLwqCgg71tp2VAZSThlSsIhIYtZxsfoCqEJtNRMhznNXaIGjxNy6DxZ6k74K1Z/fDgUKZXQRlqF
ARf3+2DBEKN0j0wEEk4JLcNCwokp1ZypvI8L2kjCITkxITacLTpBVmaBFhOhUSPiJSm5RJ2AdJYi
MU7yknJWzCNqKxizTtqQZgufqU9KmzjyppB2LFScLvSgqJsbqziEKRUixLRniG532W14SAYZzxZ/
49KySu9DlM0xz1mJmIdwpSARIJqKUxrakaMrUYCH9aq4L2saldsI/ob7KzIpH21pg+6XiO6eBSOz
UGTA0o8fpLd4XLdY+O5IIkDQf61zExEiaxtUnG4UGIU5xG21AMGwF1wMV+YcVyRUDlVIdIty6li2
WQKiKdFTuldHRQRuySExaXmGzwGEZqlFPq5qpZSo7ZCvkFFQisAAFY2fTLTp3MWdQVIgzFTG4VSH
Co+4B7Wo0Yet6Eg0iwiZuAp24FYCUVMERCqCub1569be+s7ap48On25tba7vH2xs7EjN39q71RWA
SyLMUB6k3D6J/TBwu4e7oipNl6m6FhqLRW9M5aKRB4iAUv9SSydqKaKF0wBbqIg4jsSEmBd3fbKP
OjRdCcxZJEaVgDnkT00a9iED3IHTzUsBMlSYAAcYLpZ9DoCcqLQTDNE9e0VFlBKEat6opBxMEqOn
MF4AXG60OQo+tIZZPJYEdTyrL4HTaNTNswZeECxjBGMjFEQuQvaAwzKlwPI5XfB95IcAFRYvlhh5
qjRTpRJEnE4RVSs0CX3GlS38bFFunlMU+3rZPESTeMIU2MMdnYPY8JslEjjQTL2hDkghHVOwS18y
VfIeFrR4VODUtSYBNguDl6juY7TIbMREaQnRLqQLvB6GwX01dJEV6TIXtFaE7k2n0sCoCYcU5eJZ
RIOVEpkhFFoDJnWq8GtQAIwHPBK4L4jOHPSXwpB3tANnyQOSbEFqNkKWlSqDCUDic9nloI8I8CmO
ytxFFL7I22chC6BOp7YFPnBFjQtREUiLtKKOMjRa7oh2IutSsAJqMChI8TFK0xsV21ljWnEwRIhY
O8oYURGF+kgynO/W3Qf5G94migaLlOrADSmoJQ9/SUbr5VWqLMzh98NMwG6RdrgzBxci65vyb7LL
FGxASzYE3LIK68jaO88hTdxSUpVbFuWNidWk/fbW/t7+9q21/OzwPw/W13d35Pl31ztV2YrgjNRY
dXfV7/n9JLHa7gd9MlxiyCpMEx1fqJfDiGkgUSvs6Mx9t+cDtbed+HhKSJKMmuTe07skvwLSS8xV
MGxZT2r4rXJBrowjIrkyVAaMNkEGe4liAkJoJShtdzKP4CFSlMtmOAyCSoWlUBtKnb1bBafsRi+i
MKsK0egjjsV9kRIsH3OnrTA2FUW0/DYp6+MpBK+yGZ8pt8HUJ0i5HFFJBi7wWy9AmTVZZrLgyxTC
PC5xWZxEklyg9L6zePbHGBB3PHTEE3MZTuTrcuqMXCd611R9P1zWFEVyViJvwTC6Rbi4QYczOoHM
LrC/Biy0NNJn+uJrsItBVDSL01orGXPWbFlYWFDPFCsm3Vch5nIRvjCVYnlvVceVUINyLF6P6hbc
ED5Q4Bw7zyZ/F9fusFLwkKfFhULxEjZjzTBzUGy6vD+d0bi5jNKWg8sWZA14IBdEICJ30Tir+ZLF
eTN/OCy5qAQxlDQKhoWJrnbcJREpi79Q/gWAsLBinH0FJpgDa6YJ+j4XKUfkcG8jSYgDqql3SQsQ
IihZjKR6WFR/J/mwxCceHrrcllBunE9CpchS9PAVPhEkNFhZRTUi6Ut4QbczmBwpNYnFnLmEok89
thNrGESv6rdM9aaJs2OVEf+LSsS9nD6l8HBymd3tvd2djb29tQ+Pzp69ONjY2lQ91NbWbocll/CT
WUx+EmVosoaovNaqAi93ADfVLBdGAq2NH6XAKgEV+yDg0AWXakFRKs0AGIBmPa2UX3UkJLzSYbUK
wOQHWmtEnpxGQDnCgeq7O2bP3SbDBRkekhvsN/kh6DolJU8u1oH1izCb8lKNRKYDSVMSWEH8jaub
WbKDyu0UtkWSpKbFhkyhdUnlK7hBkSAH21CkNKnXMAVZjBCiM8xHziIXlD76og8Kifvz4YQOB9gl
LydnOiPQChVkrhugIuaCC6PBV/Qm7SBF7+TJZWoQyJkZ1iU0hMmnSrUGkuN8ASmrpyXJdHlAyQmY
22ErRUm3hXBJZ4yghsRX5A7dWoqKbMg6s6YESOIulZxI4FjgQCwcjqRD3Q7OBS3aiI8myiJrwou4
d5a7nyzKbdCXqT5QrATUe/vJAuFFF7oiJFlIHzIKVY6+cVkkIE1UitoRFF38mGYma9HEmCYRNFm5
Uoy2n8xkSAiwCBYiGTUK8tLBGRnpjxccxKN9oX9n9AUnksFKKUB2EIY4r05tuuhcxM/yI9CDiVZy
A99IbuIWlhRUiatZiWOSxl/EU4UJWIm0vJgio8TjJprhApgHHEORzOl8j3h6jxImGpD5MskQgGNq
HUAMqi/EVXa3H2yufZY33l6frzja2dvdPNjY9pKjjd297fXB07vrVnFQoUQmB8cuU0JUkJlRDkAw
QcYQcyPbFXHhvTcxBEkEoYmM8CnvhpFNwetLsR3x+iZcIO45qQOdIi1rkfChscRLKfM2PVCURFXm
GV7foqsoIEVOyHRfJdxsUSKaoBcAy34JRUk3+wC1DF+RYlQCd9WTTc5Q0YlbDOJYxDU6uTI6waMi
CIAzEP0SiYXgdFL1SwlBFfzqRh5PEULzGrxtlqHBhViIwEqxg1nSQyHJJE4wW5F+loblCzbc0DrJ
GoXzyH2tVAQ4FnfICk1j1SAUyw4DKAr3jkRyqNEaLSui+amCMcbY75kms6i3SHMMfgBQwpnhWQqi
mWecNG8iBgo3LF2EcfKQnKQFQDtRENAF94XfzGIPWvAWmxPcEKABWO6OC6U+UoBG8wu/RSSDZZSF
4GHwt2LtVT8jGC7KABkewyArVglPNlP9YUpyR6yODEttijitlq+o9EQ0f5dz3FgJhSShtIYOjNot
YgNmUTR/FwSEIKOJCg0YT95N0gpaChrB8z5lqaGFp6aFUihhQDo8stwrhFByBTl+lvQKc9aIurZQ
BVJZNL7WulZ903gViRbgYeoZsuoslbEKTaMgG243M+gisBlsYV9VByVlTzVE1MfGt7b2tjd2lz82
vr25oxrTnc3Nrb3B4cs75tBQJIYL5UdU8hENXb5EgtiVcHYSqTXYSTOtlNpsXPUDKSaHnwIOp3D4
oG6bdvh4pcBuZuR8Uh1FYPJlkRcIh08an8UAMZYdDJUKJgQu8TDuh7C0omykqG2fTH8IJhy+lrFk
yKbU5QgTg8MH2bjkJCBkIVcM61IEgJ4zGhx+pC+QIeHwW2FUkFyzDl8u0peoPzUUDl8lHzn4ed1u
FkwnblKbAEMl0cECG4XDzyOHL6ktkQ7L4aexw0eilYsRHhV/znIBNqRw+FlZmDyQ8DTqjAkhKLZr
epqTYP79SSmxoJJ9Rg4/a/CDak46TKiGh/TJfyycOTBYrAjfBpYVN+Su2dE5J6LFTjn8bCHcJUFj
1niJuoYOumRw+DElCLlKALJ4l/Docnbh8LVQAiWVcPj0D3coZYKcxAszcI4pnLXAX7NZaA28GoqE
gziU1pHDT+HwQ/UYdN3m8EsO0xs7/ISAWMQ1QFbh8MvI4Zf/8Q5/f3P31vr63EMFeW9/d725+o3t
3c3NEbgn908qAQA5kTinuAvJliqmdhJJyjH9YQ3l0mI/KCRwsmjeVuEx9vWA+xLgPmslpLQc3Gd8
fWhFUOZkbYOvxwdlUXre4yS29VLgvkyDewsYUsa+/hLgvpQV4D7FLYb3WADuLUXiTumMkN0ScA/E
HoP70sB9kl8WuAfTY7gtlQl13iD2YWanwL2tAvfhxmbAfdAnUkxmwL0J3EesmgH3aQTuLcB9Wgnu
ZZLz4L7MgftIrq4M7gNTNEOIeODMVkcNFr8pyM2A+2DRMhzFFLgXgh35+mXgvozAfW6+vkyD+zQC
983XYwuzvn4e3OsKuGVE2UXgnqc3RKfPgPsS4L75ersQ3Ddf/98S3Kdf2Ndv76/f2t9f+/Bf072D
3otvhjvf6gZGXkLmwMgzfCNG3q9NfurMVzdQ8gnOk2riVkwIwZoWsPLJCzLGMy1Eqaqf0hZ/tihx
9sUMl9pB+kE5YvrKPkXMhopS4E69/yrOqkVuMjyUBNgxjFJMdGl6JvOvhWCUCBAYADbepsi3TK4u
SllcqkeUDIcJ+Uh6HHVFpQStrFFomuYsK68oJX7FQkYwG7HyjfvIoqPRo2q3tYjAXNiYBW0qcCaO
2MKFBCtvY1Ze0UQhOnKDFLFGbrBE9tWZqI0kQTcbVqqI5iuesCqmnRXjXepCRSTXUUYpkysRKKDd
W2aXZIVdEabObRoCiUutTSL/0VCKVOlkMVVFRTHWWPkkwbWJT0kVFxLGMeNOKk1bBinSYBLUYOWl
moTEWxorL8lY3RKW0jLMWThW1sVQIZholzmVxmFrUIGSBglPQ/aYghjyeyY3T9OsvBItmhD9x4I2
SQKNlW+ZTKjts6w8tm2yGwnvbmEKVhGftfhMZLfEAjdhPEhm/oh/0MojVt6C584jVl6+A588ZuUR
RxAVCCWhISRVxchQLKxI0cOhCuCB+JyC3yA3co8hKK8sUcS36IeoTrFWVCACHu/BCOo2SJW725t7
mzs7e2vp4Z2DjY3d7bErV1EdtIHYOybLl6gKRdUwK5EKBlUbiKlg9kKrKHJgYL4cHeZ3/6+TwDJ4
XlOCjxPC6EuADvqlGyudWNjAleNybkFKgJgicZA8DuO9Uppey9KyKnEoGhDACqFGNSHqiCuVtClD
pWATglSOWWKw+yAfL9XCdWI/okROREArUPa7h+hI0SPUXXcu0CMlYGEW+ZD4wcOYSi9DG3Y4W7vT
RfUYuInvB/+WR3doIV3Lu2RceQ7aRSIdU+kzpzrVHPcfyMtbEcRVRXkmJCM4R/Uk0VrhmBuW4XTM
fYQHjM+aeDDwXfQ8Boa1OSp6is5kap6DqmOggKujUUum520R/IpyuFIi3JSiiSjt1mTYYAJIwowr
l3rcIDSpJPCuDPA7y/v5PJNdaqxoYWRhngSryaL1I3fXkSMVGSlIJrNUMiUIktO5EYrBoDtSl0Ow
h2Kxhv1zGIfiWEyhRiHHeo6TY5JgNnUDSVNJUUuSWO54lIcXdH9iVCP9oYiYb1vCCSaiEm2UE9C8
Srojvgr7kP0kmW3WAHXhkAApyLDMp3ICrpkbDYzlFDnAtlSk/uao+UnSXFSdUZh9UETteIKSZfHJ
1MAU3K+cVQ52khQ+d7d3d7bWt3fWHp0//+64x+V7jWZZ7+RaKVvJRJJWUMgEtCKWVFTXI0vwh61N
mJMcNdiawtpThhu5DYkEfGJH97OQrc38zJph+TkIKYwzH3YqQHH62TQzUAfJosuZ5AwYIF/mmRUF
4b40KamB6af61ATkS5Tw0G3x10mPQzAmzKrGTIzPVHE70EMeNlZLAU4rSqq3OQh8bAbuCZgLg5G9
WiYXfR1eSP7QomwmgTrB6aZKSt9ANse0hY3L+JLgIo7H5xMBQBXCXoQqByQmPaKIFnOKioLmnhk5
h9sd41mHQN3zaJeAfUQfC8ohq1yQzML3VJc+IrpCg2RRnsISMNVsaPC05ro4ToOEA9EnwXHEcJRA
C0bUUrkcdxWMnvB0FsQQawGrpFjihtcNpTvChRZht10zonLK+sQ0qKqKSjmulpnbsBFE/RR+jVie
FMKipKoZECYlj+t/K4ZESLOcRygnKca5hcjCgmwZ6tspS8eEoBShELPuR6tUAhAAZYAtbbaKAqkL
sHFj5LAtF2Z4RkKJdJfU7ojiAhMFU0qEPi0PAfU03D+fsEhSW1VDGgtooAvyY4pbckIlpKiOKXUu
WTeeRY9xeApNK6vMUKuoRK1hFrBMcsLyBdiY+AYKHGpD3e3tjZ3tzd21j798cX6wvrMzcuRGmRoU
rkqzklicxN1J+5Fjh1Pyq3dJkp6GUcR6A7zhOZPKFwC9ymU6aYNBsHgGk0tEjIFdlbqVo17Drywn
BBiH3RHJ4yQXpbG5ESwl6Jh6SKdpypHqWkCj6bJH+SzxyL6ujK1VPNWjBswCBBYleRC9fIpfHkoS
ujiixNJIutWIhWPX6uPhc+M/+T2Hb4UFKWTt5KRRWJjkGgRPWZiRieCxPXs18EDkM6IUwZ7+K3ms
JxMhJRQFxVCX8c7qd1GmhCeUXphVMQXNkAOJpaIsLSphVCoB8FWWyM7BheSmOeLceI0cOX0IrYLp
OGAIlnCbpPUeYC3YNQjzYkL9SgscUco8yRPkM4QEobsaNZFUqKlV7zPWibYXKUJCKAlWZKz+g/j2
dRPOpYukXEYDyRWF9hhPEhNrFiwYdpi6hgQUgXSl0GBFzUIBgH8jq0y+Q7RDwoDD3IBskqFIsLgp
5JhGenRZJV7KL7AtE8+dw7DDNWuVFJMz8MGW3Yq6jABQQjYGHfg1U7ArnKyijgZPUuuhcgAFvJi1
1BgfHgiAmWaGIwHECxeNv+5L6wxPxIBp4YgYCAIIskR8DrpB5OkYo6foYDxPnRPhxBnCIlqq4PWz
kLZvO9177v21R98cPX16sLG7cWvw5FoNproPE0mvMgxchqo/UE5E49RZIgwAHjhGCpevl6jW5PwR
VQ6TGQXs8Eux0immxbEncaxmwvgeR9xDd4AgxjYFxQ85RFWLxfMEiViSNNo29uRRtSAB0MJ7kbbq
64AH2FCXRAMWqQQFYmvw5Lnx9ym8XyzKDhqixQutppE92BxVbjDYUOWEW9UJq8Py5CxNAwGWEIrJ
4BGpLYkMQ4YR6OGe8a1BoWqalMY4v5IszATZB4HAXfzU7NH9JBCS2LreSpCwHkKCkFf27nNk0X+k
Gu6QN2hwn0kBICuyWcAlOXlT7Y4y5YQnF0+JbKg1RnoutNLuQvypGxAlu1FGhQ9S4mBySirgYC3I
mZuMUKEvLEgxzBon5xFGTKUFBa7cy0/OUnBVEjhQ5dLKi+ALTiHAkvtTrc0E75t1/YzHTZIHRDUb
qCmsnMdCkqhyHH8qWhuFMzNkL0FQabZFjbMvHerawpmXcOYWngO5waS8BlUuwCajlTNX2ZvgTs6K
Gqk5c0auk4ShVsP4cGNZIVzXxshx5j71XUiyKVaF2JgkDSGppCOoclPakI2dhApRgaIRzbvuhBIR
wd1QScSEW4kN0Q3wjp9SBsPfYiRNgIFHTvd7IL726Nnh84P1WzsbI0fO0CISKH1RIBEWhowGmcTt
4SY7I0ALAhQR+srQUF+TBBlS1kAilRuTLIn6JZDgtItcVxtXAftscWOFQWDBU35TNL+cIKI8jYny
CCiW4oHAJN0qqkEKtQPkJ60kq6hNxbCOCgaMlxxGxIAbErq9bg0JvRmVPxwX04pc5EeIlxIJmXUp
biNu0osq3KaE2Ol3CVfOXJmU1lCHLOpM9RoRxJugnQWRdBAgwhAOVLKQZezyNBiGtYBaAq8pTZNP
yeLnEdJ9kzJ5O0wXjahofauKRNMRfVFo7MAGwXmEi8+ECTeVgDoK3jmWL+wrwUdMkEEEUKQSJZJY
M5JuNuVZ7AJFdqmuGP43hddNJifrSzFFfHfbU7dN88zSjLsQ+wYTSZwSW8DC7jC8iKoheZF4yJCL
WbDkBCt94jQYWLapTRHlEyRPGUbYsEUlT+0JCMPPJxMXbUFQZwupAp0lSWRiVelkSXojllwGhnnJ
54bMhO37poNRQaJoKkyi8o82qywBUSUODnl3W0jyWg4gjyIrVFGdbBYTIHSXTim3MIlyt+BMqLbI
GIlKUJIKY6SqF/EKJiTZwKggjXuiyJyjJoTIZ5HchnNKJmqdGdaVUmjxvSPf3d/e3txfe/jD869e
HGzcapLn5l6XYgCAsMLVRUsc7VpqbBKmY/0nf3UD7FeCQMgK+si7Sv0FVhWWgibxh8qLgHMOEI2I
SVAkJYbZyWLBQjXz9wWKjShSOQbVM087cw/cqi1yX90JC0E0CF8iRhTJajqFEQjV06npjlw5VE+m
qETJeZJmlpUQgPP9b/g8Mlsx7UMYaiJDVjELGkuWr3E2RK/cEQZihKKYTtQTnJZEGIsMNbm1Ryps
TfXUUgtsU5TLE76V17oDE6HU2EVYYekKOYWdQ8JZ6KME3nrlrHDtI50bRyl9LFQ3qZ6FDrtfZK9K
ePWs8DsQt3A/ZGsDn40z15VLUWw0tSNWNwZWcpWNWFGwk6shmC0/wIBIZUjhG0ml8IYFR+FQrUhz
4IZzEj4J+BbTnELsDP9oSJ4kRVmlDRiWhdhZxGOqO6ayBD8S7TA0ILGtwcxnWXKS2CniWkpbaBLR
bWu0conmPZkInYKG4MPrkqPbQmY6oomdwkqQPTbSYaEP6m4LTbbJTdRn2UjslAmQRTQFx6tXBrGz
hBuH7EKqEldDPiOVCPHOAeWM2Dl241mjgG9D2gk4TumLSAVRflKHYQEpApDMRWGfkB4URUdGgUIT
aD9FfLTIKwoDEW9Q2Zzfun5zY2/rYIsXqOzu7seGAZvbep2aBEF8WNGAoGHoZrPClRI8n8GuBHJG
B6A8Urfi7Wh2sPEijOA/sNO1hdgiRJ9bng7Hyao0OPIGAYwrE9QUlYQ+FNgyC7hNm9+a8hBHnUV3
XbQsCwxRDjlPJJ3cgOSrDJQhWaO5kCRLVkcQP3EgKRZwyAEqGKCuS+InxAZEV0xFLnHfhdQ+QdN5
DAw1imVEQC5UUElhQWOWxgTnGe6mDExSY+ghIJSIFchT3S1j10n8Ee9tQSpqcWUNt8dSvGuRMRU9
UhniZwoJTHpVgFkWS2TfOXi/0oU6JR2DTgOK5ZtJnpP8lhhtaCJxfAW3xBwltAEjucTcBB64Y+lI
HfK2qhuiU1ZCUfT1CFRhOggEGS8J8iebcLCDSJKDYBQNjHoTFhXslg+lThZeCbU4DbPCoEn8DF3B
dwyRnqEgmkRHZYyhFVJJ6dP689Y7hkRjkBVyGecCm6iMPOh1If1kPNkt8VMMYlhziJ8pVHlBMOZP
RtLQFbxGRgGmkLSpbSMJnqzTwW+XxSzlFgp9ZAsaCUtXGp5U7CJLqJVxSf6+BLk1Fj8tKBcVIMkM
JWd3RtWORSLgS8TEWggSxFxZUIsCZ+xKrrpJJYLK6k0BTe64muPt3d3dnY31tfQo5YON9VsjF07+
kZrEAKFSBJmCS5XaCicEIq/hibv16QlRIlmgtayALbYyQ+tneW1HyQLnaWBacsxLUTKLHWahMtN4
dA2cm8YumJYUWhbqqc2VJBrSlOTmYsqGDHgAQYpHos6rSLCXX+2sMS3ClClFKq9Jx0lb5Mo0JyoV
eSwHBaBMFzMhyIGzDA4CHtgvnCVQWIo6kawMygT+xHR7tQMgiRl01iFqeIqkZZF5khMlwBjZfuTv
4u556gYPbqWNk2TUxvtbpKbDHICelCkUlUf4GnEgZrIMRWHMTzEGqNeJYAHFRc0IQK1wzyrwYCSL
lD4PzV2WAGBiZ8QwhXxDti/iq5hMVp7YSTHT8oKBzEIzzndk9QtTy8puQBDFy7SttWiqGUBDyXKG
ZE6BxcmeEoWBQGTNRc7KA6KOjFlkSRDpBM7rVElDwVyHPAKwIyEEkELdjaoM3Hi60jgWYgO8LeA6
WYiigt4ynKQA1okVSimOwDka0mcTeUL/pmYFZyWIzRAGOBfOziKgBC7xJCibgtHdAM6TomEJHicD
ihlLTiOIF8klejFvCvyuVkvk+CXgtDibQh7p3oechIQbe1Wugl8w/Q/TgbmDpEF2LUBBH1txPCJC
4AgwIXxtd7s+/L+zu3b368OvDzY29zZHKDzoA1auSDQ4D4OkEOudm7dCj0X1JOIkWDH4dYoRRDKW
xpYriSQW+QNCOTXGUTyzBAkcxWwJi4V27EURKA0Z540mhiCTBCynSlhCIvSxE1/jzkdpjPIDCX5k
g/IqiCUQ7L7RrFjaIOGh6N0OYtVaXByOLGdpnd3g9onzwMqhHgJyK4RPhl+KDS+vLmI/4bWI6BbC
pzw5ztkzIHBx5Sog6UQqSOgGyUgHE5APsk1qhQ+p0tBWwmLMHDYHwlVGWJSJoeRn6XA5SFOkqpTE
vbFajYFwZG2D9ulIsJOIqc5I3Q6vi+NHUW3zL6hTU7UkUr6Eriv6mcTKUtCuRv4hIwXPdQMYKQqm
bp6sCHn/pFuF7QzBKfME4oBq4NgLLabQPslgwA7qGhSLZCuzKM9KRTJ56MWKPq10JUQj3qVUVFWl
FHVQb0oMX5EspXuUuFW73SBJUvlVCJ9J6mMjc4P2FWKWfiXO17SiBJ1C+CTDlttLqsyqrfJiNZyH
pH0FX8p+hBKytS6UJMsxsmqIzDSs3mYnYEEbem4ybr+mnuVKTaCOKhYtUhxNVLE04VNuSPS8SHlm
WgYrC4Ch1kzRBWX00MQks01CbxMe8Fm22qPwra2N7Vtrdvh178K3tmOL3c2NThwJtkuOKeSZIgSU
MpobLB/2vGNQSHUH35SSgra1QK3ClZKjUiAeooN4kdYcpI34MmE1ERmid3D8HUp2DhScGpmVFW+U
sgMWwDd0eihcSdQCZKmbcodCcW0goGFUeqadFxVNTHyYRc2I+2+Gi0WPz2uE5lCCKNo5kfDogzwm
OKJwJQvh8ewfEKoBOSiqEqwvTSVVY5EfkmBGXTX0KuqJNFYl1SVuvnFOciRltPOiUGNwkTh2AVUQ
VAov7WZRuWoAiqYhq/xPTCuVkmkoXCmDqkkRR6gbuEjVHmHn5HP0P2iVLDYluRoFqQkYLCr3gjmV
rcqjc/NDE6lg2yWKasWchYWVRi4HU5dkCbjbrkDssLblAlKSbajUUkxrLKu2ZIOfFwvtUw3n4bDQ
MwzyGWJoCYOrH8XOwMqC8P2pFcSXKJrSOUUoGu7EcaEhcgn7h/aQBGFEBJCHtjzVV1WRNgy9PFSt
INIWGbZ4JShBXV9bD0UJIgEKfSEIbSWeFtGM9YmX6azVnKMUFfGFblVD2asYckUN/tIeB4pu1BvA
XQm8lxgqvzW0BbBl9k2cMepMoBErb4JfkaplSd7K96HU9KBfRHUyIfxXlvgTBH53e/9WfX/42p/v
fWYHm+s74b83brEJSBSRiApQ3qVMWmAX9iv8pQHNVIeglJF4Hrq9smTJIzkIbQUKH/USx5NXqUKA
MBkQPKtQhbAOgeivLiiiEzS6ogVtqvawuXAbcL2eJlYsVpqn2o+CyZOQm2w8gXYxVh45N0RQE8JN
oeX6JVW5Y0oI4FR8KXZKTLFBlBZoDZUeRI2TkkaiqClyKIeBWVOMJHBiyJLiQuLBpYjkdzU7peBR
YG+VyHM0oLcEiQL2AEF0ulJqVA38b9HwuP82Ad7AuJzg+B1wzsS3DJrikJALIDJGGrOnMHo2Mw04
R0xjjqsFaxtEs7V+tRJyi9uTwzJ4cVCFySNA3BHrcAqgc2hdUy6TZND46hK2QaqMI/Vh7kKgVz4F
tuI38nhoWBPxGh7J2+zkfooEF2XT9N9arYduLQ/e3wxgU0RMNU4hi31uyVsjlvkQU6qg0EcZn0la
TORUR1jp4rsJhJRQVU+sEvIMA22p8RX0U25TgD0F8SRuHjI3wzkjeLF+c3Sb4gJZq/IvRkdcBoDM
5IChB90jaiVHMqv7UurkzoA7Myk1iPxMGIKxoJ4peyCJdJ0mhH1IMbEYynGZauK4rCpHVl8s6d0t
1r5QOq6xI32qyOv2zt7m9t6ttT8/OT472Nhc3xuct99LnnXhJVx4MKGSDCzct09it8x35+a8sziO
5rzT2HlrUQb+xrtYEgmVw3mnSzrvRClSbkMhz3B1590I9nDeuTlvSiCS9NpYkiPnneect8Sd5rxz
CTExizxpEUPQl3BVNJKpOe8857zTyHmXaeddFjpveQkttcyUybuVuPJq552WOW+pXnPOO6SUKAgU
DWPK6DR4JF0R/n0eO90RrQ/MkZx3Gpx3u0AeOe+y1HkTInHerEys3DMId95JjmhgTprzjtKNhDan
xHzGeTtn3WwjBSFnow5mcXSivEEDpUvhvLPYAiGahc7bRk8BZadAyDudRWCcFP9Sc95l2nnnkfNO
Y+edtQaQjQfnbU3EWei8U4DvaeeNPuX3G847KQFz5maZ806qRkH5Z1UMzju5yBIp37Tz9mkdnHdG
FlUZZBmcN5MmLXjaeefmvOGjZTlaH3qFdrhvW+G+8Roj9+0vQVEinVe47/rSiq2d/b3dW9vLtr7d
2dlY39Z+iHs7t/Z3w73vdCQQRTG6qJZBwyuVDBBljTJIypWpcCFpQ/NBZIImjEIKOf8U5iEjQUKC
7FOhDjcvptqYcJOII/k0iy3KCBAky1oXYvlDJpBgCtuTSRfcLXUmPhsoZ9ELoJxDHcui4smjxbPV
qewiU4BgjNSJQBauU/3i2nhf77zeTEPf8IZIX6bxtWC9SxJZCg+BLYkDi9hp4teUg1iwxOCMOB8m
OBByHq5GWDcVtZi4RikRkZEaIqWYJLOQEIJXtZghKFqlyrg7AkdnQuQ4NRHyWcGfWXHMZhKJMDZY
6o4wl6U4kaoEC5KsjXx8GavPG+1ShC3yAHT/BkSwcKY/p7jxyE87DAPqKGhd+pla7G7cqzoSPIc/
gmVRj9NiGQ4x9PCBEGbeuIrvXGkqBlL+YmJXG46y3BYDpSQmW8I8peqgXCcLkTORVhqoDXlHQwIQ
6KDBA0vBOJlihVEpGaA9UeHiRIjPfOiMjUonHBaT/BGrD1bSgi5jaXTgf5943R4HFdEYKRLWsEEN
u4qRtCY10SGcmFiErLuPZcZYYwWdBYAy4TZuLstNFalANMXF4Z2rG6rS497u9s6tW/s7+7PeN9/a
3N/c3499aPc2B11yY7MT5pAvSAycTzI0kySCXBrGLhRv1tUM4EmyRLIWr+STPED2AgoVdA/63eUX
ESyD7uETgmSlViGfYFC1XI1KuSLxFZmcEJ3CQqWbQE+FlE6YS10olCzfGEUMM5vKzUqORxnl4elP
p9XvaRrQnvqBrFiN5alKJPiGgpvSE2tZN8RtJx7mSJI1kFAZRLM2Ji7wmWqW3FUKKpVw5kJLAui+
8ILj0OvrU1SBlEHwzGVgV0IVF8NqoRF3IVMIW6stT0BKSHPgcnIZ5eS+OjuydsMmSClxOCUFbddK
4Yxxhk9NjjFNmD3FTKQkIIxoBjIuuqwmFEVU5oo9UAajKiJgbUyPVDPTZRUvCoEqWbiwVqnpTYYM
LY0qhV5Vv+6yYIwT4YyZ8jC/mulo6HSBTtaKb40Z/HarmMBRFWGHtjSKWE/3VW4qXYpbagbOxdBo
FCq1hgOcwFDxKGciQgrXNytnKXtNt8KM8uaERTLaYBfVlDrFpVnBPAeOJhhzTvanbiNjaCRPYH/3
2SXG1qSnyphrZ3z/VO9aOGqsTbMVkL7IuSl5YyXQbc2HrirBFwoLuAVVrdIQVk3inkUKqAaT9RKw
LIWLwTSYec4kESmAFhM+0/FJ5RL1kZ2t/b2ttUdfHj7d3bt1sLXddkRZ3++KwktWLGB6RYeIn8dB
56RHUVw4M5XyDJiYZyusKSdBAzfMDWRGLarEm1uNpCqMQxwGwTogs1/XB8eUpQkyZxVIqevycCUr
nAiEkbaa3F1dvx0EJWJTEXwxlasV3arkUwvILLfhrwdSOVyKat5R6ZJR3u/ZLt45BWoo5G6lEUfU
+CfxpeRrAN8GxKOoyk9mMSlRIl8X0k9CjDn0N/lJVK167Y7qlUGrKCSRwAwe8oBa5LkGP7IINHZw
SUG/JUm9WDRVUIlhD4UoawKT70uctOZSsF4JbcnUVeUgAfYjA3IJVjhdREaOIKinrQWTyNMlfTha
rM2x5TqmmXOrSo3KLChkz5WTHlESSvbKIeHlKMLIiuyNGk0xPRZWWBSzPM3WkKi0LYVynCKhi7oC
d40IF8qYSbNNsh0NWXjowqpIShDjpuD5M1I7C9FtEFyJiZmI+kbVaGEFd5WRp3KEIuqXiT3cmx7R
Kg1wUFYlzk1Ukko1hBBC0rV4EKBILYkL4N46TYQWvUzKnSHFBImCBiXcsJ5itCKLhRoD8wp/FYRK
xadcRDeoCNLpIcpRcoNwuFILLGro7VEVLlwBr+eFzRagCU+P3J2QSy20dOXGil2+Totz4bKBKA3M
zVZ7t72xv76xv7l271/Xt/a2Rg4bzGQSiBFcICpU/yEKR9op+Sn8U4evZ95LPHDgJB80LTQUVRu0
l4qKRaBzAVTx/EtUAWpRYZHB5mUYHzjiDiwJCwjXxfgYfhTzgssrRP94fiV7ORwKcQ4JF99F8a96
Kja/SefCD52UbOCgQIA0R+i6HLS1BgC4646QLTBTKgGLs0VbSpiBJ6E7K5OnmD/gCFCGtWxRU6yx
x8hlkIUl61PRibMqwo3iX4to6tBEi2os9EgeWSoPDYq7DJLfipB2CrI1x2M/4kiVG1H2krV+c6B9
5hXa3U1Q06TfGMLUxVwA8EuUbYDJ4nkTUhV5ohTUMvWDeh4Dgwzi2AlYDXvOouabX/f42Y1mkrNG
1YhZ6q4SYg07R0a1DsPK0ogBJ3kzxQZo4cJDDURp93zMWZbragRkVmAuId/lgNrEGB8OL/XJKNWJ
ghF568izilhS5UVUXBDsujysUifp8V5IIzCPYDhFEKydugrfTLeUsDH+DmjsVsL4FUnuuaU6XsVd
VB6cZcApkowsAbZQvA0xFPUidK4roQ/CPKUAgCnKdrAD/0IWL61cGQKPB4jdGKA2RhJjLyefeHbI
pZTOp0wxGb6yiMbTJKckiFvkjHGO5qOdebbM9GhQ1n0RC25vba9vrG+uPTx8fnawsbsb9Mj6bgfT
3cSVEAuG6qusaSbMUPWU6I/zlQICETEBtS5MB0nUCmdytOkkmLTjopoRJe6SHXH5kh1z+O9YVFk1
oqxtkSKMhpQVpgHgklTknCJaICY1OSGp7l4uHLPxu3GOY9Ac/Z59RE0rFCcvroRQFTT2tOYItcpW
7u7dRpqjVCacd1AaUcJoWkCxx6bkSCUehYvA/5TgVOW2VP3k4SZ2NHWYROTIreY9YBRXdqPnVHfe
WRtnmGoyINKiEiQ3zVEJpjsE1q7behfpcYKfVsNZvFZRSA/3j38q8qsdsk3YIAAt4fiy8qbIK8HR
ootqO10MAwlQ45WV2xBrtZCSwk1WXmtdHhf8yYwGp5LaCBvkTUlhkcWfKoZrIyoE+Zg1ozHrOYUL
B1hiEGwBm7UuJTCKF08yxBTuJQoZC6iIty3LFhufKcEut0hahsFLyoAKgV0oGqJO4C0LEYdUQSzL
igsmI2dX5pB9oTaIsxlqoLAgk1R2cB5BKWnXBxwprDt6DAi0MMBFeiDpYGSNzkunOAoYlXSalgfj
ELJ80Zi7J2aje1iqALwh1OBYuD1ywBLqG/7D6RHsX/URUoTIlQvuF2a1KJKEkukADvtlUSJx+K10
t7c3dzc2bq3d29jdPtjY2h85bHHnQkbqo+GplHcUMXcUiClaOGGo3hbFlhi1BAMFMQ9HIoYOJaxQ
4yLyNoUulyRWFRL1VEaqGMqSZEJXp2ANOEcZflMNs7JecsEicRIkTW6gqw7MXhMjRNgbbjCLJvdM
OiOqAfbLUM9B2oVTUDUB8SnBApsyYOtEEHvUNfGhhXhUJJ+VSPBNDK5WkXY6itBjFq41605IDUW4
kf4hUReLR+VNsi2FLDAUkuUjYRUgEF0JYEFUc0AEJgpAYkIGprtRV1Ioa9QDhVgoisnkM0LAFO9o
igT0PEupgYBLglAijnPIPLkZC8afQyP0ZvQaTJGSEoMFxySFOmRM0luFIBiLLjW1URxkIW1l6Xpb
UAZ0JOzWnQkPrAj9hHUn5kWcnYqZ4H9Co5VWwxkyKaPmgOEiysLeKM+RHm6S8yySLgdTTTmUllmk
CJOj5yZqo/zEY2CDTonzIc0vQLWimwm2WzcfrwlDapReGjpZlpPV/QO/inrmTLr4EXgnZP4kmi9E
WvdA7tPa3eKYeDxU4kuBmJBhy41q5i2IemkCiTq9HAwskKI0LTNIOq0hyym0AS7MaEOnhDwNQxU6
HstiAIZRVOBIRrI3Wq+SEzB+d3tjY317a3N77eHGrd3eY+/ujjw2hAHWOiRvpKtKyVL8CqgXqDN/
ngkXWTR+gfw9rxM6EaD1NUv/3CI7FWQJlnqIo3RC5CzMDKm4GEqIR0+XCfmFohCVcISUL3IDp+sQ
g2RKeXgnphwcjx2a4mxWpR0EHzx6JDTeWbZDEhEWaavpNqPQSTlJkZSXxQzwAIB47sQ4SqpE7ZXv
VlqfxVQrXnUpGC/5VaXV0IkMIxRZUr4wwveo+o4YLaRc7zfzH/MlRjYHOmCyu0C+SQsIXFCk/Wfx
5xF4jAaVAHZRSMVomQ3kJmww3kC3JC2wqPIAQkYqbIo8OCc1H9KgUkJSYWJP5jlUB1Ie15i26HI4
Pvw/HEATRXy0ZQY5piGFEEU/kdFK5IRZt4aTjvftyNwF8xoKy7GUFFJSa65eu/P2ShESk1Qaao2E
uxxUXS5BMvva1etSg45hhRXVkuZIDyGzGVDSJ59pXleA7ZCAtHN1q83Is5oexETfSwls0SzKrV1+
Us4ENYgU2JT5FO1aL6VCIyfrTI3RxObEYpF4FDID3FSRdJVDgyQRzjmq68H1yi5yLEmNrRJPVgn2
JlJYDmU0qKgAfuWcWzyG+G/X4z+tBGhzXBgCCXwhCZfhmCJGOsTe39nbXPvk0ebGzsHGdttRqmHs
9DvG/h1j/1IY2y6LsfPvGPu3i7HT7xj7V8TYW7e2dzd21x4enz09WN/ZGjlsmEeJouI9C/Qmao/o
lwRx4kMj5Ji0ZUAaoWvQsSpCYCjbpAk0+y0Ym8ObqMkinFHC+YEjVL8U4MbCqcOfSwuX7xXxN+Oh
hTE0+W7EGin/P6yupagvmfbQPg5uLj6u7p6o9w0PLU0cWyvEYBwG7oMyEFQi4nmS5yilVZUAigeI
rbUt+Ql34zIHkNfkngsJBzaoLET1C4S2AsCwBrGjaTwe9wfERnEgG2HViG31eNwpJU6uNWmBwvXj
SXBqTJAjE0GFHBS45NAYGguIbRIaArqaojOZgENsoqsRAFlN7ogyTjoghcZFhKMvW8qUjX5mEwcf
SqCY3kg85H8U9oytvhJg0G2+RDiiR+gMHrLQ2HJkOgGxU1boZw2nEcQmbio3SUUlRWDD+pAUTG9o
/jhpnzwpu4OwI2jD1UcQu4wgdhG6VRpqkNbTENs73iXJYjmHNZN3DBB72kkrMzAgNgUwWDWaHd5I
wTI0Tp2Ocl5GEFtDjENvyQiORzUzEVxb4K//dGLwFcFHHjqVACDhoVNUpAwQ22R96By6WhaA0pcU
nWh5J0JBBfeyDDd7FwiEq8mo8JVadNkaC44n6hQuc6CR5r/r28629rdu7a79e9nZ23Ae+3eM/ZvH
2PY7xv4dY/+Osf+nYuztna2d/Y21u1+f1DcN7+6MMLYCjwW0StQ/SLGFiKI+x6dRKMd7z2Z3RZwi
JUgWiq3wbYt8wlU+zgX4RFzBepOqHVQYECDCPVJRvAc1+lpEOpbyLTHcoJd8oCjpUclgHmRgH2YV
ZDnWCsIM1NQYPQq6EuhGBX0GcLPgz6xImU6w2aA5ixOM8oeMjk+NXQFstyoOUViFqtNWqSE0CY4N
Lq/LSm+o7nNbaTUiqRVXRnVAUikIBSwdrLPbb1LFSBFdL6JYPlaFTm00nEy2wINJ6ZXKOqjFNNMI
pwAKjZ32IEU9iOR6NYGRqCgNJ0081VdYjhZyiRyIugj9z+0fHBxQnNIT1dt13I+7COFVekgFgc99
0Ww7HEqqCahjxNZSKkojYYiOFfxsCqRSJCikHEWTnQquYPSjOiTozUAifv2kGWbhOfueo8bBdQKH
HLDHKgsxJqK5AO7RxzLHq5kZ7IQdC4bRaArsJPAIqe0jzE7ZQBgAZ0YGyCq0CE+L1qGKFZkjbxrG
W0iisSjw83b9LHkGL0UxazUcHcbg9C5epjR7krX4FKiSBTBD8szb2KLeLcSlUZmIIUcwTdgcpD6r
ysLlFNXfyDsnjLkAp0n5FMIz7fEkITMo7F0Eq6mmYy4sxhfNw6U9r4IqpFCtjgYY74n87Z2dna2N
jbVHZ6dPDzZura+Pdcc0ykWxWtX9Ct+qJIVDcJBETH94MTIc3aNT9OTnHpJK0CaBoBRjeZVnooIn
bkrBAzjbsJcEMb7kSAeMwlu0EN7agtsI30KgE+ql/S7Sl4yvpWJevUgE4hSCkEgADVKJJ6XcauSq
k6KOAdYjjhcNiuTOerFO1odDEbMgSkWpDLE7kkqLetWkF6fSNRwzBVmAiPCkhKehvkqBNHeEcHJn
XR6oB5rT2oYYACCVIgrDd45pIDFF6Z3ICi06gJtSMeg0kGg3qj+kEQ8Q6EuRzTA3YFdBFka7ADcB
FQnzU+s5q/8+mUkRrEQtHuYpjCPzydYSn6IViMykTE5EGw47ofr6AlWFnmYra+qY4RTRJcWU+EMZ
SgkBtcK14hhNjElCf/UmcoBv3xRE2JPZK7EoLbIQRVbB0txoyOh2FuNBKkZQw66zeikqE95C/jBp
z19ZZFLm42mXrxeuho0ivVkk9cOTj8QEeC3yj4FetWYvWdV5AglduL2sPFN2Nkqq8EmaYX1OOgaS
SaLZipYxrjSpT+JN5HuYGZ8Q7T9hZWzdxMTC1ONhYGrkyw0zYwvb4Jpwj2JfuUZSaFPeLa8i7qLD
P5OBqXrYB6Xuv7q9sd5D7LWHjw/dYY8KReT2te4dzwp/qlpAoCDlKLgoWBEEpUGDJJGDRWmIhaGI
zVZmmbPmQdYRqAIrTsr2uDQRNyv39NzSrYBCEZRV+WeKFYJky0NX5K2bBu/X0R6VrZ4kEYNL3HQh
+wMAUGFQVLBb5V3T1UgXIQhtGB85CEZS9LyysMKbli0KRdKoUKSMWOxMAYWRuVmYAu+igQcbFYoo
OQEIAubkN8RiZ1hscW9tzNGwi2pPLSqBIZ7GLHbx18kUajoMoi5Yt6FIRz9JX8/jQhEWSyurV6FI
UrYuKG4mxCvITclSBwh3l+D2VEJ3D97egFBDBU0RSvFu05NRoUhOMW2YfBoKRQoFU6y/0qnKKYeY
UMRPqFBE2XZCdxH1QahMKhSBdQArq+QHciJHXbeFZ50uFFFRusXCoIOwxlCrrRAlCkXKDIsN2WNk
QMFix0mNxZZxiu7z3RKnCkUYOHwkg2SjQhFlO1DuM4UiKuhIwoMlJq0VioQak8VipygUUZ13aYUi
Nl8okoJac+Ta5eBioR9aoUhZXSiSxiy2RaGIhV0ljES1aq1QRNyN5zuNQzKl9sw0DpXRKShOELI4
LNI9vVtEQpKIAPW2PqR+a29nb+3h6dcvDta396Yq+1iGQRJkPXOC0ARmKaE9Wnhrd0YCySW4QLTH
jL1o+gT60FpUiUdM6Zo3VwFT1GkhO2oOclT2hXCBsXShCY4q+7TyTJ7YwmHnkcP2AdHeIgqZImSE
mEh15LCxXhs57OpC8G8pki9SAwcqQ2VfOOyMwwZJFXakV+FTSGcWzxCMZUc8Actc3LRvBCf32urZ
ZNXGqs6RquJeHNXByDTZMYrDYKhzOCqfgNAg5BL4270A90iaJCraEFSLDVIWcGnksP2bLrJRcEc4
bN1bszU5QqHaAmTO2s5YCdSUw44Baw6brCwPg9eJPpGVaRjJZFR/RtavhUO0UiDzbFcVeHQSJjQc
tlwVxgfPhM6D7GhQyEKjSUsgp1HkwqgzvOHgsPXasuawwXklWGeme1TZJ14mqwizOWzJjoSdscMW
CSIWYKjsY3vbFg8G2dFE2ZVllX1l3mFrqEN6DZgZDluwWaBktrLPGlEZDjtSBVmRSZIloYzKPgu0
WsIkUxGwaw5biRhrW0tSVMoC2THnZpV45OHOHeZ1VONFlqesQ9lDEYmLZWrNZco1a5vgCTKckgVx
UXy725vbuzsbO+ubW+wUNf8uyd31je2DTd8uand9fSRKIj5ACWVhQfdPWJzxrTIP+QQpMF4wDRE/
kCZlqB/1MSgSJhRDHatDX8Lt4CyzepFCajchb5Q+ZH4+8ZtnlwESEqbDbTrg76CGDJ9AVNZLdMpf
BO/FA+XSKrxhNksWb1SCWawNdjoHpVN5Tg4lEV1cFDI8aGoZmLYcQRVDEUyq4MFE5NpZnCaIr8IB
4yFOGBdJOY11tCaaokdkxjWioBkbBRNyoLJKJArJWkKBJxGhq+zcw2cOBdVHWtpT1C0jeCSpghRX
ZBgysL94/yTNQ53k4qENW+tCiRmpp3QsLGiHFCZRdKslh+SXhJjgZhEyEZ9zUG/EZRvqZsTymmn1
8kMhLc+uaJbU4GloVYlL0eHEbZF3RuB3dr2U4L6wTgsXW0xmRUgnX8opPF9mdwYROPLXhn7hd9lM
IORD5abEbkbblHEofJv8MJOZ1VvmoQjTeErRCeUWrcQikUlDDemvqdJpRWSob/whUa/EMImV8LQq
6LUSEIakzyjH6WQWUB0WMDkF9gZmWThduQk5Jr0NzEIe1WRYsCOiIUO2hMrA/n0XNYXrHLRGGlFt
grAloKaKJdCFpFzTU3GZpJVZYnjIZzDNWviiujpI5SQAjyRA2Olub29vb/QI3E6OXhxs7G7sjgv/
8qCrOPrBL0iIS0EFKHNLAvg+X+1dQTpxiPDwmHKhqM2OmoLSNJ4bpaKGbEhAIfJAEYvAuBK+JwUC
j6I0v0VT2i4EalE2orxLVE8W2QwCx+NGQVuBE0bhM4l1IxlVCLzk4T2Ilgf4YcGLBGWSgmz2ESIt
MlEmoEGVEi1A4FncYpEUTiqufWiV3k1RJu4QQY+NMpGnHlEmxQL4JGXGWQqpcFYSZVJESIwoE0W6
yCkt9ONEepgtSBaqf2Yok6xSI6UTyEvIRrTYELjypRzsGjtylai/S/ScIMWAEcuCMsniv5nIbuiR
+qQkIJKgDI8hcJl185Hv5CTcAkkkhW6EwHPIF95HpjIQOF5T637u2ZoxZSLpGdZNicOIMslD5QhZ
exK/plAN4cgxaUHhn4x5OQK3QOADZRK8yRQCT8MgiRDJgcDdUDqY04Eywa0qSwwaNOWBMsnCWYHA
FUjJbX2khIazKKsSCDzPPFuTF1EmnuxOIXAsQ1ilIfCwvjywmBhJAg57ZwcELgfrY9Ehd5pK/1RC
EisKzlTJYxYDLKrSyVecOHasRUxQ6m7vb23dWt9YS9+cnh2s39ocI2xcm4XkIYWCkIDntRRlMqba
GzRovVY0qWZF/ilJA6ZBiCj11USi+hhgl1FiEtl3OCgZUqN4EdHb77E/cmRQOSQ27560jFJUVIPh
imn2N0OTUxUKgcR+5KjbgZAiEBFBhRfrpTvlU0PXouYvqRALOUcQn+o7WTY+V1oyoyOQFx5AvEZW
yYQhdGaKX8zCLhA3iCYpzMJMsobyH3EZxO8ukpbwES0NYVET9i3YDA40qcGcHCxKqK9FiUEayv4C
2iRNZL3GuOwPnd2NSOkt1AJyEuAJSCtI1kFsScF3pxTVcCIxwJaQ1xQuyZK5cg7lXwSOSZqUOkX4
kPmz4JRyQpkIOpVIibG2UlrnMyJuVvmAeMFOg5Gj7M+s6Z7SfHOU/VlqGQ74sYuVY/LKgYeUAjdZ
UwKNx0ssijf9JuVWpZXkwS9II4e0VbJCbmnw3V2KDEdkEkkeqosVWU4zGICGdNvWbStNgAsLg/LI
kQhZRFa66ObZ5gR9rMj1h5KZBLAYNWF5NxXxpyYYW/DxFvJQoznAKoTbPIiWlP0B5iNzNaXTyaK9
bG09Z826IxynQJHRUhQoWriRMAgmnyVDbY93pXQqAFRViIRBdz1VlNzf395c+/f7Gzs7Gwcbo31Y
o46kGW5I1tC2NB61lu6dRMn//HUkaXkdCb2RtjxTR5IvqiPBT0BOkO10Re3YuI4kbjVJU85KUMyC
3/cl3XErtrCOJC2tIyGnWlxHktNQ8ndBHUkEs6GORMgtDXUkDekpI3IPQ00F5FIKw8EX6m4kuWel
kVN1JKEzK2HM03Ukojo0Jug4c3UkVhRAVL6BqyLCkRapgWDF/DY70cn4z5V1JAnyutWRmOkVL8mi
6qDVkZTIUVPUkWhcIQDHdSSmDLzVkSTVkaTldSS2tI6ktDoSWX9WNlpiLHxhZIs6Eokvkdm/ch1J
oY5Egwj8FrYXSnNqgjRG3IBKuBLxJlR2qMxgdFVHkkZ1JIriRXUkKp4Bgpi4spSDhmUWVUdC1GL1
jutIcmlOP8t/qPxmeR0JGUSKOhJbVkci1mS6jqQsqSNRTSTBUqEZycXkHnPkBOKfJJVcWEcS4ZCv
bm9tre/t7a89Onx8fLC+u3ErHPZWlwTbw12nKPsLOwEIEIXk3fDQWVukxwILh43l54bKww8lQpWP
ccFho5mkcNhZ9FEOQEiuIIcdScfIYZeRw86BqJvDVsmfEDPzWJJeDo5yGQ5bvAbadDhsGZvnnTg4
JyjhqeTxLIh6Lc48dtgap5HDRiIjWbNsMUU6TAPFVKfmsC0GjNgqo8QXTDnsAr2lqZLSwrocCv9S
GHNuDjsq8qDd1VA47ExekQYYkyBTGzNJa5Gc84kfWii8KVGmMxCLIRhSkQKXxzkmmSQpryjySix5
C4eNGA2djMNWOQWxSA4btwLklkLcatNMZHmWzkXGkwaHXWKOcckmv6g4je9WHBDQMPXCg1RB4Xb/
6Y4iS34kcmWV7gZlL548xfNXKVi+cNh8mpLSyCQ106Jl2pPPZdVpSedhkeXm0bDmIodt48I/7jlU
ibaCs9IgFrM8lzMZWpJFYxxoIJZmE65EQ+ck+xXyTU6vaq2l0IsJ9xqhYs1hm4icLKsweEqLJktj
BDRdCTWzAVwjVvu9qAilBEQQIDCJKmI60bsUrMX+Ol7rBmYAJJCx+RT+2QLphRgbmXxlC6Omw/kd
zLx+4Y9Cbt/aWEtfnp7WZ9envDWJrjiyAkQL9rsFnLhT6uTUeJcEpgNgIXHzK8GySGMpQnvwHLWd
8NZZ8BomuIjnE03Z4LVwnDt8eWuLBUyqa0TMUHGUXRdBeogRXznSykukkcFp+W0HeYy31uIoWfad
CrV3Cc+K+TSSVKvYqNUi4/NVm1QMlJr+xiGKUBiPMl/xa95vWHrkMpVpB+5EPwOrEbKU6NEimCcr
OCRXrFUkA4IsEfdMPEgKxGdWmpdnvfqrcVL46FIiFMHnlTwYhw+glkWR7ejBhZgJslEIBClzwoFy
A3JcUKE+YBi5YJ4mLB4fyXG8cjXiAmloF0gzKylVilvoSMi5rBYLvgybd2EaT0P5WtDRltLUQoQV
BTrDqXqOXuQbfQ64VVRT6kCMpKlErWW05HQuOYrCLjcphUdSsZwK98SHRWEWC8tSWQmLEd5I6ACb
8h8pxtNH3F8TGHdTiuyZSqwkolZ5aQm3rUDJm5tA4PQnmoe2lnlZuG24lhyRo5PHg5ospVFcrAr0
Fg1Hwf/lkA1LZ1FolOWb5SUwFCJwBFo9FxQlLqljBXG0+hNuj8WQFDGB+6LRiF2dkEIOT24CqDoi
K4hnKDUWgvoK52aMSyYv9pnpsfX+/s72/tq9WjSyvr67P3hris2TsHUqIbiF8yOA8y0BSXJd8Rc6
Z6FX9QZRLyhwxGL8rZLh0nxzB8OhuIgk6LyDz9BIghf16G0nRADCKUSwMcziPNDwxMAxTnB3AkbO
X8trIYZAJpqIXeYmifUqIfBG19r7RjRt0uHpZGoxNUepCRbuLhv+ugR/rRWDYM6QquoBq8c1awI8
l+tE9UujTlG3EkSjuFySSBCj1IEK2GCoUYdamoxvijkyVZE0XSoLa3UoRVmySWNOuU8ZiwhiofWo
J/LXZ2dwKWBTrkuQXsKySacvEZowyAzaU1EbaV2hb8GrIo8RKkhyLdJutidg2RepoqZcP6tcmtRn
qMFQ9mG8HFBQPFGhDscXBWVN9DXNYBZ5nXg5YOOACYfy+8CpIK/FsYpNTnFyEQyly7LWhtpKDl4x
iusMc2JVpRF5jYSWRsRckTRsLDLI6wC24nByHpYdiAYskoPIpCMZhQXo4gsjhbacxG+VEXktV03W
zJRKecxed0XeoSqtXCKj99mm4occW+R1kErV7REDVTBVrAn1kvDhoogx6MMpllVh5/0S5DVlaHJS
kCCwglrMMkw8ei2RKxK2VMDW+BLyjwwXIt8UPl/YtdNkFNOg+ERV8npna2tze28tf/rhwfrm9tb2
yFvjZ7OyLWsXSy0xUKYWYUOB19q7MqWNNeIx+MjU6MWkKqNY2/USnYU+4PNEU8UaAxECZ87hTLIq
kIobpYkHE7kvqSb+tdIIT5EhQDJPNXI4QC4getRMJBlmNnBmAozeFqKuBZ9NgkNinFUZjfYgAFJa
5lcXRlfENWL6WdlCFGVYCTZAmN0V9SJI2VkUwQdCN8E6iJZElVXI5sg+SS6hE6esxKUozreUFCwA
QdPky6D+OmupYpb06SasHDiQfGSUJVIAH6kuuEaWEulxJAUl8GEuQS9gRBEn3NoJDUXpMUHELHgQ
LZQUyadbcB3QTmRiIH3dmvyTFljIvOSDUQqSjM0ckhIgueZsQXc3m5N8lRoBVg0az2VyMDTIw5wl
cG4RysbUhER9wMknmB7DdJF5rSURRXgG0jixDnxousiSTBeE5QUT4QvJgYtqLbRmyN4sSkp1Y+Ss
0FVZyYa7ldxAbw6SqlPpgyhAkzZXiigQqADlQtTsG0Er8xQRFlcoZTGsLWgkXGJBZCaA4US9rtzE
mpg8TWnSQEqxTphiw3cl+SiXpItQaSkhaSVxg4LvKUCXAGFkcnVrp0LhDJq0ww5iXAi1UbUUYFLs
T9aSxPpUNiANXt56vffW+WB9a2NvZ+Stk4U/KuGaU8OPAmniNqLwiLl2MlRDCMsVMMeDGkUPjXQP
ecaRVe10B3YMeGTKDjF7CdeEdvRxvB7ImWop0YwWjLaTByr6ZtErNcrEXQocSnvho1GcJPXU0z/B
cAqLkNCpfZKtOgpwv6ZllCEmZfluq8XEfSj5xladR+xCLsU34VaodxDlpIiexiMicqlrISipShCH
IRGrBEoWTIkABxTl1VAsI+AmyaoFulQ2jlcoJmkFdNZpYQXrk3F9RFK5pGE6gmfVIu8URgTYWPuy
TRXSAaIsUjKpVPWvroRrVI/EyoWcQUCI2gZlnFDBXlgSjopghjPKYiIKGoHJ56okIGlouhwBJklp
0vAKXeKDgCXkQk1E4A11mghMSvDFNODeWehcgfxAYNkfnkziCUqILUIuprwpK8fETopimZdOWo6F
Iu/RdA3jLAV4kVvkFQyonsgzwQcKWoML8tkWzZ3V+8C5PtVdkEosJ5ipsApJUP4rHxZJfbnhviJI
arjJHHijhGQ75IAWaDkFjAniQjwm96vcogW4xiRZCV/qU0XyUyT3mw01C8o0coq1Muhu4QykNAW7
RJgRCZAbM2Lqd2p/e923ajYS5+IQu9ube1vru1trj2x9e2djcNU5z4uMQeeX1Mg3oHoWggBx+6Nl
up7J+0Ft5aCeyR6CjiRrEffVCdjO1YRwJq1YGbkC1oy7kXjDuD4jL4rxLZAtajSlQWIk5HSl9ZTc
poSBFt0YmV6bctPZlqiKVZhU8B9qQoAzTWK0piZgwTM1IQIhKTjQJjFawIIZiVG9DRmlwVu5LGFD
EIN8asOwS2pCBJ0sGFf5gyTlg3wj9hZJKYjokBiLNVaPRxlS+Lvg/XLuYt6S0LY8lUj3BlDDgpIk
Rj+rawUYoVAog2PIIvFrNSHF2t1RE0JEC2Zd7FwqcaGEaGIpJEakK0NizFqRIREO15XeWALmKOvk
FltNiEVNiMd7SYxSV4OVi0EVBDXXJ1PQVDnSZ92SaeLIQnV9mSolDmCqZDG0gbNk1xq8QMpAhSYx
6mR5mAz7C+cQXkiuxfMB8j+Jy6oJEaobu3hsKe4piavOKfS21GpCkrQqwfWMocGXA9VtkBglGMCJ
xiQrbxWiSYyu6w2Qwi3nN8HpeF9HaP1AsCKSCkQg5VyZVqTSmf0MsGOjPeUPEk7wHSXwojCH6CtE
7CBFFBaRGPfWd29t725t7fPY48bcY4+b9WHHg21/7nFjd299c8SUyDJUecP8DmRpq/0oFk/jOS+G
CXdFqr+0JVUKlyjEzxFtU4j6xeIavP47jBq2N6N3JH2jJFio3pOmnKg271QXSltJmX8RpjHWMQ8y
xLM01ry89E9ZcxYayo3xdCyUA3H7Sozon/3JxxwEiePXISDhGeihCTqlNDwEIWq6KegJrlxKO6wd
HSWkKNAALPzNhNyJQhOPAeEUZK4ohKMooYwowbRS34F/lR+DIBGxHopnAqVJh7OoBMVCIMiLqgUo
18ryagyh6HTiUUfGm0XoCSKLXJSQmZRlh+vJesyhdPHMHuGN7IZHRuQflR5laVMWKAkX41lRDq+d
1fnQ6KSfRzGIXA4EeQcdoxRESX6SwetAHiuJq4b8Zg5WME/GJofmFsg3l+DB4DS4fT3S1oXdcW8W
zXsDLd6GqhGphBx3Jy9deBYqUhqBazhx3Joe6GPwfcIMZxBClCnPjdOsTQEPfvgVsMna2Y5bQN4K
h5mKnsDSHWDBDGNb4w5iRST53BcFDz3IQ4CAUaIyAxInvAqPBMP46fSIK+0mTQ+wgKQULiL2loAk
OcSYgtjh3oHlnZUWyBMxuh0f8OiM2N5gWGS1kPcMRg4Y5Wukk0UhDkLI1x53tzfW13uPvbt27/GL
bw829m/tTePvHLBNWo8/26SLNxVQulekiPUqneJckmilDBucJYynW2j4O4tOc/ytK3j0EcQi/pt0
Nfx0w9+BbtkM2qCIsDcMW/0ks43nKU2LE/lUGMUEiU2clAV8U7ZnQdlliUWBv3Pg7wAEGuh5/A2X
ESCKEr8WVYDn0yV+LrgIfyc5B7+LbGnA37oIolcRO1KEv0uJ55lcrZYtpRH+juPSYF3IfAnlRVkA
+NsW4m9npCIgUBQzZP6iW02esNa6SDeEE2HCgUEqTIEhT9JEZNU+410JD5VS6LzKU0qzPVxqZOca
CLewISFgreYAP5JDDQNP4RUhe1KU+EkzykrqA9IFghOfmSLdC9HR/Z5yIJ+kwO8sIOjfcLkh4Eiq
IDEd8vQihgvczkpFC5DDEZFN33xhJBA+jDfYL6EUU6sBf9NiZxKTnOM9apHhN5AlBSOnKDqW6cdD
illLUiV+kqDFASbQbbF4HrFI6AmaNWs9K04CWmStCWJY+rtqtIrUBkEWeNPQ8fE7GlKRNhAkJDJi
aEzJG/hbaX2KhK6QYLWEw9GBKZkT3PFDeV+e/E64Sy0R72MR8oayDnKV7/0F4BZrouH43l1vre/f
Wt9be/TpZwfr6xtjZ92KpSywMw/682BsGR43TTZEK1+uvsdvRm6P6RSpTowbEksrQoxZDoj3JAgg
y9c5omDBknPgbXmSFUiE/fr+ByALaDnJKVIsMSMEdyvxwHU8eY0KaUL51NlmxhL1ktw2R0ltCSzn
o5Jig2CkY3RqcGaKEooi/lSaqtswo8nOMFSqcAb356dkUc9FbJq07hz1M44/UhyYKRMUwo5OUo8s
hJ0h+N3Yu6hyNTkAuo19J6mLWRDE1xCWjahWyBgyd5iFsEOC5HEIooSaJtN0t+slVDnSL0mQqcTI
FcELIJm8NkBHEkRRGqEklC9ywEuhUO9wg/5k6Z2ihAYmnrnNMPfcIYMR4wdc9K87C4IHHh53gw0w
t+EiyNI1+dxRF885y9gFGanZ8L4xgGB+pqaI7/I97jAS/ZdkvGW4Yaht9UjY2i/WcTtcoYipAwIw
RbDSsk+qeZT9iVEH6DuGFcj0e4FzDUlNA5ZLmx1/rzw5A8WcmI4mLqtsStmhFlYqKugtwtbSZUuk
90nrLu5ApUklDoaEYsccrVJlPJlpZ/0SzotGK/Lw7NOBbQtby4bFu5mWcEQ3MiXRW56XRxl5yuFx
NPrSaThBa5qRGmxOyE15v3t2d3A9tt7Z39rd39jZXEaW7O1t7GqPqO29vf3mztc7kVOyF6HOprc7
EGGQkNtN0MsjSCcNE668qDAG3UOcrfgq8h4TwKrnmt5VrzKhImGbyJgEY2IaqDmQA6+mp00dUgk1
NCltBblakHM5slh0BOHQkiKvaDJlkXrKV7oOuC8H8adJGHHT4Gr5ciwd843aBBIaQ9/yBwvFDGYV
coacHPKJoeZJVLVQnBLcUhH6MxvSE/4jl0zKiLU4BJLqzccmLZED6cplUCpJq+ND4CCsO89sMJV4
Af5zNId2KEDjKDMHQuc2Oim9JrI2uFkHMaYqvDQIio1Nqi110k8JNEk1NkZBn45qSrWkUW+tzq0v
c5WkKP/GXHJ4FAjTNpZcnRXcoaCZ8G4hiQJC+1wmqlLI+7Ipa0I17gRmMrRLyipdbNi1VS60IQQx
elV/C4w4Fpy5xQ1bCUYYnBjQwCkx6/C5cvqaVyxI9ASJQ9aMymD8CkQh3XQpJQIBorDuVolpEVCD
LCuZtygQ6g13IngUwmXWkhEzFJjJE1vnEIvU98B3cNZj4TLB2wS2KmKrOpmVqUgCysQCBjNEJj8h
SqmI8fIiMdZkCQ205KDB5dpwBH4LrGspj+zxCjurmrAStYAB/R39QxRSApUlNXCytWIQhr27fe/j
9fXt9Xm6e31rc1109+bO3vqtwYWL5xNDmuPCOchM1TOUWGIYDcewqVgmPSP+k3sOt56ZDxvfOkl5
F4wziRmRHgFcZRYiYdSmUDOeoGOeIpEPpYFejsB3w96J1eeJaKchLhE13PEk0rlC+ZUoswb7Lchp
VyAzmSGgRs8goKuWUFDk0AJ1Q8Z1KnsNgpveErQbbdkIDWsu0we90zwV5R0hyIDOhb+Tks2MswCZ
pKT9l+DKwKwtGwqG2zlrZXIkfHK4HjzgiELzaApGYj2TWupuCUhK8VQXmEcCkC5YNOzBSlDvl1qM
EImB88njQhGwS9FESBE1ZMFcWh1glwJ+412QkAQZctLVuZHgEpI48GqejXaxFGy3wpJQqcmhlaCc
GEIHW1k5rWCBoVkmeTHNUCohY5ZY6wSPXKSoivWBfg4SBi6qGYw8sv/X3q4SlqDEwEwFcnytpK6V
sOCiOqleKooO9GvBHLB8kSYKKsoASTqgZUmBvfy2wlTBWMp4QUtNkMwldiEBQQdrV5S2inZ11EN6
O1DGZk0SYTZz6Om6F00wwTagl0o76hh3oXhAU7H+ReA2OC6GKgw4VFQ0JNUKmFjVEPxScIUs/iKO
d6BaO1ZJ5Jfeo9Td3tvd2tnb297ZXgbHd/qvdw420C73tjfXd17Nm9v/aG+eL+PN82/Om+efz5uX
n9Wbp5/Rm6c5b55/Hm+efhPePP3dvLmjdQ3O2Junv4M3t4u9efnv7823dm9tr99as2+PTg7W9zbG
rtroNYoWonyhcqYU5fFebxRVeEn1x9mfyUsqMRD1bO23zOCWltW3nJ6RN23KRKVFETeHOOOnFmUc
jTlmtRKBVU+qwS8kXG6+Jk2YsohI0SXooxDjeUqQHZD5lBbh6gu3SlYMxak6G70ZVZ1UZp2luKZG
wlrzvFC0Sbo8ewcOakyBv8HXZRHg0arsJUexZpek1OpksUZZUkrREqELXNYw/ezF5VEelVUZJ25S
nSS1p4oImsL4NhtPypOzmwRKNYTfC3xqIVmoTMCnvqONFDywmGB6j1dmCTL07h3UJza/8WAmgruo
Lgb1w0RbqwOSRkMgEjuHgQ4qiQianMT7immWSFNUB8c8u1lBuUNLB6YghKHbNdE7yS5y7sRkoL6N
yjhhfgl1KaiNLL02SegNsSBEFVNKF3pDgU/JWmUwjxJZOpZcljGSmjuRJc1E6lXh2asCmGCqeANl
RgzQX45auLhEAA/wxClvOiMF+XvHYQ8IUznSoYQWHrw8pAVEUox540IdLFKUnaWyqfLL+61SBK5p
Gh4vJ2IJIOWI9209NPQHAQAGA+2xEs8YsxYQS7CJ1azQ0HYUs7OpwJnHn3DgpoUb8jDUbwrWNsGE
CR3hbjp3se45w8XWl9ts7N3aWt9de5jubG3sTmFqSHxKJt2FibVycT1if2kuCFjjEyj+MocMEvIG
0jrKuhHo/C7lqLmh6VfbGGjO+51gvHAX8heq5ska/+HVNjg7k5aHLy8WAvYgNSAquLF77X6S/3ZA
n8T+EZdSeGnFQVNIdzvtAh5p6eMK0XizKliy7NWIn0nlWeONtbmsvHSsDGUcKktEdxCdlnl4yR1q
ho1GTTURevIvOHklQkmfl+FdZKjoiK6piTwIaAm1qeSWkzCk7Fws9j+h1eUi7OW9yeF7RSSFD4F2
tZCnkijIwT+bchliuhe4UGVAzZI/TpO4HxWUmNpWWDNTGTA9DrDqFuNPyGaYXOKQ/CTzVRjygv6X
qfoMilf7Kfn5od+pehT1FejlHZJswJh7NOiy2MQkzT5HzRJ1GSHu0yNgY0t0uuBTcUMl6gSQclJk
i+Glk0IJayY21sZLa7bDS0sOSTwZoAcYqIvysRi92mY0jqQGg5eOyEZdlsp+9Gob/Ky6l7LwvjXM
AiVPdid0593rSGqzjCRcNKmhEJsBNcScW9TZsD8aYrHJOzL6KZKPcNGetgj2S63h1TaYDbNlkojI
KhDHIoUzS1EXU0hCVEmSMYEiDlGxMvyH0C9ygGC+ZZEK+I2k1d7dfrC59lnefnue5t7e3Nhfl1K5
ub63sz3lxJOkLwe2CSyWJRmnMJaCxIhsKJfTmdIErUVGHd+aIlcgRZPLNK1dFydKkCKhdpqyu4TD
SASL0spiCxdKvHU9uoG+jVpK5KOGiIS94OJSUmVAYndLE5qV0eANCKFJpQhSpUyJWSY+8NQ+9Wdq
tV3NhnIVNPki4KDglHnTvJXQUhBxMHVCcZSrxaC59ZJ8UBGu2iigGz+bwo4UrFZsYhGKM2JjFDGC
a2XuBLRAoWFr8EkJ1U9FGJLlUtAjIsQoMWCM3KWRMzXsyQPZRAyLTJ4eG+qykLaQAcUMKYgRehZ5
Fasrg1A5I+oTCP66BYiREmMQWJ5QQpWagAmlEJSo0R4ZRSevn7R+S465gSRjbMvw9IffF4WGSQOW
VTikjJrgUUTCJ+IibFAZSimypioMQFAqNdQqA6UXUXUSdxSVOklVJ1lDpNxZMh1wJliVopCv50YQ
UIciD9PXImrIXiQSpyYw+iN7gslMCXbaVpVyiDg5aQVBSwYxohlVfwpxPYw0WzBZWV7Vs2ov/3D3
IgYyB9GUdJMmyioYM00yyWCzGA2a11NCigHbWa5FhBRhwCio6ALc4IP80ASlGYSricsSqQQy9186
k2NimKkhQKnc2Zl/H9n+7ubB5i3I7a3t3fUpwkRlR6IXFK9UNqKCCVPO2zRh71nHnWRlfklrHEyt
akMMQQ7WBKmchcNHUuGoVuK0BDsWSTyYQYvdG4xn4iGyPNEysYC4HQyERUB0SjnQIhtyCNKJciy6
M+w9KGWoOdgcg5/Q0jYlbAJy4siF9nWyphC61zvfZVlCifJzLSQsVPUiiRITGFSarQ11bUUmlTqC
LC0HgCkRFaIsgIzeCjv9+iiBHZPwFh6GYiclrNyqysSw9M7EQgIo/TaTST+jqCcroaSMJLfokfTG
VlM0KHJhJhaEHjKKJl4iUrGaC+chN6aHMWphMPgujeKQoiTRlmJyctgJg0yIC/bA2gVIfPyELop3
BMXD9lMwPiDoZgjtmHoURQggeBj4KDMJVw+oEWOujgm0Ak8IQ6ZUIpafNeY2S0DIMX3wmJ4ac2NZ
0CpKEeke8296iIWAALPFySwDrdesaZBt5iT7Jl/NrRLRYaX4jLhzaxELurxoweIlMgbKqqrBg6QB
ol/GYfg4NQpnSU4k0g98bJ3s2SdQKXRLCoP2FyecRSGZaf2QoIqXNoWAjNdTtpBGg4Zox+AV9kn3
rwryinz+uBdZGVzjoUW5GMXThcRIUkQtHLxV1ce1e/bxxs6Y2U4oNIBMXGSCasGvqeAly816NORy
vrFoFtjWrSZ5ffEBfMwCxdWynNzvgExzCl6PxDUkD/J9koBgyGDM6iWIkKTzeJisS3v/k9IvyTo4
Aj3b6Jtsss655cLQZUEDU4Vn4aIiNUREqk6tpW04Y6UnXCEB+2ApYQUz5LMzntQZg4zgKEpWrE3N
xVrwYFkgod6onk7RZJR2EZ1u5OEBumg9oQI7+5dDmiFqZYEfcUONeRRsYMJwrp0SYxJWwkXTV4t4
F/U2hKAQpynyLsqDkFtyy4vweDCtKQdyAZg5dJQRYJ+CcNCyCYmKxaGUNup38S0dsWg8pKN+kuYC
zSzmrOUn8URODv4ciIirtiTFDOpTfbHQJqgcxj0m9UG4EUU6ulCkjJJkFC3zDlibw/JAkKb0i4nB
tpPSrRwJhAuvJKii1ky5g4lAkcdint2TYxRuZZ28MWuLKJHEy3FxnFvRdDJgkhg6UrkYfTqEMwF5
NGaWQJ5Tk8Xbi/KGomw9u6Ia2xHRLO5bFlqI5iZEJ12lmYesBc4dCgohDcuu09EFSiapD9FH3EjM
j4qKExS1ACZvfDKZrHCeaZEq3iR5ZpNK7Z4AGNeFd2ldrAZfK7x39ja2d5a+BXhzfTMKSrZ2tzb3
fy8oGQobwmNeVFByqfLA315Bia/gVyooyf+QBSVJR4mvIpqn+YKShGm7nwgSPV2loKT8fAUlufwd
C0qS/V5Q8isWlPx5Y+2zvLG/iP3e2N472HMvvrW7ub1+1cd06EYp5+eHXx8enn8NlvkNPabDomFa
6RAZPeG1RVTPFH72x3RssHn7/TGdDLUbDm/6MR0QbG41QKbBsZhgGAytqd8f0zHKSrAHOIiiGcIQ
ihDE74/p4LhMtqDgm0ronaaJFa1QRA/9Jh/T2dmff9Byc29rJ/TLza2trf8GODz/Q+Dw9HfH4fl3
HD6Fw8slcXj5HYf/jsP/Djj81q3dva1ba/efHj872NjfGj1RKcTGAkXMxp8Kx7aKRWs8vfvWFFKj
uKVswqRZ9GWS68i6CeSeLDqcl/jgWhwn4bayCkBQIRyTZIqWiswdt4m3pXk6ZEH359DcJC1JHSIC
ONvpZbcKqxRbGX6+mLC2GHUHhL7+jGTb2Cg4ULEJLCYLIFViMlNbo1Rs0VNeRKblTcBTT0D4poLT
EuBbsFW6RriRjKQS01JMXVXkttBhsrSoNLxpXKowmjT3SizMSl1gmsMqcF6+BaAIQEensRqQGsDN
8tvqEivR77mkpjG6xgfuTAG7UqwD8bEl4nk9qAvJCZ9SNFRFCV2GrxVYjlEzgURVwgMMuVAolajl
puH3ECREmlCccieHimYuaN3mvIFzC6TFlxqALqtXKiIrMkZJAIlLJqUweUAk9RhK/qD6w6lQwYQr
ILMVjT0MB2V/XU4hSSalrkxOigRMun9oRqoa9B75i8hYuC3XFV2dg/5LIRwliUAaYPRoEylvQyKS
svLnotFrcl7kuZ6SofJoSZO4ULFQ5FByihyN1YK913Y7cCbqJdVZob4ybaTllDD5FcMrkySY4Ipo
Br//pgmLbhf0bvorwl5XLCCZwKocF/kf1fPWhMGi1BeH18ldhOIg+ntvd3vn1q39naW7wW7curWz
c7AF7t7e/f15yleA3fm/J/39O+x+Nfr7V4Pd5XfY/dNht10BdpffNOze39/f3Vpfs8/+v/cPNm6N
XbXpmBTXSkFbE56SihYJsOBJPfWQtF9djoHzuFRkaf4/PTRQAhVmVVC45/ElQTkBQQfg46Mjwr3w
sJaIXC3iylMjGkMH+qrjaUjYUNWeFM+HcW4lJblwd5hFvVJ9Bhxkho8rbWES8vkCX+ueR4UXPrsQ
w8j6nnlwYFHViZQD/yTz8BaPepSsmoOkx8NKyyFk8EIrdDBzsvyfkIxYyKw1GFUnuVWdwIkWwG8U
Ayl+tN4Go1piZTIiJUWtCFUnFo+lxbAlhq3I4YnG81bgXIP8hJS0HA9SqfCDS4lkHUZBoNKbzWm8
8ZvfPGSuiiEMeSG8ec5STQqTw4vsxCGK9y/iVhP5ha/DMDuDe05iTrX2o2KUElqZL5eXc2B1EyA4
kF0ghQ6LnkTI5CxZNSbcuSp+hkIRB4JJ7bvTzNFDqRfeEb8Xemc8XeWQmGc/yAApcPFR0qwmwSrU
nzK4WWKbzFNVa5ll7Isxy9GoQFR9Nh5vIAikjoWSNaygPZ94rJKFBaQoiiBJpYz1cawYJ27SV3nh
kFhJuBTs36LwJBeH3ego0F/gp6zamRyIPpFQpjYLiXiO/lKgsCOtwKXgiIkEuB+LEOyrvCOiJaox
AdDuMoxyG7/jJNSY5VL5Mvt74oFKykPrRet73ntPfWvt/rPD5wcbO/Euys39vQ58kAaGJNLcHDKN
PJ//BtNfVKXtlf/BkCQSB1NOGYlGsAPUyfiNITN0KTQc0sVgSAS7iw3fZCE09CobMSRlxJAEQAB6
TzEkHmc0kln7fdvAkChJlVxDlzhZ2W1OgceGNwenEt2cZkhIv4IhSVJlgiHR8FFOGmxUDoYkomPL
9bNynfq7XlxFS1FcThEu3qyMsvWsdJjc118fpWnEZyBINYbE0EZFFwwMCZzRiCEhBUrBkAAx0Kca
PaL1pBBKIpjDDkpMYohRBUxXRD4U3TgXdJGLqdUI07BpMTR6hNpFzaPUZld/lQiqTG5EjzR2KkE+
JKGjhKzZKUGGHrGYcGH3SPU14Fle3qACfY9zlSALwkjGJAKnRo+UWXok+7YOcJEjeoQsWuIX+lrQ
I5A3cUwXI4S/DsNLrURmmh5JIbf48u6gNHMs/yzuM8vrtDUU1wixvXZqUPmjw1BVQYyJ3gpWicQq
yRQ6gPdgIWmaHkkDPZKn6ZHiq4rk12O2lmHQI1FMIXpEoZcbTJ4lY+4idSX/m2a/BI8jMsuC5WVt
dDZOLUX1Km9VT4HL0+wIwaMjCorL8Rup1YGbm3s72xvL6JFb2xu7Kg/c2drYHjlybjioajnGiMdZ
DAZ14y0jN/3YjX4uQRIl4YvIBHEXSssiS0lJLxWGhIb+S+F7FEJASiQBSQhNUKErQx5QLCiR1Mjg
EN+VmA2ZjjMcMg1gjIiqkgO4UJmZgw9xmMjxpieAKdbMwwOUisfgT2EpdyhImTosKXmlgjeDfGT2
ReWzKnigIBToVSL+6BeloZH/BEgvjQ8RHRB0iINUdZuYKFoIeK4Y7rc73BqpOpgH8jUrgYkacQgM
qr41TupIlnNmFNqWHEL/Ai0lMDfoih+KamgE4FNQ3dSKFAFWTaZnIiYuRcqEY83I1kgJC6UWbsMi
XQCFzqTECIHsvMNOG7ltk3qBx1U9E50kZwiAWbgACM/YVDumh4fFKBcSe2MBh3OKan8yhsKVs4KH
xkOjmkYWog8zhJXCeG3Fk1GZRI5LsULEq4R5hUCLb3aHxM5AGgA9ggBDQ4vu7CBNdfOoFW6QncxA
foLVI4LE720gjJQfRrVI9o0IovSHOxIqVotZaD+1YhuRob56edt2zoMBRRe4kVyUBhYxJ449WzKq
iSfsyVH7OhlxI0LOCAWZ1YJtkw9r1eYmQPntEC8gwaBQdL5/2bGa4zjz53F2bm3u7G7v7y5969nm
/vp6FHrvbu7s7IxcOT4Wwy5ZlBTpfQYP0sEyHCZ20p9BCWay3tLXX58/fnz+9ePD86/7v/z/h/Wz
rynHEXFK7PS39SGKaEHB2sEMYD6+Vsk0QXCCjVEjOOJTxCWJhGtrj1WqUqOi5KrL8HmRdlPkJ0a8
pcslKNeswiQu2YmesdYklBDcBI4TNtCjTA6moSiJUZJf9NCNQC7dcaNCVBQVScbql0ydXLH4D1xl
Dl8n+hf6E9u0MC4/WQOWcarB+5q4GV9wFsPpaWkWC8GDDXJKAQcjcYZXxkspaiiOKhPn5fIFup6V
HfevcGZR78UoxhC6XZBighWyCCHNDdcUj9n8Yg6MmzM0WWNqBu4WrxudFKnCD6FDlHgGJZjx1uEU
Snp8OBhW2DJRl5nMw5N2RrTP6Ebu6cg3AqlireJeiUlaD2UYIQNCiBtLFr1LzA7mieMsOYKdFSEa
sE4sWixEtHYePddP7zFysjUrkdRoqBnIHCvId9EqJTgPlmSWWlJihRVBWG8/yU5cfIK/zOKTisXd
ShGJMt0m4fhBKYkmL0mYmgxZKZtyO2lm4psapMQP1ZCdA+pnSTiQAgOkTEEsx7oHiqiQwHIUdogL
UuqmU1Osu6ZkiaDsBOWKxRx3t3tHvbW9vpbSi+cHm7fWdwc/negcWlqJOJulysLIlACXJOo5wh6P
4iv6meQxCjyDFKJNE+vrvZRw1TVT1ZoL9+rMowoScvgpuGcT+1X8cQn/XqdB3WLI+LzSnHRBe4qD
eQpeGUUuYkjGTjoFleVOOtywPEe4sxknXZRvlghdooPx1oRo1XIX2A64JS3WrAy2rTF0vwBndR46
AnjJYVAljGEATqAEbr6ICs1SbawRHrD9qqgRMSd0T8ajnD6S6K6UQHoEDZOFReiF70yRJIjGJ2H1
F4lGCs6Q+JKU/0vhOWERxFsru+JxP1wBN59UFe2XkVcewouYctltp/EJgo7fFPokpRCQWjwsFhRn
Jx62eWULURavnFtgJ6JZbhaqlx1FOTa2XYKsoFqhEaAF18C68zvq1CBuTPSwSfYOrJIbhomkEmDY
Bb4IH0y2yqoBlotgIKssktcdm8Q9a22lwWVnKR9CUoBlTxXDxRU2LGR5Q+Uwb4EHPB8GTRQLy0kD
kcrQNcdOdpAi+RPtULRofRWIT3WCHYIMXcsCyBaI36zShCzGMSTCFKMN5ePzqYzbeW0PD7gciCxY
l4JGS1UFkRBHl0OFL3oCgZ9zDm9YCttGwL51WXy+jKjuGbi3u7W/vrGxtRRyb+3trm8F5N7e3Nra
fRXInS6E3JFAaHFJkIMjvQByB633O+S+KuTOPxlyP65/zuv/+yyp/69mSn2e1P946P9fCLnTLw25
0y8FudNKyE1w/cUgt/3dIHe5EuTO/9Mhd7k65F6MvhvgAn2nhr67Efxm7fLeha3d/f3tZQ/J7+zv
r2+rSnB3c2N9cOQ4kZTBOhFNxD9LrWSdqPZVAbYOQydSOpRxUeWYkI8cBRCKtyHkAC/YN4jlghxh
Ufon+hqoBHPsKI+a2kTewXiZtD9UTZ8aMFDE+FZxkmBqaxSQ9IALTWIJkviCCAXCydGz8cM5NGxJ
VGRRickcKw4lGBg+USVoWec7fDGF5USoMtQrqjxyVAnCiqskIQlrU/Xga6Kx4rnpd1mak/ygwfZJ
r5W3zYJGYZwOICkM8Hij1IEyLHwWqZvsXbcBj6JkC10TqGsCekxKFvumOccdErlixgNMCCsU9mnJ
SI+BamKVq2IvhX0iWsmDeWsdIdosGlX8lhdTitKux5T6xSCj4S7l3FLYuExBkuwAvUmjuINOhQxo
LaZSONVkFarr4sqyDaJ1UQjxNNawlSJwh9CahhLBxmcqLaprphNsMEUOgGcYSVIqWJqNJa5EdPBy
Oz1IgRmVmLasaoWmYuRIdyPEdogkJsaclhnYAFgKuSU0L7ckX77U3iq3b/3RlbGeBtqRgVMkl7wV
PAnlST5tZkBdomnkQlxFX/SwMtDgMB5F4qjKtx3jM54qj0zS7WvDXQDYEiAVY+NqDL1qBYpU8MBC
5lW/OliY0cu9t7f3trZvLXPkG7ubu7Hfyc7O7t76pBanbG6tdylmncXghqdFE7loGip/5ZAEprsL
0XQaoenQCR1u+UtFr4ymIYGKx1EaXIymbYSmtYSV8mRlNSWEPjEYIzRdFqFp9d+Vf1AjPlVoWppl
yGxFmXzW+lSxU9az7llQXmi6DL0Rmk5TaFpguhujaYvEXAS2RbidJbAFpju/K6OsIuVAQyM03fz7
CE0LTHupxAhNEyAscKgbhNC0DQQ2KI8ywMjP/PssZBl3I4gXaNoounJD7MJfLCWwA2WMCWxfPBSa
ZhtGfwZNpwElyoSULnnQ6lTNP01gBwcooNzQtAlNM6LIc4vQdBmh6TKLpoGbpVNEy0iIA5pWAegc
eaIl4H6mUzS9NJo2ofzaX8LmBWg6z6JpLQ3tLLAMTcsljMgTCtUA090roWmB6e6VCGwh6G4Vmr6I
y+7SCE1v7Wxs7K6vfZrv/u/e2W4GVt7ugqOQuEwJAfUiqgQRBjbE64aiq4tlbZUm8qsYAhhEEFZJ
GM2DADMYyrSELSuTLuFjgjm1rFivqi1l1llbQ4MDZaMlikv9SsERZ+GpogzVJ6sjtYbGCNpUrYRi
3WoZWgSR9XTwshYeOQfRJeBXKHAJbsOEJB0NOx2Ls0qS5VPLUMioSQ9YMh6s6JIZhH8cyKoXHRT1
UkXURxDBBrDyPL8TMeJeVIsy6akNGzwV00BgSLG4WnV5EQ2q2kLhvsZzp6KfwAy6Rid3CDXUHkRh
SC2goAgw8QIydwcjyXGtkc8U1lEehaBkARxHsKro+SUAuAJ2CZ440F7TOXALtCOQ1SXljAVHRvCg
K75ecqMeW0S1GLZuEBSKEJimWk9oyvladM2CnU9Z255DZUn+oacJ3J3xmHowzy9MgBDczcLlEDJ5
RKQPgo/chUVckKQkxrcEu6hwhtMx7Gdw88yl8JN12uEG03X+iDIjIkXwFIq8OeyEKupOuYW8RJaM
lXAXjIPpcZUivo20ufijZgIMirOyMB+DFtpx0BBsRHp5EpktmDUsTHUloibdRUeYp+LXJ7czsjYL
VkGDmzRVRqFq8PzGCtBUdxRtU4BTX/y+t7+zt7+39vHh4//4YOtgY3Pz1vrgp0lbAq4DmJQzeI/p
d4GwV/ZdC5i8aQ+cSq8VHYVOC0yJVMqktVckapqCNrlnijI7Ja8qMNayzLFOs/IoaLQkKDfU4hHD
k1LvHHx7UU6SwgnxsINxSzhc5ixctXL5rKiF0gCJ5Ou3WWjz9xLIRHKWEkRtLsIxOr9jDYAGUiMh
xM+SlsPDiOYIdSo7p1F0nHyuaAwdLaKYpRHqli8F159yfCFoj7E2XU+1dwmMW0LIxN93kYLo+RmR
l1EPWYZJVXqXo2Qre5kbECVb8OW4e7ShMDIWv2RnC4qaAeO6vnJyfJBjztRrRpi0nx62bvvax8uI
dPKhz6kRwi3kRy/Y9ilhzboghlVKaIMBUVSgKHXd3USH8ZMdgnxjeoTli5JsUgi1RQrcwYPkEkVu
2HMZHqHJTX7IYj1ElPgODY010oxxKHWqzHwJNiA3bT6D4ElgW0QsRkdy8JwKg+TirLwSbieXcKsW
iVgRiIH7keQnEw4l28e4a/67RPIjf2MNK7VgrTuC4RaAIYFVA0qQQzdL8kdF0MsHxTQpIQk72xrw
RiQ7hod1kgY43yf3mn1LsyyQiA9sRytxyEGIRZspOs6joqGXOaGxs727ube19uc7dz+8++jBJwfr
G7f29saeWlXjhTJ/2JX4yQMADzHquCQ10TNhGF1Pbw0yVsSMoglfQPdMEx+8esqslZWCdil49Eev
LJ4JYR3EU2YkE+yxxwoLrl7RhpQsKZdg1uO5HZBhVxrsNSVEFkkAZFgrnQ7tS+wXkAVTTcr8mAJA
pKJnFnEJXwDr4uigy6aSEO4UW9XXCtWC/ymypUIqm3jKgnBioTJxkEg8LLBQRg/My7p3zyJ4pIWa
kUAcUtEYhpz0JBTuX5SbZ06kPXh1n57kLZuSHdAgDxrFBo0J749Q4OORZESmrsbjB6ga7V64bZ98
NjfXkomGwNtU1ChfF8UYVsozC13gIyUlEhBlrS3F9cQ2BQ+fNBo8NQn90Yp6UFoCUIlOSzIhApaj
u67pi5pMcC2P3AANSxSI4BVT2IPpgRj6nIKfQQrIlD8nrbHxKJDf19KKMkAd4h+3miw2P2GiQ5gP
2t0BYglpJQf3zeD4rGI47TEz8By6q9dlQCZD1RmoX88ECcWxSGR2saT8Zjs5BxcT4Ll5nkJIXEs/
LBfLi4fUuyjJKMKG4YPgB2OGlPPAN8jbmNP1PLljaMzyb4Bh6nfCB2p9WLg93sigZ2xipybCW5EP
DPdXSPLkQn1Ou5DaYOK721t727e299b+9fHzg42tvVuDj5byG4mqsuJsoQYDCZCtVB+UoeBTYh8R
YrqPizxkUnIHcchSVXUAclqdSyrYsoWGF09LBSMkElN5ZhLRqdSkC36miHMr4gosgAV8DZFUUFpM
jOrIgk4uSVDL1044L+7Z+y9OLQFOPHRL/VC+TWDIIjGo33BDgptG1HO31UVUZziLEJz7xRTrHWzr
C7KwDpjcLmnoMFOAXeYW+K+E/E2TQ5bhD8ok9VYZOGmg266knqR6g5QaQ+AL29PK1EAquWYS3U5a
oAUHkKN3WlC8jAGvgkyA8BKsN/Epq1KvJbOaWG1CIrFAiCepWIMxgyJJGj+4mUji/WasKJ7EwBap
EjIbFdwoeW3Lzw0QTjRJPzQ6LDNj1i0CJ3mFQ2BefGOq4iLbIuWKXgAQQBBZyzpy59CgvaM+woGj
xeEz1cAD1DSyC+ceEP0k0iWF/8gwdDfIYMJKzIWHHcALzj4lSTkRxQrkGASKxBXAsBLcGGqyck7m
6piSmzFEXmOA/JNKLHOfkcU6oeNzmyKbKeKKSi6hnUT603noDkJOtl3afGLvSYJWUqTSzKTOdGWt
BOoDi8ISbCbMasQVBF6jjCnD7Qb4T2IdNE3WKKiGWVUi4E8KktngxHN3e3d/fWN7zf51fWcoq97u
pG5Qw2B0qEQVI/x0EnsuhaMEw+Hbs0h8T0X0g0ktkGOTH42lGEyA+ROZjeuQGC0tIxy0MjIlJzH0
LRl11IHai3fWuOgrOIUirGEACfcEnRJZYFGkhKCP0syd/EaPBEvvj26LmSD5Bi9aXByqVby1KA8k
g4KDziKvpI0nMR5ys8S/UIrScKXUNRLRlDIgPyiBTJHPlaxokk0MUB6/USQr4xV9XbR8w3ZNizSX
+CGpLqA03UL1MeJJFJfJmOWGlXzUkzrGXbm/Fm6hylW5pTxnQi+X/wRZSHDEblJDDFmFC1QNyGck
hHfxh1SHZUVHLNa0qDVt1sBkjsUu36GHpOXzWavAEtGb4fBS/FjAJL5sfNO+EnX8/EdHRKc1H80o
lmB/fBY60cxZsi3dj4wK6ALSQx8JJO5TI8bCdLXmneNMWZR40tDBFNQ6Aryag2dR3FVq1czevTML
QOqlF1JoeBGXg80rAblDZuHkULYz3S6Ruyke5oA53HbJCm9t5fBjhhAX3REcbUqtlKOM9EhKu8Na
fUTZ6EW3XcR4s4BkYcTvSNbh/8E/bBmXZcsmFhq6B30ti6I3lhuzJDDXmZQGGML6WrHd7f29tYcf
rfdwevDQKZIQhczAuaY8jDjLNKgyS4HF4t3ePswKLdFPwm0KE2ZMdRk/3V1dIhoVE9trSmwahIbh
QzRsHro4nFQcAyeGLyIXG7n98NCQu04p+eamkSxoyNDecxiI7HHkoRuE7BKpPHV9gk0CKzYQquQ6
AAw5iYyWG57LZAlAaMH35DFPHhoQltFeeVxc1BHEGmPma8FEK6u6hQWC+bCkwN/k6UrIgKCxhCOS
B3SRIOK5pDx0mfbQIjLCQ8sNaspRr93VFcnLzgdB45Cqx/SAkTS9YrxZdClkHUJa1q2aUrwcBiU2
OodIh+jViXX0GQOFhQvJgS7pmgi+ouqXanmd9gEJABa0TfwAfxYeOoWHdlt2D004kg3j2TPaOgg1
tSsr9pGbe+lqCjlUMJLcSgAP0yRbDusrGnw8NIWOIpGSYD/JS1H+Z+ERLfBzGnvoGAxLIw+dwkOH
M4AxD/duYQ2KRJ5iBUshKSIkUwOVNOJ08NC678DPSZm5g7AkEZhEC8Q8eOgUq1/8RgOY8gSp4Wes
NTUPbSMPjdIrv4GHzmRKptw+iu7dSCxiKAQcFFmJBADeKjdV2RRvPTZ1GT+nFrrbW/sbu5tba3bn
03sPD9bXtzZ3Rmx0uPgkLSRplSmbTSkImmSSOwBZdRGlIBqpa1BTUpES9ghQZXQ9iRF5z/bCDJAF
GaVkmMgWYmCK+8niHJ2xM1Y4XmcI1UWJSA6AABSI6atj1+lL/zgHIyn/iBnB5ODsgSqRVHZBHLfi
AsK3ECXgRwsgA6sEfLLvuaC1qNnX6nU8QoJLHyOnz2EcpnUQSyS3TDq31MXHoS3n+Der/jsFMNPY
SmJotKDqRkwNiTF2hN5hWSlIxUzaBdMumCBooQWgTL72gVcwIg1kALhBM2HQsFdSMBDqQ8zgCVgW
vkRVlmyiaNnEhgOMLChzVQzVkxVSWDtSxoOakIaJpTqfllvk9ZpD8hrx21kGTPgjjw3qRs7XLMBD
Z0EBiF5vkTFBa2bJFNCEQTICDGJ7JhqO+4lpK5rsAGbIaKKCcBxZ9ic6S5UfFvSCfAxMnVY7hR1B
exIFifEmcVHJJFCcOBJSAtPmTAeuXMdhooaUkaLAoQRjBG8B8GbHVBJujFvkrYiopOmBnA4jK5DL
gB9LMp4UC16jgNnmHHQaKZU4au2YqLgFdvHrZOlXreZKhQluLbLJIgfYaOxIVsUvKD+kHdaWt4fX
6MQkKfHrbm/ubd7a3l57mD7Y2NqZFgz9lpSM5Za9CaYkhRZleGQ8sHtJ+1KnHHmhiFcTrZZEMMpP
+EdZxK4vXS3LEi5aX2pKmxSspWQmvoUnooGVSInKj/EUouMoNAeHpxFGznrNtVx0VnwU3IJFFyZN
QVdIsMpspy0nklV/oB6UcNFAiBR2FMin+nfdsmlEBhedoijQ8Sfqjly07j3pzW0QNNFtmF9olUAr
QH94KCjBwtP+hA/psUlYVxOGo7TIaZXjIQh0WfSASCk3VfHuWVSx5RhKiTUQD4OLxsvmIrcCwPDl
m3O4aPK4KRct6UZKvOhbYh92Yc29ptJctBcfdXEwkUjhHKtK0kDoeTHRuSI6areVOBKbyQKy5bBv
U57DFEVPRJEaL7rNory9PVzVlItWLIa/DC2Y53sUPghZuGgpBQXQrSKMKI3Djcg8o8YP7kiSQRTW
CVSU0NO5JsjWt9/LIxdNPAhsq4KuTAIoIIGd5eaiJdEq88lR8DntorG1JLNyF83vOVx0acPu4SRF
dQuSdJFCR+zirdemy8taYhQQZ5uLDj82ctECYt6lYNFyLEIHFEkxVOOAfKOw1KK0zE+AEbVFa9xo
Q1qnOM6i0UZDqyV4G+vbm7f29tbK//PoYH1jd3t9LBkazD1YX2xharFc0DUFhwkjzlR1rEyIC1+B
QUonap1Mq33ESPNMzQJGuvEdAR7h8eAqLJ7u8r9VpJQUk0iGbMRIlxEjnaf5Ds/MMmghgo+sIjYH
nGakM4x0Vo29ln6g6MWMdA5SXq6Ro2GkUy7x9YjvyKpPyCNGWt0yyYBdjmhGPcxCRrrkYH/FAKuQ
i731xEinMlxzISOtmleFx9QFqh/4DpWgzTHSjbik7rCw25sFIw0aLiQbUCgjRlqp4eBOuHIw0sgm
MNIhZJcSNHUR7A1uZMxIlzm+owxRARV6YKRRDDocWxYhDHvovUgBPyLVRDgtSk10ZXXKgpEuJW5e
/E4q04x0In+Ekc5Sm3zxcYfKxxo1HHS5yG1Ibsod85DnpaBlom5PPFHTm3E7eZaRzrLOMuI7WmpP
8WMrwGPyKEsfGOkSjHSWV4o0F0Y6hwTXGGkzIVguDiMdYpDj+TTw/KFABN9BuWAZMdKQfo2RLtN8
R2N5BKSjxsZURJuk1GVRTWJFRUwnwZIOiG+SQSXVh64XCzgL1bOa0DMrZEsRytw3dbdv7W3tbGyt
fXzv4ad3PtrdvrVzsL435jw8rygmniERbJW2ehAQRkpip3IJcsSfwm75yOig9hEJBoQexB2Jqq/2
rqj+Iyrw4Iaj4gUkEi5MCSRjkr0CD38sUiMzNOGMBRrR7vOIY/DjO3HmSkRUzyCShOAn3MCjJ6WE
rJscdHimpSK/HNiBGSPPy0QZC6rcRP3FO4WxH51p6p6AkQTM3BKiAiS2VoEHfZJFhWXBlbhIzM2Q
0Tv8iCdLxM/jPXGBSrnl5pNKkUTxeFepvrc4S0kNxauiHbKClCqWWLUekTrRtSX4U9CcsJraU+QA
fVKA6qbYAX9dyXLKAZAWCXCck5ROES30BKBCslkEHT2Zl0sjjCMnMSEqxtWvA6wtpjEXBVHEgcha
shQptQ7VTDyPcJejFgcoo5oXCyFYt5qQl33wZJ6Sj2OYS6S2OYpXLMWvAcmziuisBH9FuMjRMXfV
QpwUBaJlqLLNy54GHb0FyYDcJUBKkblkkZ3VZiFMTFIYPGYK1TaFnuEdUEVrmBAl3hmU60dKFOb2
mqslSckh77aCsk7ERJLdaslZ1DwRmJvgVpp2DJNokgzF/1CnCtErPkS5cYAmQ5Ezr4kPr5Ygfa3h
kCHVy5I8cIYW4JA9+Q18Xjez3tnc31vfX3v47/njg4399a318ZOHIegrkliYhw3KiwKh1EM0xhqM
ws8BqIVDkQCxJoS/BqgtAHVpgFrRRhMpQK00oyjFwhOVMaBOMh652gYdU7DranQQEMVhFko8LAC1
gI4ExOiRJLgRoBb57W9+z4T7ZKoU8dltAqKbRZFXwOopTVWCFiUezOcYUOcA1GUQEAFQs4B6qsQD
EJgREMm1RawlZW0AagWrZEFxy+kLUKcBUNsIUJsExADUOQB1XgCo9dCJ0gjvUBcqneiOLOHMoxS1
0yVWBM4yokUCUEc0MVUdtBKPohqi0gC1YcsA6kwZXZYO7n1LwZrnaQExh2aWm4BY0qjEI4XWKhQj
HTg3KTEERB9mf8FhCkCdIskh1c8N3ytfVHtlHlCHqjYLqC24MUHdGUDdRKkC1vI7iYAu+UIaEVY0
LvHIoxKPVCSKFEVDJyMSTxGIv8lR1ckTLbAUFmle88YQiHKvAtQyyNqjLgXAb+SCTgxAnQEAAbYw
LlsIqDOrQuytiKGQmJcA6nGJR1JbOQegDhhmSg80B6UjUS6BHKy1bOR/nFhEj4m8EGPWieik1qx0
tze3d/d399cs3ds+2NvcGBd5QPqL3GDhC+dDzIoCL6KCHcG677RumvPABZoKCHDRknEH0iPP+Og0
56NDctFiVuI5Ij0MHy0XOvLRkG1UiMXMYGESTlLz0WPSI6gwPUZelpAeeUx65FC7x6SHRZFHjjK8
kociDxuX4dki0gMfXUakhwjcwUcHrZPnSA98tAXpUcJH49O6pmxcpgwPH90y7MVleI30CAqcnJf5
FpgefHSQMEo0UgTshCgiPs/yqMgjSA+BbmyR4jBme470MCuDj54mPYBJeUR6pPDR45rbhtM6eYRC
yJciEzJoUmxIKiPAe+GjcyvysBaQRqRHRGFcc07ho5WSeH4qmIOPLjM+ukhZxsEBi/CxZSjDE79k
OUiPtIj08N+Q5Kd9NEhZqUZR3pYR3/HMzsLCPetkkQ1BeuQR6ZHnSY+FZXgjQ5aPJrA2H52HMjyw
9mIf7f9O++gy7aNFUIx9dA4fnYviCxZGvpZExzI80O6GuM1dRSAgUInQE/Xk3kDqa/BU1pVWE1iv
3N3eXd/e2tles0d3PzhY31nfH+mHWZ4KTARFrGAixbDESpREWijKcCANKWKjmoQUdFiSGQuhefai
lusHUeJhUeKR5CJwyRLwnWOBEldRnMePjixjKPEQmBNsLq3EQ7X50NukrfVRA+U+eVzikUf6YYqH
ayTzhn44W+JRVNOAbxMCCc8L/Izs3T28rNqCfrGRfuhrKpcmK6USWQXQiFInC2cyVP3QfcNX2KBO
WgNedcAAZkllR2S/0g9xuqr1yDhYgW1iPP+kIv3Q07iGpOIszZJF3U7WKrBWikEGW1SWkE0ahTQ6
S4OMWUT/dRb6FGOQA1MA2Ch9y4yg4fZLSGNe4uFOgjIfzI/Wk+Bcw0ipSUYZOdHf/Zej1iDqYFKs
cOnSUZRpUlFZok0/VHlAfF6kHtmMfqijsD8Ky0IUQJgWdZajxMNUSx2cp2TQ7A/fRfmFT7WgbdH8
5ZkSD9OihWDscsj2mtAU+qGFVSdhm0E/1CpGP8wYR5sQbgVsYUPvExUjJIg+zz6KRBSoDtEfORT4
LEBkfMhImko8lM8N+mFJUeJRMGbph1kEbOiHhn4o28GGE0lctmAwY0VliUaJyya9ER77l6P3YYP3
BDt7Dm6huygvSEBcaEDkR+tub/VOem9j7V5+eLC5tTvy0WIhMymQ2G6qhsAQSvlSVohlgbjL6ojD
zAAknGgd0BZ+hIybu8gh2GS22CYiRBGFgIWkBbm/yAlVIQDeImrj76UhMsslaXyCtSBFZiV4N/xp
xZwahROUR0JkBuxAXSYVvBXiqiOaLgAPQx8Yo4jga7ohDRSoRl8lJryQG7SQYyNZd69fQvEIAS3Y
XhRTwQRyGyk0wnVFZqks2XRAA9IlriOdZCA7sjhXzZO4z2XqYYiDJl/G6Gb5XUk3JcCmT1UywaMc
NMmgnppIcndZwi6wpO4XtSdECfJIZEsWZlWnEjCVHNh08aytpZLI2oAXWqnyPAsExDzo05m/Ilpy
b+4ak+gXkUVZQFoBwnf6N4kPEZRSA/lkKXng1/wOg9vveOI3mHh8bIpsUdJlxrwSDAq/1kmIwvbS
CCeIONUm5NxQeBJaZHR15ZTEoOPdsCcT2aTkLxIG/ywBqoveRkYQIzdMSleKkgkYHIJgGL/hWYdU
LYh89wy5ScE+PggecE9ajw7uOsF5biONBCVlG2FtYuoCvoL/kdmixgjaT841C62K80gSwJMoafaW
ps4i6YscnHOK8JJU5Yb/KSrxcE3Og4Wqf+pD4du9Q95fe3g/bR5sbOxsTu2xFFUXHkySyW2LJLOw
kiQpKkqESjFh4cCwVEdn1aMZKWCBJmlsR1TdzZV4WMtGR2xHgssJtgNQWhojjXv2e2SNUcSaLRpN
JXoFUnfDESOdmnvOSHZyjUVhQnmBFTHSCbYDVseCmrMRI51C2wnZA3AjYjFre0IhUsnOU2wH2GrE
doBrNIuZkhk4UxhZoE0KYk8ESk5ipMFCXuIxZjtSkGaNkU7TbAfAIqSSEjtWWwkDs8YUhcJi+LJG
fpgQTLAdhPxGYrWMMAfbYY2RlmHUS3RZ89fclu5H6SzZjdgDZUXhULKcdDPgspyRTuGhlaGyM2zo
TlOMNB6uBEuUYuBLPNJSoAHdeFTiwf2T+6UWMIJnUVhvjLQF26G0cprtYICV9QZSLBr8YKST2N6h
xCNFkaZbqY9SZIClERZp+qHDSFuDmIDyDCEmR0WhsnXWQVbc8Kkasx0wqFnZf/PQhM4uRYkHQTFH
H3Mw0j44NGotOOYx26GelyhhTe3e4U3EzFFRWebZjsZIU8khCs1KgGhZoUWgNu3Di9GRWTJrqURO
WQTHTTqXPJsjWr05x82xrtDbm5u39uvLyu/8P9u7WyMPzSKjCDxr/gjDSeR5lpGQoRO73Va6YCqC
BsEyoAhIKlJLmgwmVdciihg5E/qHuwNfo+I5SDySfok0zAA5Iseh0gkfWICOdzvgp2zKceRWuyyZ
WkNmmh8QnIFCxXO0sYY84SkHUgZppEUEAmBfETbrfEiOLNTu3farQFJAJhaeymmMVyaylsZjiOLl
PdIpuBmZBiVXGDE3EpcX3wHy7qJGTP5X3JMRRmKaOD+S6iwuwzfsIT31U5AXovzGItUXG2XRP82F
3h3BdyoVkosHCpB70f8iLcqkc/iLjAm5FgwHvY1SX79Pk6WpbAHGMFFfoaUq1KZK+MjB8YXg0BQy
vrghRjJ4J+OstsLcWAUfjMjU7MK6GJlmRhZ2YGGGwbPAN2meDZJTI5yjjzhMyNSI3k1Mpj6DNKCQ
tKSsRlmfkk0bk2UMnmmBSnvzRx+UjToNbVErAIohpESJiA1wE3zVaeUFHSEhRffNQhLPoXmEQ3R7
6MIZN+oQ9smERkTYiZFg8QPLEjtoCwZrzsjyhSLIeiCiWOURpOpp1IiWEr+LXxbSVOKSU9BWYlhB
0e4MUhh6GTEISQm2KXgzWkJCBam4g/AUz1e39V/fXd9Ye/jx+s5U6Z0pDTMTsR7MivB0BGLBX1IK
Luj7dljzOOQwzSuWpqikFRA6B9otGstYSdaID6ZWXII1CC1qIiB0DghNzioMPxR15AFCS/QYQWjN
p4XRBO+Ug84bijpgvHIOYbRofVIaCYRG/2lV0oS7rCrpNmdRk9rEOyuBPZNmBaI4NcFQJFpAaKFX
0qgQDJNKHpS5ZSB0lE7hpls5Z1ooGPrCKARTA0LDgML5wVYItU9D6DQPoeGj4DmA5GMIHTBqXNQh
wbD4U+GCn5HuJrjcyEKVJZHPajaK4lKXpiB0kaIkYs00/BBViaxYdWeZp8KD4RP4FmhrexGVVEYQ
2o2EQqoGoSVE8x/4c7aoowgXjIo6AOMqbpC4hf6RqAJV1QP4XmwHENpLb00efgpCK4BrEUowbFXS
A4S2AUKLqIQqKqNKkSYYNgidZ4o60jSEVkpLKRI1f0Jf1pC/VAO0JIX0QrwA/6ccEJoWS4PQ+HF6
bq1KWuYZhhsQulXChFSiJWHBz5miHjEJriLbAKFx944QOvBSGaHonLXmFUNVgCV6xwIS1tAAscwo
9Hdxe2d/a2Njb62/8r2D9d39FVUdxDjli2TDmmXoiktUdah3aWFVx2Keo4ycdAmeo4itm63qyAur
OkpUdeSrVHWUUChLq+pILdtKs1Ud5TJVHTaq6shReQfhMFTepbmqDpurvDOgjxiVUVWHXVTV0Zz0
UNWRRzxHbsTkEicdfKnIywVVHWMnjSOdr+pIraojL6rqoOYrLanqYFaHqg4qTVNSVUdZUdXhxrK4
qoNEPGhSOemZqg49ypKGhksST3npqg6Lqg48fv7pVR0mOgEGUJSSWGm5eU87sJi4mtk0z5GXVHXk
hZV381UdDZs0Jz2UhMg5igVYXNXRnHS+uKpDTDQh1sZVHblIArAxzxHuZkFVR3PSeVlVR6T9TCAL
KMuFZUE3IoZytJmqDiPvDPYfAFIuX9VBffTWrZ31NXtwv/fRu6MN8FRwUcSOJYmXfi1NnoLdwHgy
ef5ivKIDxEsMBYWkckgOWB8IU8xa4qFwfFtU/bL48MihLpPB6VmHoMdj11kMLAUwId2W0mFDKTwZ
B9x9zuzbkQeyw5slPUF9CgIxNZebVNSBag6cUjRtpq6nFYuiOAuFjA6mdLaoAyyrh1biCSVsQZy7
sgo8uK4M+M3BBIHYM7bQPLsou0QdUc2jvO3C7SUy1UJIoi1VU/jwCnMC2bI/viocxfLIyjBMRR0K
qTlHV00EG0gFARKqiOWeQoIqpqKOYE6tFXW4Ieqh8KI1kUTdkuGGCpOJEEZ2XkqkhbENBcMFix9y
WCJRUaSiafJleFb27WBJcQlfgqpBKEK0uEKcFebu8ze7b0duRR05SmFSmwXGLPSrFLsFTRkWKRQr
PLBLFHsxrhlyoaB/SfDLUi6ScAEitGbYWl1z0kS0fTsQULOGUwqv5ECcOhAoaV37AOnRm4KnQu4m
7TCtZ6XGZEtQlkqw9HLpLNeWBTAzT/1Q/qb4GYm/RTxUUQcVvxAxNl3UUbSIQ2LFVZRRUUfrkvBN
q/CRZKVEHc8iCKLSg2AYATFi0QjYKkRJrBZpagRU8L/4XHeTpbt9a/fWxtbG2r1H/5Z6L72+15D0
RqeqkxHfovTVxChnURZKz4NzcefvkEv8Y2BqU8ZboG2SyMFiei5CtETuQrLxlZ6Up7OeSziXAOme
+7b8AW6Wn4qwW6wlT17kopXui9SKVc3u8ylcdIDOYsExiLlI8WSXBgZBuYPDFaTK0hblokk3stwS
HPQwT+3tWUkJUiChHGlS0M0ksMoKg6okGVWWB/vmbsfvNtQdynPcspHn/v/s/VuTJVdypAu++6/g
Hzgpe8c9XiBitCat+gzR5LCtS6bnJQVEF9Etp8jqKYKcI/PrJ7d9qss9bolMFFCoIpMFApkRe/tl
LVtqaqq23BEJL3XLDF+oZJF0SGKBWqgst4RlFTh5iAx1RwibQrUXJb5EWuowQj7aBdFme6ZkyLfx
