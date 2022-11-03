#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;	# for script arguments for file input passing
use Pod::Usage;   	# for usage statement

my ($help, $oligo_target_file, $allele_corrected, $allele_uncorrected, $corrected_50nt, $uncorrected_50nt, $output_dir, $sample_name, $cutoff, $qc_pass_filter, $remove_ref_lt_alt, $only_need_primer_in_raw_data, $debug);

# command line arg handling
GetOptions(
        'help|h'               			=> \$help,
        'oligo_target_file=s'  			=> \$oligo_target_file,
        'allele_corrected=s'   			=> \$allele_corrected,
        'allele_uncorrected=s' 			=> \$allele_uncorrected,
        'corrected_50nt=s'     			=> \$corrected_50nt,
        'uncorrected_50nt=s'   			=> \$uncorrected_50nt,
        'output_dir=s'         			=> \$output_dir,
        'sample_name=s'       			=> \$sample_name,
        'cutoff=i'             			=> \$cutoff,
        'qc_pass_filter=f'     			=> \$qc_pass_filter,
        'remove_ref_lt_alt'    			=> \$remove_ref_lt_alt,
		'only_need_primer_in_raw_data'	=> \$only_need_primer_in_raw_data,
        'debug'                			=> \$debug
);

# debug flag will print out lots of info so you can verify the script is doing what we expect it to do
$debug = 0 unless ($debug);

# if the user running the script uses the help option or forgets to
# pass files, print the usage statement
pod2usage(-verbose =>1) if $help;
pod2usage("$0: Not enough files provided.\n") unless $oligo_target_file && $allele_corrected && 
													 $allele_uncorrected && $output_dir && $sample_name 
													 && $corrected_50nt && $uncorrected_50nt;

# unless you gave a cutoff, we'll use 0
$cutoff = 0 unless $cutoff;

# unless you give a qc_pass_filter value, we'll use .7
$qc_pass_filter = .7 unless $qc_pass_filter;

# unless you use remove_ref_lt_alt (remove data if reference counts are less than alt counts), it'll be set to 0 (or false)
$remove_ref_lt_alt = 0 unless $remove_ref_lt_alt;

# unless you set only_need_primer_in_raw_data (remove data )
$only_need_primer_in_raw_data = 0 unless $only_need_primer_in_raw_data;

my $pass_dir = $output_dir . "/pass";
my $fail_dir = $output_dir . "/fail";

# output directory name must be unique
if (-d $output_dir){
	die "Please select a new output directory name. $output_dir already exists.\n";
}else{
	mkdir $output_dir unless -d $output_dir;
	mkdir $pass_dir   unless -d $pass_dir;
	mkdir $fail_dir   unless -d $fail_dir;
}


# declaring hash and arrays for storing data from each of the files that will be read in
my (%oligo_target, @allele_corrected, @allele_uncorrected, @corrected_50nt, @uncorrected_50nt);

# Tell the user what files & parameters we are using
print "\nRunning verify_primary_alt_allele.pl with the following parameters:\n\n";
printf("%-30s %-150s\n",   "Oligo design file:",       $oligo_target_file);
printf("%-30s %-150s\n",   "Allele corrected file:",   $allele_corrected);
printf("%-30s %-150s\n",   "Allele uncorrected file:", $allele_uncorrected);
printf("%-30s %-150s\n",   "Corrected 50nt file:",     $corrected_50nt);
printf("%-30s %-150s\n",   "Uncorrected 50nt file:",   $uncorrected_50nt);
printf("%-30s %-150s\n",   "Output directory:",        $output_dir);
printf("%-30s %-150s\n",   "Cutoff:",                  $cutoff);
printf("%-30s %-150s\n",   "Minimum QC Pass Fraction:",$qc_pass_filter);
printf("%-30s %-150s\n",   "Remove Data with ref < alt:",   $remove_ref_lt_alt ? "yes" : "no");
printf("%-30s %-150s\n",   "Sample name:",             $sample_name);
printf("%-30s %-150s\n",   "Only require primer in raw data:", $only_need_primer_in_raw_data ? "yes" : "no");

# read in files and save to appropriate data structures
%oligo_target       = read_file($oligo_target_file,  1); # use hash as this data will be searched, everything else is in arrays.
@allele_corrected   = read_file($allele_corrected,   0);
@allele_uncorrected = read_file($allele_uncorrected, 0);
@corrected_50nt     = read_file($corrected_50nt,     0);
@uncorrected_50nt   = read_file($uncorrected_50nt,   0);

# need to save returned data
# need to declare objects 
my (@allele_corrected_pass, @allele_corrected_fail, @allele_uncorrected_pass, @allele_uncorrected_fail, @corrected_50nt_pass, @corrected_50nt_fail, @uncorrected_50nt_pass, @uncorrected_50nt_fail);
my ($allele_corrected_pass_ref, $allele_corrected_fail_ref, $allele_uncorrected_pass_ref, $allele_uncorrected_fail_ref, $corrected_50nt_pass_ref, $uncorrected_50nt_pass_ref) = remove_data_if_alt_allele_check_failed(\%oligo_target, \@allele_corrected, \@allele_uncorrected, \@corrected_50nt, \@uncorrected_50nt, $debug);
# need to de-reference references
@allele_corrected_pass = @$allele_corrected_pass_ref; @allele_corrected_fail = @$allele_corrected_fail_ref;
@allele_uncorrected_pass = @$allele_uncorrected_pass_ref; @allele_uncorrected_fail = @$allele_uncorrected_fail_ref;
@corrected_50nt_pass = @$corrected_50nt_pass_ref;
@uncorrected_50nt_pass = @$uncorrected_50nt_pass_ref;

# send allele corrected + uncorrected data from the alt_allele_check function to the filtering function
# also send the 50 nt files to the filtering function
print "\nRunning filter function to discard allele corrected, allele uncorrected, 50nt corrected, and 50nt uncorrected data below the minimum read cutoff of $cutoff and/or that conflicts with the design file.\n" unless $remove_ref_lt_alt;
print "\nRunning filter function to discard allele corrected, allele uncorrected, 50nt corrected, and 50nt uncorrected data below the minimum read cutoff of $cutoff, with reference counts less than alt counts, and/or that conflicts with the design file.\n" if $remove_ref_lt_alt;
my (%ac_above_cutoff, %ac_below_cutoff, %auc_above_cutoff, %auc_below_cutoff, @corrected_50nt_above_cutoff, @corrected_50nt_below_cutoff, @uncorrected_50nt_above_cutoff, @uncorrected_50nt_below_cutoff, %ac_missing_in_auc, %auc_missing_in_ac, 
%ac_above_cutoff_above_qc_filter, %ac_above_cutoff_below_qc_filter, %auc_above_cutoff_above_qc_filter, %auc_above_cutoff_below_qc_filter,
@corrected_50nt_above_cutoff_above_qc_filter, @corrected_50nt_above_cutoff_below_qc_filter, @uncorrected_50nt_above_cutoff_above_qc_filter, @uncorrected_50nt_above_cutoff_below_qc_filter, 
@corrected_50nt_above_cutoff_above_qc_ref_gt_alt, @corrected_50nt_above_cutoff_above_qc_ref_lt_alt, @uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt,
@uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt, @ac_above_cutoff_above_qc_ref_gt_alt,
@ac_above_cutoff_above_qc_ref_lt_alt, @auc_above_cutoff_above_qc_ref_gt_alt,  @auc_above_cutoff_above_qc_ref_lt_alt);

my ($ac_above_cutoff_ref, $ac_below_cutoff_ref, $auc_above_cutoff_ref, $auc_below_cutoff_ref, $corrected_50nt_above_cutoff_ref, $corrected_50nt_below_cutoff_ref, $uncorrected_50nt_above_cutoff_ref, $uncorrected_50nt_below_cutoff_ref, $ac_missing_in_auc_ref, $auc_missing_in_ac_ref,
$ac_above_cutoff_above_qc_filter_ref, $ac_above_cutoff_below_qc_filter_ref, $auc_above_cutoff_above_qc_filter_ref, $auc_above_cutoff_below_qc_filter_ref,
$corrected_50nt_above_cutoff_above_qc_filter_ref, $corrected_50nt_above_cutoff_below_qc_filter_ref, $uncorrected_50nt_above_cutoff_above_qc_filter_ref, $uncorrected_50nt_above_cutoff_below_qc_filter_ref, 
$corrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref, $corrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref, $uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref,
$uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref, $ac_above_cutoff_above_qc_ref_gt_alt_ref, $ac_above_cutoff_above_qc_ref_lt_alt_ref,
$auc_above_cutoff_above_qc_ref_gt_alt_ref, $auc_above_cutoff_above_qc_ref_lt_alt_ref) = filter(\@allele_corrected_pass, \@allele_uncorrected_pass, \@corrected_50nt_pass, \@uncorrected_50nt_pass, $cutoff, $qc_pass_filter, $remove_ref_lt_alt, $only_need_primer_in_raw_data, $debug);

%ac_above_cutoff = %$ac_above_cutoff_ref; %ac_below_cutoff = %$ac_below_cutoff_ref;
%auc_above_cutoff = %$auc_above_cutoff_ref; %auc_below_cutoff = %$auc_below_cutoff_ref;
@corrected_50nt_above_cutoff = @$corrected_50nt_above_cutoff_ref; @corrected_50nt_below_cutoff = @$corrected_50nt_below_cutoff_ref;
@uncorrected_50nt_above_cutoff = @$uncorrected_50nt_above_cutoff_ref; @uncorrected_50nt_below_cutoff = @$uncorrected_50nt_below_cutoff_ref;
%ac_missing_in_auc = %$ac_missing_in_auc_ref; %auc_missing_in_ac = %$auc_missing_in_ac_ref;

%ac_above_cutoff_above_qc_filter = %$ac_above_cutoff_above_qc_filter_ref;
%ac_above_cutoff_below_qc_filter = %$ac_above_cutoff_below_qc_filter_ref;
%auc_above_cutoff_above_qc_filter = %$auc_above_cutoff_above_qc_filter_ref;
%auc_above_cutoff_below_qc_filter = %$auc_above_cutoff_below_qc_filter_ref;
@corrected_50nt_above_cutoff_above_qc_filter = @$corrected_50nt_above_cutoff_above_qc_filter_ref;
@corrected_50nt_above_cutoff_below_qc_filter = @$corrected_50nt_above_cutoff_below_qc_filter_ref;
@uncorrected_50nt_above_cutoff_above_qc_filter = @$uncorrected_50nt_above_cutoff_above_qc_filter_ref;
@uncorrected_50nt_above_cutoff_below_qc_filter = @$uncorrected_50nt_above_cutoff_below_qc_filter_ref;

@corrected_50nt_above_cutoff_above_qc_ref_gt_alt = @$corrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref;
@corrected_50nt_above_cutoff_above_qc_ref_lt_alt = @$corrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref;
@uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt = @$uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref;
@uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt = @$uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref;
@ac_above_cutoff_above_qc_ref_gt_alt = @$ac_above_cutoff_above_qc_ref_gt_alt_ref;
@ac_above_cutoff_above_qc_ref_lt_alt = @$ac_above_cutoff_above_qc_ref_lt_alt_ref;
@auc_above_cutoff_above_qc_ref_gt_alt = @$auc_above_cutoff_above_qc_ref_gt_alt_ref; 
@auc_above_cutoff_above_qc_ref_lt_alt = @$auc_above_cutoff_above_qc_ref_lt_alt_ref;


print "List of output files:\n";
# write the files
# need output filename (using basename), array, and output directory

# files from alt_allele_check
print "\tFiles from alt allele check function. The \"fail\" files may be empty- this means your input files were properly formatted!\n";
write_file(\@allele_corrected_pass,   $pass_dir, $sample_name, "alleleErrC.passAltAlleleCheck.txt", 0);
write_file(\@allele_corrected_fail,   $fail_dir, $sample_name, "alleleErrC.failAltAlleleCheck.txt", 0);
write_file(\@allele_uncorrected_pass, $pass_dir, $sample_name, "alleleNoC.passAltAlleleCheck.txt",  0);
write_file(\@allele_uncorrected_fail, $fail_dir, $sample_name, "alleleNoC.failAltAlleleCheck.txt",  0);

# files from filter function
print "\n\tFiles from filter function. Data above the specified cutoff $cutoff will be in \"pass\" directory.\n";
print "\n\tData above the DP cutoff $cutoff\n";
write_file(\%ac_above_cutoff,               $pass_dir, $sample_name, "alleleErrC.passAltAlleleCheck.aboveCutoff.txt", 1); 
write_file(\%auc_above_cutoff,              $pass_dir, $sample_name, "alleleNoC.passAltAlleleCheck.aboveCutoff.txt", 1); 
write_file(\@corrected_50nt_above_cutoff,   $pass_dir, $sample_name, "50ntErrC.aboveCutoff.txt", 0); 
write_file(\@uncorrected_50nt_above_cutoff, $pass_dir, $sample_name, "50ntNoC.aboveCutoff.txt", 0); 
print "\n\tData below the DP cutoff $cutoff\n";
write_file(\%ac_below_cutoff,               $fail_dir, $sample_name, "alleleErrC.passAltAlleleCheck.belowCutoff.txt", 1);
write_file(\%auc_below_cutoff,              $fail_dir, $sample_name, "alleleNoC.passAltAlleleCheck.belowCutoff.txt", 1);
write_file(\@corrected_50nt_below_cutoff,   $fail_dir, $sample_name, "50ntErrC.belowCutoff.txt", 0);
write_file(\@uncorrected_50nt_below_cutoff, $fail_dir, $sample_name, "50ntNoC.belowCutoff.txt", 0 );
print "\n\tData missing in one of the allele files\n";
write_file(\%ac_missing_in_auc,             $fail_dir, $sample_name, "alleleErrC.passAltAlleleCheck.missingPrimersInAlleleNoC.txt", 1);
write_file(\%auc_missing_in_ac,             $fail_dir, $sample_name, "alleleNoC.passAltAlleleCheck.missingPrimersInAlleleErrC.txt", 1);

# even more files from filter function
print "\n\tData above the QC pass filter $qc_pass_filter\n";
write_file(\%ac_above_cutoff_above_qc_filter,  $pass_dir, $sample_name, "alleleErrC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.txt", 1);
write_file(\%auc_above_cutoff_above_qc_filter, $pass_dir, $sample_name, "alleleNoC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.txt", 1);
write_file(\@corrected_50nt_above_cutoff_above_qc_filter, $pass_dir, $sample_name, "50ntErrC.aboveCutoff.aboveQCfilter.txt", 0);
write_file(\@uncorrected_50nt_above_cutoff_above_qc_filter, $pass_dir, $sample_name, "50ntNoC.aboveCutoff.aboveQCfilter.txt", 0);

print "\n\tData below the QC pass filter $qc_pass_filter\n";
write_file(\%ac_above_cutoff_below_qc_filter,  $fail_dir, $sample_name, "alleleErrC.passAltAlleleCheck.aboveCutoff.belowQCfilter.txt", 1);
write_file(\%auc_above_cutoff_below_qc_filter, $fail_dir, $sample_name, "alleleNoC.passAltAlleleCheck.aboveCutoff.belowQCfilter.txt", 1);
write_file(\@corrected_50nt_above_cutoff_below_qc_filter, $fail_dir, $sample_name, "50ntErrC.aboveCutoff.belowQCfilter.txt", 0);
write_file(\@uncorrected_50nt_above_cutoff_below_qc_filter, $fail_dir, $sample_name, "50ntNoC.aboveCutoff.belowQCfilter.txt", 0);


if($remove_ref_lt_alt){
	print "\n\tData with ref > alt\n";
	write_file(\@ac_above_cutoff_above_qc_ref_gt_alt,  $pass_dir, $sample_name, "alleleErrC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refGTalt.txt", 0);
	write_file(\@auc_above_cutoff_above_qc_ref_gt_alt, $pass_dir, $sample_name, "alleleNoC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refGTalt.txt", 0);
	write_file(\@corrected_50nt_above_cutoff_above_qc_ref_gt_alt, $pass_dir, $sample_name, "50ntErrC.aboveCutoff.aboveQCfilter.refGTalt.txt", 0);
	write_file(\@uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt, $pass_dir, $sample_name, "50ntNoC.aboveCutoff.aboveQCfilter.refGTalt.txt", 0);

	print "\n\tData with ref < alt\n";
	write_file(\@ac_above_cutoff_above_qc_ref_lt_alt ,  $fail_dir, $sample_name, "alleleErrC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refLTalt.txt", 0);
	write_file(\@auc_above_cutoff_above_qc_ref_lt_alt, $fail_dir, $sample_name, "alleleNoC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refLTalt.txt", 0);
	write_file(\@corrected_50nt_above_cutoff_above_qc_ref_lt_alt, $fail_dir, $sample_name, "50ntErrC.aboveCutoff.aboveQCfilter.refLTalt.txt", 0);
	write_file(\@uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt , $fail_dir, $sample_name, "50ntNoC.aboveCutoff.aboveQCfilter.refLTalt.txt", 0);
}


# read file subroutine to ingest input files
sub read_file{
	my ($file, $is_this_oligo_target_file) = @_;
	my (%design_file, @file_contents);

	open (my $read_fh, "<", $file) || die "Cannot open $file for reading.\n";

	while (my $line =<$read_fh>){
		chomp $line;

		# for DESIGN FILE, use a hash as this data will be searched
		# everything else is in arrays.

		# DESIGN FILE, columns of interest (indexing by 0 here):
		# 2: Chromosome
		# 3: Loci
		# 5: Ref
		# 6: Alt
		if ($is_this_oligo_target_file){
			my @temp = split(/\t/, $line);

			my $primer = $temp[0]; my $chr = $temp[2]; my $loci = $temp[3]; my $ref = $temp[5]; my $alt = $temp[6];
			my $key = $temp[2] ."_" . $loci ."_" . $primer;

			# if we already have this loci (chromosome + loci combo stored in $key) in the %design_file hash,
			# check for errors. Errors include same chr+loci, but different ref alleles, OR same chr+loci, but different primary alt alleles.
			# if there are not any errors, we continue to adding it to the hash - which we really don't need to do (it'll overwrite what's there already)
			if ($design_file{$key}){
				die "Reference alleles for the same loci do not match: chr:$chr loci:$loci ref1:$design_file{$key}->{'ref'} ref2:$ref\n" if $design_file{$key}->{"ref"} ne $ref;
				die "Primary alternate alleles for the same loci do not match: chr:$chr loci:$loci alt1:$design_file{$key}->{'alt'} alt2:$alt\n" if $design_file{$key}->{'alt'} ne $alt;
				$design_file{$key}{"ref"} = $ref;
				$design_file{$key}{"alt"} = $alt;
			# add ref and alt info for chromosome + loci that is currently not in the hash
			}else{
				$design_file{$key}{"ref"} = $ref;
				$design_file{$key}{"alt"} = $alt;
			}
		# if we are not processing the design file,
		# then just add each line of the input file to an array
		# and then return the file contents array
		}else{
			chomp $line;
			push(@file_contents, $line);
		}

	}

	return %design_file if %design_file;
	return @file_contents if @file_contents;
	close $read_fh;
}

# write file subroutine, to create output files in the user-specified output directory
# input:
# data to write, directory to write to, prefix of file (the sample name), the suffix (what describes the content of the file), and if this is a hash or not
sub write_file{
	my ($data_to_write_ref, $dir, $prefix, $suffix, $is_hash) = @_;

	my $full_path_to_output_file = $dir . "/" . $prefix . "-" . $suffix;
	print "\t$full_path_to_output_file\n";
	open(my $write_fh, ">", $full_path_to_output_file) || die "Unable to open filehandle for writing $dir/$prefix$suffix:$!\n";

	if ($is_hash){
		my %data_to_write = %$data_to_write_ref;
	
		for my $key (sort keys %data_to_write){
			my $line = $data_to_write{$key};
			chomp $line;
			print $write_fh $line, "\n";
		}
		close $write_fh;	
	}else{
		my @data_to_write = @$data_to_write_ref;

		for my $line (@data_to_write){
			chomp $line;
			print $write_fh $line, "\n";
		}

		close $write_fh;
	}	
}

sub remove_data_if_alt_allele_check_failed{
	my ($oligo_target_file_ref, $allele_corrected_ref, $allele_uncorrected_ref, $corrected_50nt_ref, $uncorrected_50nt_ref , $print_if_debug) = @_;

	# de-reference the references that were passed to the subroutine
	my %oligo_target = %{$oligo_target_file_ref};
	my @allele_corrected = @{$allele_corrected_ref};
	my @allele_uncorrected = @{$allele_uncorrected_ref};
	my @corrected_50nt = @{$corrected_50nt_ref}; 
	my @uncorrected_50nt = @{$uncorrected_50nt_ref};

	my (@allele_corrected_pass, @allele_corrected_fail, @allele_uncorrected_pass, @allele_uncorrected_fail);

	# send in all allele corr, allele error, 50nt
	# also design file
	print "\nRunning alt allele check to compare with the design file:\n\tAllele Corrected\n";
	my ($allele_corrected_pass_ref, $allele_corrected_fail_ref)  = alt_allele_check(\%oligo_target, \@allele_corrected, $debug);
	@allele_corrected_pass = @$allele_corrected_pass_ref; @allele_corrected_fail = @$allele_corrected_fail_ref;

	# send allele uncorrected file to alternate allele check function
	print "\tAllele Uncorrected\n";
	my ($allele_uncorrected_pass_ref, $allele_uncorrected_fail_ref) = alt_allele_check(\%oligo_target, \@allele_uncorrected, $debug);
	@allele_uncorrected_pass = @$allele_uncorrected_pass_ref; @allele_uncorrected_fail = @$allele_uncorrected_fail_ref;

	print "\n\t* \"Passing\" data includes:\n\t\t- matching chr + loci + primer in the the design file with either a reference match,
	\n\t\t- the expected primary alt allele,\n\t\t- the expected alt allele in a non-primary position,
	\n\t\t- and alternative alleles that are not the expected one.\n\tFor the last two cases, modification of data will be done.\n";
	print "\n\t* \"Failing\" data includes:\n\t\t- no matching chr + loci + primer in the design file,
	\n\t\t- or wrong reference for a chr + loci + primer found in the design file.\n\n";

	print "Running function to collate the list of primers that conflict that the design file, for both allele corrected and allele uncorrected.\n";
	my ($collated_failed_primers_ref) = find_primers_that_conflict_with_design_file(\@allele_corrected_fail, \@allele_uncorrected_fail, $debug);
	my @fail_primers; @fail_primers = @$collated_failed_primers_ref;
	my $fail_primers_size = $#fail_primers + 1;
	print "\tNumber of primers to remove: $fail_primers_size\n" if @fail_primers;
	print "\t\t",join(", ", sort @fail_primers),"\n\n";
	print "\tNo primers found that conflict with the design file.\n" unless @fail_primers;

	# need to remove data from 50nt if @fail_primers exists
	# TODO: need to test this works as expected
	if (@fail_primers){
		# remove any 50nt lines that have failed primers
		print "\tRemoving primers (that conflict with the design file) from the 50nt files.\n\n";
		my ($corrected_50nt_pass_ref) = remove_lines_by_specified_data(@corrected_50nt, @fail_primers, 37 , $debug); 
		my ($uncorrected_50nt_pass_ref) = remove_lines_by_specified_data(@uncorrected_50nt, @fail_primers, 37, $debug); 

		my @corrected_50nt_pass = @$corrected_50nt_pass_ref; my @uncorrected_50nt_pass = @$uncorrected_50nt_pass_ref;
		return (\@allele_corrected_pass, \@allele_corrected_fail, \@allele_uncorrected_pass, \@allele_uncorrected_fail, \@corrected_50nt_pass, \@uncorrected_50nt_pass);
	}else{
		print "\tNo removal of primers required for 50nt files, as no primers conflict with the design file.\n\n";
		return (\@allele_corrected_pass, \@allele_corrected_fail, \@allele_uncorrected_pass, \@allele_uncorrected_fail, \@corrected_50nt, \@uncorrected_50nt);
	}
}

# Function for comparison of the design file & a given file, which is either allele_corrected, or allele_uncorrected
# We need to see if the primary alternate allele that's in a given file matches
# what we'd expect (e.g. what the design file contains).

# check that the ref allele was called -> PASS
# check that the expected alt allele was called -> PASS
# check that the expected allele was called not at the primary position (in secondary or tertiary, etc.) --> modification of data, PASS
# check that a non-expected alt allele was called --> modification of data, PASS

# things that need to be modified for the last two conditions listed above:
    # modify ALT, total DP, AAF and %nt QC

sub alt_allele_check{
	my ($design_file_ref, $file_to_compare_ref, $print_if_debug) = @_;

	# de-reference the references that were passed to the subroutine
	my %design_file = %$design_file_ref;
	my @file_to_compare = @$file_to_compare_ref;

	# examine the file of interest (allele corrected or uncorrected)
	# see what matched up with the design file and what does not.
	my (@pass, @fail);
	foreach my $line (@file_to_compare){
		# split line on tabs, then save the relevant fields to variables
		my @temp_line = split(/\t/, $line);

		# ALLELE CORRECTED\UNCORRECTED columns of interest
		# (indexed starting at 0)
		# 0: chromosome
		# 1: loci or site 1
		# 2: loci or site 2
		# 3: reference allele
		# 4: alt allele(s)
		# 5: total read count
		# 6: AD (ref, alt 1, alt2, ...)
		# 7: ALTREADS
		# 8: DPREADS
		# 9: AAF
		# 10: %nt QC

		my $chr = $temp_line[0]; my $loci = $temp_line[1]; my $ref = $temp_line[3]; my $alt = $temp_line[4];
		my $primer = $temp_line[37];

		# check that we do have this chromosome + loci combo in the design file
		if ($design_file{$chr."_".$loci."_".$primer}){
			# going to extract the first character from the alt string to get the primary alt allele.
			# alt column could be only one character, such as A
			# it also could be multiple comma separated characters, such as T,C,G
			# and if there is no alt allele - the column denotes the reference, represented as <*>.

			my $primary_alt = substr($alt, 0, 1);
			#print "Design file ref and alt: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $design_file{$chr."_".$loci}->{'alt'}, "\n" if $print_if_debug;
			#print "Input file ref and alt: ", $ref, "\t", $alt, "\n" if $print_if_debug;
			#print $line,"\n" if $print_if_debug;

			# verify reference in this line is the same as in the design file:
			if($design_file{$chr."_".$loci."_".$primer}->{'ref'} eq $ref){

				# ALT ALLELE CHECK CONDITIONS:
				# 1. ref allele only
				if ( $alt eq "<*>" ){
					push(@pass, $line);
					#print "PASS due to reference match: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\t", "and alt is $alt\n\n" if $print_if_debug;

				# 2. expected alt allele in primary allele position
				}elsif ( $design_file{$chr."_".$loci."_".$primer}->{'alt'} eq $primary_alt ){
					push(@pass, $line);
					#print "PASS due to primary alt match: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\t", "and\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$primary_alt\t $alt\n\n" if $print_if_debug;

				# 3. expected alt allele in secondary/tertiary/etc position
				}elsif ( $alt =~m/($design_file{$chr."_".$loci."_".$primer}->{'alt'})/ ){
					# Need to modify:
					# ALT (col 8)
					# TOTAL DP (col 9)  = ref + ALT
					# AAF (col 10)  = new ALT / new total DP = ((col 8)/ (col9))
					# % nt QC (col 11) = new col 9 / col 6
					# position of match in $alt
					my $alt_tmp = $alt;
					$alt_tmp =~ tr/,//d; # remove commas to get position
							     # used to pull out correct alt counts
					my $position = index($alt_tmp, $design_file{$chr."_".$loci."_".$primer}->{'alt'}) ;

					my $AD = $temp_line[6];
					$AD =~ tr/AD=//d;
					my @temp_for_alt_counts =split(",", $AD);

					# pull out correct alt counts based upon the position of the expected alt allele in col 8. So in the AD column, the correct alt counts will be position plus one (as reference is listed first in the AD string)
					my $correct_alt_counts = $temp_for_alt_counts[$position + 1];
					my $recalc_dp = $temp_for_alt_counts[0] + $correct_alt_counts; # ref + alt
					my $aaf = sprintf( "%.5e", ($correct_alt_counts/ $recalc_dp));
					my $percent_nt_qc = sprintf( "%.5e", ($recalc_dp / $temp_line[5]));

					my @end_of_line = @temp_line[11...38];
					my $end_of_line = join("\t", @end_of_line) ;
					my $modified_line = "$chr\t$loci\t$temp_line[2]\t$ref\t$alt\t$temp_line[5]\tAD=$AD\t$correct_alt_counts\t$recalc_dp\t$aaf\t$percent_nt_qc\t$end_of_line\n";

					#print $line, "\n", $modified_line if $print_if_debug;
					push(@pass, $modified_line);
					#print "PASS due to: ",  $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\tand\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$alt\n\n" if $print_if_debug;

				# check for alternative alleles that are not the expected one (e.g. not in the design file)
				# modify alt, total dp, aaf if so

				# 4. alternative alleles that are not the expected one (e.g. not in the design file)
				}elsif( $alt !~/$design_file{$chr."_".$loci."_".$primer}->{'alt'}/ ){
					#print "PASS (needs modifications) due to: ",  $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\tand\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$alt\n" if $print_if_debug;
					# Alt goes to 0
					# Total dp is (original total dp – original alt) or AD=(probably easier)
					# AAF goes to 0
					# % reads that are good quality (%nt QC), put this to 0??? Not sure yet
					my $AD = $temp_line[6];
					my $ref_count_from_AD = $1 if $AD =~m/AD=(\d+),/;
					my $percent_nt_qc = sprintf( "%.5e", ($ref_count_from_AD / $temp_line[5]));
					my @end_of_line = @temp_line[11...38];
					my $end_of_line = join("\t", @end_of_line);
					my $modified_line = "$chr\t$loci\t$temp_line[2]\t$ref\t$alt\t$temp_line[5]\t$AD\t0\t$ref_count_from_AD\t0\t$percent_nt_qc\t$end_of_line\n";
					#print $line, "\n\n" if $print_if_debug;
					#print $modified_line,"\n\n" if $print_if_debug;

					# save modified line
					push(@pass, $modified_line);

				# catch for any other condition we were not expecting
				}else{
					print "\tException in logic found in alt_allele_check function, due to the following line!\n$line\n";
					push(@fail, $line);
					print "\tFAIL due to a condition that is not covered by the logic in the alt_allele_check function.\n" if $print_if_debug;
				}


			# reference in this line IS NOT the same as in the design file, yield an error:
			}else{
				print "\tReference listed for chromosome $chr loci $loci primer $primer does not match the ref in the design file.\n";
				push (@fail, $line);
				print "\tFAIL due to different ref in design file for chromosome $chr loci $loci primer $primer. DF has  ", $design_file{$chr."_".$loci."_".$primer}->{'ref'}, "  and this line has $ref\n" if $print_if_debug;
			}


		}else{
			# warn if we do not find a matching chr + loci combination in the design file.
			print "\tChromosome, loci, and primer combination $chr $loci $primer not found in design file.\n";
			push(@fail, $line);
			print "\tFAIL due to: chromosome $chr & loci $loci & primer $primer missing in design file.\n" if $print_if_debug;
		}
	} #foreach close
	my $input_array_size = $#file_to_compare + 1; my $pass_array_size = $#pass +1; my $fail_array_size = $#fail +1;
	print "\t\tSize of input data to function: ($input_array_size)\n\t\tSize of output \"passing\" data: ($pass_array_size)\n\t\tSize of output \"failing\" data: ($fail_array_size)\n\n";
	return (\@pass, \@fail);
}

sub find_primers_that_conflict_with_design_file {
	
	my ($allele_corrected_fail_ref, $allele_uncorrected_fail_ref, $print_if_debug) = @_;

	# de-reference the references that were passed to the subroutine
	my @allele_corrected_fail = @$allele_corrected_fail_ref;
	my @allele_uncorrected_fail = @$allele_uncorrected_fail_ref;

	my @megaarray = (@allele_corrected_fail, @allele_uncorrected_fail);

	my @primers_to_remove;
	for my $line (@megaarray){
		#print $line if $print_if_debug;
		my @temp_array = split("\t", $line);
		#print $temp_array[37], "\n" if $print_if_debug;
		push @primers_to_remove, $temp_array[37];
	}
	return \@primers_to_remove;
}

# TODO: need to test this works as expected
# input: (1) array to remove stuff from, (2) array of stuff to remove, (3) which column to search for stuff to remove, (4) debug
sub remove_lines_by_specified_data{
	my ($data_to_clean_up_ref, $data_that_needs_to_be_removed_ref, $column_of_interest_to_search, $print_if_debug) = @_;

	# de-reference the references that were passed to the subroutine
	my @data_to_clean_up = @$data_to_clean_up_ref;
	my @data_that_needs_to_be_removed = @$data_that_needs_to_be_removed_ref;
	
	my %use_for_searching;
	for my $datum (@data_that_needs_to_be_removed){
		$use_for_searching{$datum} = 1;
	}

	my @to_keep;
	for my $line (@data_to_clean_up){
		my @temp_line = split("\t", $line);
		if ($use_for_searching{$temp_line[$column_of_interest_to_search]}){
			print "found line to remove:$line\n" if $debug;
		}else{
			push(@to_keep , $line);
		}
	}
	return \@to_keep;
}

# this function for filtering for passing read depth.
# if corrected/uncorrected read depth total does not match (aka both corrected/uncorrected are >=cutoff) -> throw out amplicon
# if amplicon is thrown out, discard 50nt entries (NOT VARIANT)

sub filter{
	my ($allele_corrected_pass_ref, $allele_uncorrected_pass_ref, $corrected_50nt_ref, $uncorrected_50nt_ref, $cutoff, $qc_pass_filter, $remove_ref_lt_alt, $only_need_primer_in_raw_data, $print_if_debug) = @_;

	# de-reference the references that were passed to the subroutine
	my @allele_corrected_pass = @$allele_corrected_pass_ref;
	my @allele_uncorrected_pass = @$allele_uncorrected_pass_ref;
	my @corrected_50nt = @$corrected_50nt_ref;
	my @uncorrected_50nt = @$uncorrected_50nt_ref;
	
	# process allele uncorrected "pass" data and segment data into >= DP cutoff, and < DP cutoff
	my %auc_above_cutoff; my %auc_below_cutoff;

	# $temp_line[8] - total DP
	# $temp_line[37] - primer ID (three digit code)
	foreach my $line (@allele_uncorrected_pass){
		my @temp_line = split("\t", $line);
		if ($temp_line[8] >= $cutoff){
			$auc_above_cutoff{$temp_line[37]}=$line;
		}elsif($temp_line[8] < $cutoff){
			$auc_below_cutoff{$temp_line[37]}=$line;
		}else{
			print "Error with allele corrected pass data!:\n$line\n";
		}
	}

	# go through allele corrected pass
	# 	with default option
	# 		for the discarded ones in allele corrected pass -> discard from allele uncorrected too
	# 		if any discarded from allele uncorrected (<cutoff) -> discard from allele corrected, too
	# 		then for the discarded ones, remove 50nt as well
	#
	# 	with only_need_primer_in_raw_data option
	#		for the discarded ones in allele uncorrected pass -> discard from allele corrected pass, too
	#														  -> discard from both 50nts
	#		for the discarded ones in allele corrected pass   -> keep in allele uncorrected pass
	#														  -> discard from 50nt for allele corrected

	my %ac_above_cutoff; my %ac_below_cutoff; my %ac_missing_in_auc;
	my %remove_from_50nt_files; # this means the amplicon was < cutoff in at least one of allele corrected/uncorrected
	my %remove_from_50nt_corrected;

	foreach my $line (@allele_corrected_pass){
		my @temp_line = split("\t", $line);
		chomp($line);
		#print "$line\t$temp_line[8]\n" if $print_if_debug;

		# if above cutoff, must check to see if allele_uncorrected is also above the cutoff for this amplicon
		if ( $temp_line[8]>=$cutoff  ){
  
			# allele corrected above cutoff, allele uncorrected below cutoff --> DISCARD
			if ($auc_below_cutoff{$temp_line[37]}){
				# add to remove_from_50nt_files
				$remove_from_50nt_files{$temp_line[37]} = 1;
				$ac_below_cutoff{$temp_line[37]}=$line;
			# both above cutoff ---> KEEP
			}elsif($auc_above_cutoff{$temp_line[37]}){
				$ac_above_cutoff{$temp_line[37]}=$line;
			# allele corrected primer ID is not in allele uncorrected---> PUT IN SEPARATE BUCKET
			}else{
				print "\tError in filter function. Missing primer $temp_line[37] for allele uncorrected\n" unless $auc_above_cutoff{$temp_line[37]} || $auc_below_cutoff{$temp_line[37]};
				print "Unspecified error in filter function" if $auc_above_cutoff{$temp_line[37]} || $auc_below_cutoff{$temp_line[37]};
				$ac_missing_in_auc{$temp_line[37]} = $line;
				$remove_from_50nt_files{$temp_line[37]} = 1;
			}

		# allele corrected below cutoff:
		# check to see if allele_corrected is above the cutoff.
		}else{

			# allele corrected below cutoff, allele uncorrected below cutoff --> DISCARD
			if ($auc_below_cutoff{$temp_line[37]}){
				$remove_from_50nt_files{$temp_line[37]} = 1;
				$ac_below_cutoff{$temp_line[37]}=$line;
			# allele corrected below, allele uncorrected above cutoff ---> DISCARD
			# TODO: keep anyway
			}elsif($auc_above_cutoff{$temp_line[37]}){
				if ($only_need_primer_in_raw_data == 1){
					$remove_from_50nt_corrected{$temp_line[37]} = 1;
					$ac_below_cutoff{$temp_line[37]}=$line;
				}else{
					print $temp_line[37], "\n";
					$remove_from_50nt_files{$temp_line[37]} = 1;
					# remove from above cutoff array for allele corrected, as the allele uncorrected is below
					my $to_keep = $auc_above_cutoff{$temp_line[37]};
					delete $auc_above_cutoff{$temp_line[37]};
					$auc_below_cutoff{$temp_line[37]} = $to_keep;
					$ac_below_cutoff{$temp_line[37]}=$line;
				}
			# allele corrected primer ID is not in allele uncorrected---> PUT IN SEPARATE BUCKET
			}else{
				print "\tError in filter function. Missing primer $temp_line[37] for allele uncorrected\n" unless $auc_below_cutoff{$temp_line[37]} || $auc_above_cutoff{$temp_line[37]};
				print "Unspecified error in filter function" if $auc_below_cutoff{$temp_line[37]} || $auc_above_cutoff{$temp_line[37]};
				$ac_missing_in_auc{$temp_line[37]} = $line;
				$remove_from_50nt_files{$temp_line[37]} = 1;
			}

		}

	}

	# check for primers that are in allele uncorrected data but not in allele corrected
	# shunt these to a bucket and print 'em out later
	# TODO: this won't be run when the toggle is on.
	my %auc_missing_in_ac;
	for my $line (@allele_uncorrected_pass){
		my @temp_line = split("\t", $line);
		my $primer = $temp_line[37];
		if (! $ac_above_cutoff{$primer}){
			print "\tError in filter function. Missing primer $primer for allele corrected\n";
			$auc_missing_in_ac{$primer} = $line;
			if ($only_need_primer_in_raw_data == 0){
				delete $auc_above_cutoff{$primer};
				#$auc_below_cutoff{$primer} = $auc_missing_in_ac{$primer};
				$remove_from_50nt_files{$primer} = 1;
			}
		}
	}	 

	# remove the irrelevant lines (those corresponding to amplicons < the DP cutoff) from 50nt files
	# TODO: need to process 50nt separately for each corr/uncorr if toggle is on
	my @corrected_50nt_above_cutoff; my @corrected_50nt_below_cutoff;
	my @uncorrected_50nt_above_cutoff; my @uncorrected_50nt_below_cutoff;
	for my $c_50nt (@corrected_50nt){
		my @temp_line = split("\t", $c_50nt);
		if (($only_need_primer_in_raw_data == 1) && ($remove_from_50nt_corrected{$temp_line[37]})){
			push (@corrected_50nt_below_cutoff, $c_50nt);
		}elsif ($remove_from_50nt_files{$temp_line[37]}){
			push (@corrected_50nt_below_cutoff, $c_50nt);
		}else{
			push (@corrected_50nt_above_cutoff, $c_50nt);
		}		
	}

	for my $uc_50nt (@uncorrected_50nt){
		my @temp_line = split("\t", $uc_50nt);
		if ($remove_from_50nt_files{$temp_line[37]}){
			push(@uncorrected_50nt_below_cutoff, $uc_50nt);
		}else{
			push(@uncorrected_50nt_above_cutoff , $uc_50nt);
		}
	}

	my (%ac_above_cutoff_above_qc_filter, %ac_above_cutoff_below_qc_filter, %auc_above_cutoff_above_qc_filter, %auc_above_cutoff_below_qc_filter);
	my (@corrected_50nt_above_cutoff_above_qc_filter, @corrected_50nt_above_cutoff_below_qc_filter, @uncorrected_50nt_above_cutoff_above_qc_filter, @uncorrected_50nt_above_cutoff_below_qc_filter);

	if ($qc_pass_filter){
		# alelle corrected
		my %remove_from_auc_files;
		foreach my $key (keys %ac_above_cutoff){
			my $line = $ac_above_cutoff{$key};
			my @temp = split("\t", $line);

			if ($temp[10] < $qc_pass_filter){
				$ac_above_cutoff_below_qc_filter{$temp[37]} = $line;
				$remove_from_auc_files{$temp[37]} = 1;
			}else{
				$ac_above_cutoff_above_qc_filter{$temp[37]} = $line;
			}
		}

		# allele uncorrected
		foreach my $key (keys %auc_above_cutoff){
			my $line = $auc_above_cutoff{$key};
			my @temp = split("\t", $line);

			if (($only_need_primer_in_raw_data == 0) && (exists $remove_from_auc_files{$temp[37]})){
				$auc_above_cutoff_below_qc_filter{$temp[37]} = $line;
			}elsif($temp[10] < $qc_pass_filter){
				$auc_above_cutoff_below_qc_filter{$temp[37]} = $line;
				delete $ac_above_cutoff_above_qc_filter{$temp[37]}; 
				$ac_above_cutoff_below_qc_filter{$temp[37]} = $line;
				# remove from ac above cutoff if it exists
			}else{
				$auc_above_cutoff_above_qc_filter{$temp[37]} = $line;
			}
		}

		# keep allele corrected and 50nt corrected consistent
		for my $line (@corrected_50nt_above_cutoff){
			my @temp = split("\t", $line);
			if (exists $ac_above_cutoff_below_qc_filter{$temp[37]}){
				push (@corrected_50nt_above_cutoff_below_qc_filter, $line);
			}else{
				push (@corrected_50nt_above_cutoff_above_qc_filter, $line);
			}

		}

		# keep allele uncorrected and 50nt uncorrected consistent
		for my $line (@uncorrected_50nt_above_cutoff){
			my @temp = split("\t", $line);
		
			if (exists $auc_above_cutoff_below_qc_filter{$temp[37]}){
				push (@uncorrected_50nt_above_cutoff_below_qc_filter, $line);
			}else{
				push(@uncorrected_50nt_above_cutoff_above_qc_filter, $line);
			}

		}
 
	}else{
		print "\tError in filter function. Minimum QC Pass Fraction is unset."; # should never get this error
	}


	#TODO:
	# if ref < alt toggle is on, then make sure to remove the lines where ref allele counts are less than alt allele counts.
	# must do this for 4 files, then send the ref < alt lines to new files.

	# 50ntErrC, 50ntNoC, alleleErrC, alleleNoC
	my ($corrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref, 
	$corrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref, $uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref, 
	$uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref, $ac_above_cutoff_above_qc_ref_gt_alt_ref, $ac_above_cutoff_above_qc_ref_lt_alt_ref, $auc_above_cutoff_above_qc_ref_gt_alt_ref, $auc_above_cutoff_above_qc_ref_lt_alt_ref) = remove_data_if_ref_lt_alt(\@corrected_50nt_above_cutoff_above_qc_filter, \@uncorrected_50nt_above_cutoff_above_qc_filter, \%ac_above_cutoff_above_qc_filter, \%auc_above_cutoff_above_qc_filter, $only_need_primer_in_raw_data, $print_if_debug) if $remove_ref_lt_alt;

	my @corrected_50nt_above_cutoff_above_qc_ref_gt_alt = @$corrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref; 
	my @corrected_50nt_above_cutoff_above_qc_ref_lt_alt = @$corrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref;
	my @uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt = @$uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt_ref; 
	my @uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt = @$uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt_ref;
	my @ac_above_cutoff_above_qc_ref_gt_alt = @$ac_above_cutoff_above_qc_ref_gt_alt_ref;
	my @ac_above_cutoff_above_qc_ref_lt_alt = @$ac_above_cutoff_above_qc_ref_lt_alt_ref;
	my @auc_above_cutoff_above_qc_ref_gt_alt = @$auc_above_cutoff_above_qc_ref_gt_alt_ref;
	my @auc_above_cutoff_above_qc_ref_lt_alt = @$auc_above_cutoff_above_qc_ref_lt_alt_ref;

	#TODO clean up below:
	# table to check logic, make sure the numbers make sense
	my $allele_corrected_size = $#allele_corrected_pass + 1; my $allele_corrected_pass_dp_above_cutoff_size = keys(%ac_above_cutoff); my $allele_corrected_pass_dp_below_cutoff_size = keys(%ac_below_cutoff); my $allele_corrected_pass_not_in_allele_uncorrected = keys(%ac_missing_in_auc);
	my $ac_above_cutoff_above_qc_filter_size = keys(%ac_above_cutoff_above_qc_filter); my $ac_above_cutoff_below_qc_filter_size = keys(%ac_above_cutoff_below_qc_filter);
	my $ac_above_cutoff_above_qc_ref_gt_alt_counts = $#ac_above_cutoff_above_qc_ref_gt_alt +1;
	my $ac_above_cutoff_above_qc_ref_lt_alt_counts = $#ac_above_cutoff_above_qc_ref_lt_alt +1;

	print "\n\tAllele Corrected";
	print "\n\t\tSize of \"pass\" data sent to function:($allele_corrected_size)";
	print "\n\t\tSize of \"pass\" data above DP cutoff:($allele_corrected_pass_dp_above_cutoff_size)";
	print "\n\t\tSize of \"pass\" data below cutoff or its uncorrected mate below DP cutoff:($allele_corrected_pass_dp_below_cutoff_size), specific primers:\n";
	print "\t\t\t", join (", ", sort keys %ac_below_cutoff), "\n";
	print "\t\tSize of \"pass\" data (primer IDs) that are not in allele uncorrected:($allele_corrected_pass_not_in_allele_uncorrected), specific primers:\n";
	print "\t\t\t", join(", ", sort keys %ac_missing_in_auc), "\n";
	print "\t\tSize of data that's above the DP cutoff & above the QC pass filter:($ac_above_cutoff_above_qc_filter_size)\n";
	print "\t\tSize of data that's above the DP cutoff & below the QC pass filter:($ac_above_cutoff_below_qc_filter_size)\n";
	print "\t\tSize of data that has ref > alt:($ac_above_cutoff_above_qc_ref_gt_alt_counts)\n";
	print "\t\tSize of data that has ref < alt:($ac_above_cutoff_above_qc_ref_lt_alt_counts)\n\n";

	my $allele_uncorrected_size = $#allele_uncorrected_pass + 1; my $allele_uncorrected_pass_dp_above_cutoff_size = keys(%auc_above_cutoff); my $allele_uncorrected_pass_dp_below_cutoff_size = keys(%auc_below_cutoff); my $allele_uncorrected_pass_not_in_allele_corrected = keys(%auc_missing_in_ac);
	my $auc_above_cutoff_above_qc_filter_size = keys(%auc_above_cutoff_above_qc_filter); my $auc_above_cutoff_below_qc_filter_size = keys(%auc_above_cutoff_below_qc_filter);
	my $auc_above_cutoff_above_qc_ref_gt_alt_counts = $#auc_above_cutoff_above_qc_ref_gt_alt +1;
	my $auc_above_cutoff_above_qc_ref_lt_alt_counts = $#auc_above_cutoff_above_qc_ref_lt_alt +1;
	print "\tAllele Uncorrected\n";
	print "\t\tSize of \"pass\" data to function:($allele_uncorrected_size)\n";
	print "\t\tSize of \"pass\" data above DP cutoff:($allele_uncorrected_pass_dp_above_cutoff_size)\n";
	print "\t\tSize of \"pass\" data below DP cutoff or its uncorrected mate below DP cutoff:($allele_uncorrected_pass_dp_below_cutoff_size), specific primers:\n";
	print "\t\t\t", join(", ", sort keys %auc_below_cutoff), "\n";
	print "\t\tSize of \"pass\" data (primer IDs) that are not in allele corrected:($allele_uncorrected_pass_not_in_allele_corrected), specific primers:\n";
	print "\t\t\t", join(", ", sort keys %auc_missing_in_ac), "\n";
	print "\t\tSize of data that's above the DP cutoff & above the QC pass filter:($auc_above_cutoff_above_qc_filter_size)\n";
	print "\t\tSize of data that's above the DP cutoff & below the QC pass filter:($auc_above_cutoff_below_qc_filter_size)\n";
	print "\t\tSize of data that has ref > alt:($auc_above_cutoff_above_qc_ref_gt_alt_counts)\n";
	print "\t\tSize of data that has ref < alt:($auc_above_cutoff_above_qc_ref_lt_alt_counts)\n\n";


	my $remove_from_50nt_files_size = keys(%remove_from_50nt_files);
	print "\tPrimers to Remove\n";
	print "\t\tNumber of amplicons to remove:($remove_from_50nt_files_size), specific primers:\n";
	print "\t\t\t", join(", ", sort keys(%remove_from_50nt_files)),"\n\n";
	print "\t\t* Primers to remove includes:\n";
	print "\t\t\t- primers below cutoff AND/OR\n";
	print "\t\t\t- \"missing\" in at least one of allele corrected/uncorrected,\n";
	print "\t\t\t- AND/OR primers that do not match the design file (e.g. reference mismatch, or the chr + loci + primer ID combo not found in DF).\n\n";

	my $corrected_50nt_size = $#corrected_50nt + 1; my $corrected_50nt_above_cutoff_size = $#corrected_50nt_above_cutoff + 1;
	my $corrected_50nt_below_cutoff_size = $#corrected_50nt_below_cutoff + 1;
	my $corrected_50nt_above_cutoff_above_qc_filter_size = $#corrected_50nt_above_cutoff_above_qc_filter + 1; 
	my $corrected_50nt_above_cutoff_below_qc_filter_size = $#corrected_50nt_above_cutoff_below_qc_filter + 1; 
	my $corrected_50nt_above_cutoff_above_qc_ref_gt_alt_counts = $#corrected_50nt_above_cutoff_above_qc_ref_gt_alt +1;
	my $corrected_50nt_above_cutoff_above_qc_ref_lt_alt_counts = $#corrected_50nt_above_cutoff_above_qc_ref_lt_alt +1;
	print "\t50nt Error Corrected\n";
	print "\t\tSize of input data to function:($corrected_50nt_size)\n";
	print "\t\tSize of only data representing primers \"to keep\" based upon DP cutoff:($corrected_50nt_above_cutoff_size)\n";
	print "\t\tSize of only data representing primers \"to remove\" based upon DP cutoff:($corrected_50nt_below_cutoff_size)\n";
	print "\t\tSize of data above QC pass filter:($corrected_50nt_above_cutoff_above_qc_filter_size)\n";
	print "\t\tSize of data below QC pass filter:($corrected_50nt_above_cutoff_below_qc_filter_size)\n";
	print "\t\tSize of data that has ref > alt:($corrected_50nt_above_cutoff_above_qc_ref_gt_alt_counts)\n";
	print "\t\tSize of data that has ref < alt:($corrected_50nt_above_cutoff_above_qc_ref_lt_alt_counts)\n\n";

	my $uncorrected_50nt_size = $#uncorrected_50nt + 1; my $uncorrected_50nt_above_cutoff_size = $#uncorrected_50nt_above_cutoff + 1;
	my $uncorrected_50nt_below_cutoff_size = $#uncorrected_50nt_below_cutoff + 1; 
	my $uncorrected_50nt_above_cutoff_above_qc_filter_size = $#uncorrected_50nt_above_cutoff_above_qc_filter + 1; 
	my $uncorrected_50nt_above_cutoff_below_qc_filter_size = $#uncorrected_50nt_above_cutoff_below_qc_filter + 1; 
	my $uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt_counts = $#uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt +1;
	my $uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt_counts = $#uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt +1;
	print "\t50nt Uncorrected\n";
	print "\t\tSize of input data to function:($uncorrected_50nt_size)\n";
	print "\t\tSize of only data representing primers \"to keep\":($uncorrected_50nt_above_cutoff_size)\n";
	print "\t\tSize of only data representing primers \"to remove\":($uncorrected_50nt_below_cutoff_size)\n";
	print "\t\tSize of data above QC pass filter:($uncorrected_50nt_above_cutoff_above_qc_filter_size)\n";
	print "\t\tSize of data below QC pass filter:($uncorrected_50nt_above_cutoff_below_qc_filter_size)\n";
	print "\t\tSize of data that has ref > alt:($uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt_counts)\n";
	print "\t\tSize of data that has ref < alt:($uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt_counts)\n\n";

	# return data structures of interest, so we can print them out!
	# add in stuff from ref/alt check
	return (\%ac_above_cutoff, \%ac_below_cutoff, \%auc_above_cutoff, \%auc_below_cutoff, \@corrected_50nt_above_cutoff, \@corrected_50nt_below_cutoff, \@uncorrected_50nt_above_cutoff, \@uncorrected_50nt_below_cutoff, \%ac_missing_in_auc, \%auc_missing_in_ac, 
	\%ac_above_cutoff_above_qc_filter, \%ac_above_cutoff_below_qc_filter, \%auc_above_cutoff_above_qc_filter, \%auc_above_cutoff_below_qc_filter,
	\@corrected_50nt_above_cutoff_above_qc_filter, \@corrected_50nt_above_cutoff_below_qc_filter, \@uncorrected_50nt_above_cutoff_above_qc_filter, \@uncorrected_50nt_above_cutoff_below_qc_filter,
	\@corrected_50nt_above_cutoff_above_qc_ref_gt_alt, \@corrected_50nt_above_cutoff_above_qc_ref_lt_alt, \@uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt, \@uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt, \@ac_above_cutoff_above_qc_ref_gt_alt, \@ac_above_cutoff_above_qc_ref_lt_alt, \@auc_above_cutoff_above_qc_ref_gt_alt, \@auc_above_cutoff_above_qc_ref_lt_alt);

}

sub remove_data_if_ref_lt_alt{
	# de-reference the references that were passed to the subroutine
	my ($corrected_50nt_above_cutoff_above_qc_filter_ref, $uncorrected_50nt_above_cutoff_above_qc_filter_ref, $ac_above_cutoff_above_qc_filter_ref, $auc_above_cutoff_above_qc_filter_ref, $only_need_primer_in_raw_data, $debug) =@_;

	# de-reference the references that were passed to the subroutine
	my @corrected_50nt_above_cutoff_above_qc_filter = @$corrected_50nt_above_cutoff_above_qc_filter_ref;
	my @uncorrected_50nt_above_cutoff_above_qc_filter = @$uncorrected_50nt_above_cutoff_above_qc_filter_ref;
	my %ac_above_cutoff_above_qc_filter = %$ac_above_cutoff_above_qc_filter_ref;
	my %auc_above_cutoff_above_qc_filter = %$auc_above_cutoff_above_qc_filter_ref;

	# run through for array ones
	my (@corrected_50nt_above_cutoff_above_qc_ref_gt_alt, @corrected_50nt_above_cutoff_above_qc_ref_lt_alt, @uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt, @uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt, %ac_above_cutoff_above_qc_ref_gt_alt, %ac_above_cutoff_above_qc_ref_lt_alt, %auc_above_cutoff_above_qc_ref_gt_alt, %auc_above_cutoff_above_qc_ref_lt_alt);

	# run through for hash ones

	# allele corrected
	foreach my $key (keys %ac_above_cutoff_above_qc_filter) {
		my $line = $ac_above_cutoff_above_qc_filter{$key};
		my @temp = split("\t", $line);

		my $AD = $temp[6];

		my $ref_count = $1 if $AD =~m/AD=(\d*)\,/ ;
		my $alt_count = $1 if $AD =~m/AD=\d*\,(\d*)/;

		if($ref_count > $alt_count){
			$ac_above_cutoff_above_qc_ref_gt_alt{$key} = $line;
		}else{
			$ac_above_cutoff_above_qc_ref_lt_alt{$key} = $line;
		}
	}

	# allele uncorrected
	foreach my $key (keys %auc_above_cutoff_above_qc_filter) {
		my $line = $auc_above_cutoff_above_qc_filter{$key};
		my @temp = split("\t", $line);

		my $AD = $temp[6];

		my $ref_count = $1 if $AD =~m/AD=(\d*)\,/ ;
		my $alt_count = $1 if $AD =~m/AD=\d*\,(\d*)/;

		if(($only_need_primer_in_raw_data == 0) && ($ac_above_cutoff_above_qc_ref_lt_alt{$key})){
			$auc_above_cutoff_above_qc_ref_lt_alt{$key} = $line;
		}elsif($ref_count > $alt_count){
			$auc_above_cutoff_above_qc_ref_gt_alt{$key} = $line;
		}else{
			$auc_above_cutoff_above_qc_ref_lt_alt{$key} = $line;

			delete $ac_above_cutoff_above_qc_ref_gt_alt{$temp[37]};
			$ac_above_cutoff_above_qc_ref_lt_alt{$temp[37]} = $line;

		}
	}

	# keep allele corrected and uncorrected consistent (assuming no toggle for only_need_primer_in_raw_data)

	# keep allele corrected and 50nt consistent
	for my $line (@corrected_50nt_above_cutoff_above_qc_filter){
	 	my @temp = split("\t", $line);

	 	if (exists $ac_above_cutoff_above_qc_ref_lt_alt{$temp[37]}){
	 		push (@corrected_50nt_above_cutoff_above_qc_ref_lt_alt, $line);
	 	}else{
	 		push (@corrected_50nt_above_cutoff_above_qc_ref_gt_alt, $line);
	 	}
	 }

	# keep allele uncorrected and 50nt consistent
	for my $line (@uncorrected_50nt_above_cutoff_above_qc_filter){
	 	my @temp = split("\t", $line);

	 	if (exists $auc_above_cutoff_above_qc_ref_lt_alt{$temp[37]}){
	 		push (@uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt, $line);
	 	}else{
	 		push (@uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt, $line);
	 	}
	}


	# TODO: make sure that if data is missing in uncorrected, it is removed in corrected

	# switch from hashes to arrays
	my @ac_above_cutoff_above_qc_ref_gt_alt = sort(values(%ac_above_cutoff_above_qc_ref_gt_alt)); 
	my @ac_above_cutoff_above_qc_ref_lt_alt = sort(values(%ac_above_cutoff_above_qc_ref_lt_alt));
	my @auc_above_cutoff_above_qc_ref_gt_alt = sort(values(%auc_above_cutoff_above_qc_ref_gt_alt));
	my @auc_above_cutoff_above_qc_ref_lt_alt = sort(values(%auc_above_cutoff_above_qc_ref_lt_alt)); 

	return (\@corrected_50nt_above_cutoff_above_qc_ref_gt_alt, \@corrected_50nt_above_cutoff_above_qc_ref_lt_alt, \@uncorrected_50nt_above_cutoff_above_qc_ref_gt_alt, \@uncorrected_50nt_above_cutoff_above_qc_ref_lt_alt, \@ac_above_cutoff_above_qc_ref_gt_alt, \@ac_above_cutoff_above_qc_ref_lt_alt, \@auc_above_cutoff_above_qc_ref_gt_alt, \@auc_above_cutoff_above_qc_ref_lt_alt);
}




# for usage statement

__END__

=head1 NAME

Script to verify detected primary alternate allele matches that of the oligo target file.

=head1 SYNOPSIS

verify_primary_alt_allele.pl [options]

  Options:
  	-h|help                 		brief help message
  	-oligo_target_file      		file containing list of oligos targeting which sites, tab-delimited
  	-allele_corrected       		allele file (error corrected), tab-delimited
  	-allele_uncorrected     		allele file (not error corrected), tab-delimited
	-corrected_50nt         		50nt file (error corrected), tab-delimited
	-uncorrected_50nt       		50nt file (not error corrected), tab-delimited
	-output_dir             		folder to save output files
  	-sample_name            		sample or chip name, which is used as a prefix in the output filenames
	-cutoff                			minimum read depth to take (default: 0)
	-qc_pass_filter         		minimum fraction of passable reads from total reads (default: .7)
	-remove_ref_lt_alt      		remove data where the reference counts are less than alt counts
	-only_need_primer_in_raw_data		if toggled, data will be retained even if that primer is found only in raw data.
						by default, primer must be found in allele corrected and uncorrected to be kept.
  	-debug                  		turn on debugging output
=cut
