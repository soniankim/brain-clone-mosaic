#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;	# for script arguments for file input passing
use Pod::Usage;   	# for usage statement

# command line arg handling
my ($help, $oligo_target_file, $allele_corrected, $allele_uncorrected, $corrected_50nt, $uncorrected_50nt, $output_dir, $sample_name, $cutoff, $debug);

GetOptions(
        'help|h'               => \$help,
        'oligo_target_file=s'  => \$oligo_target_file,
        'allele_corrected=s'   => \$allele_corrected,
        'allele_uncorrected=s' => \$allele_uncorrected,
        'corrected_50nt=s'     => \$corrected_50nt,
        'uncorrected_50nt=s'   => \$uncorrected_50nt,
        'output_dir=s'         => \$output_dir,
        'sample_name=s'        => \$sample_name,
        'cutoff=i'             => \$cutoff,
        'debug'                => \$debug
);

# debug flag is for printing
$debug = 0 unless ($debug);

# if the user running the script uses the help option or forgets to
# pass files, print the usage statement
pod2usage(-verbose =>1) if $help;
pod2usage("$0: Not enough files provided.\n") unless $oligo_target_file && $allele_corrected && $allele_uncorrected && $output_dir && $sample_name && $corrected_50nt && $uncorrected_50nt;

# unless you gave a cutoff, we'll use 0
$cutoff = 0 unless $cutoff;

# output directory
if (-d $output_dir){
	die "Please select a new output directory name. $output_dir already exists.\n";
}else{
	mkdir $output_dir unless -d $output_dir;
}

# declaring hash and arrays for storing data that was read in
my (%oligo_target, @allele_corrected, @allele_uncorrected, @corrected_50nt, @uncorrected_50nt);

# Tell the user what files we are using
print "\nRunning verify_primary_alt_allele.pl with the following parameters:\n\n";
printf("%-25s %-150s\n",   "Oligo design file:",       $oligo_target_file);
printf("%-25s %-150s\n",   "Allele corrected file:",   $allele_corrected); 
printf("%-25s %-150s\n",   "Allele uncorrected file:", $allele_uncorrected);
printf("%-25s %-150s\n",   "Corrected 50nt file:",     $corrected_50nt);
printf("%-25s %-150s\n",   "Uncorrected 50nt file:",   $uncorrected_50nt);
printf("%-25s %-150s\n",   "Output directory:",        $output_dir);
printf("%-25s %-150s\n",   "Cutoff:",                  $cutoff);
printf("%-25s %-150s\n\n", "Sample name:",             $sample_name);

# read in files and save to appropriate data structures
%oligo_target       = read_file($oligo_target_file,  1);
@allele_corrected   = read_file($allele_corrected,   0);
@allele_uncorrected = read_file($allele_uncorrected, 0);
@corrected_50nt     = read_file($corrected_50nt,     0);
@uncorrected_50nt   = read_file($uncorrected_50nt,   0);

# compare each file (allele corrected/uncorrected & 50nt corrected/uncorrected)
# with the design file.
my (@allele_corrected_pass, @allele_corrected_fail, @allele_uncorrected_pass, @allele_uncorrected_fail); 
#TODO not sure if we need these yet
my (@corrected_50nt_pass, @corrected_50nt_fail, @uncorrected_50nt_pass, @uncorrected_50nt_fail);

#TODO add in 50nt files to pass to this filtering sub
my ($allele_corrected_pass_ref, $allele_corrected_fail_ref)  = filter(\%oligo_target, \@allele_corrected, $debug);
@allele_corrected_pass = @$allele_corrected_pass_ref; @allele_corrected_fail = @$allele_corrected_fail_ref;

#TODO add in 50nt files to pass to this filtering sub
my ($allele_uncorrected_pass_ref, $allele_uncorrected_fail_ref) = filter(\%oligo_target, \@allele_uncorrected, $debug);
@allele_uncorrected_pass = @$allele_uncorrected_pass_ref; @allele_uncorrected_fail = @$allele_uncorrected_fail_ref;


print "List of output files:\n";
# write the files
# need output filename (using basename), array, and output directory
write_file(\@allele_corrected_pass,   $output_dir, $sample_name, "alleleErrC.pass.txt");
write_file(\@allele_corrected_fail,   $output_dir, $sample_name, "alleleErrC.fail.txt");
write_file(\@allele_uncorrected_pass, $output_dir, $sample_name, "alleleNoC.pass.txt");
write_file(\@allele_uncorrected_fail, $output_dir, $sample_name, "alleleNoC.fail.txt");
#TODO add in output files for 50nt corr/uncorr pass/fail

# read file subroutine to injest input files
sub read_file{
	my ($file, $is_this_oligo_target_file) = @_;
	my (%design_file, @file_contents);

	open (my $read_fh, "<", $file) || die "Cannot open $file for reading.\n";  
	
	while (my $line =<$read_fh>){
		chomp $line;
		
		# DESIGN FILE, columns of interest (indexing by 0 here):
		# 2: Chromosome
		# 3: Loci
		# 5: Ref
		# 6: Alt
		if ($is_this_oligo_target_file){
			my @temp = split(/\t/, $line);
			
			my $chr = $temp[2]; my $loci = $temp[3]; my $ref = $temp[5]; my $alt = $temp[6];
			my $key = $temp[2] ."_" . $loci;
		
			# if we already have this loci (chromosome + loci combo stored in $key) in the %design_file hash,
			# check for errors. Errors include same chr+loci, but different ref alleles, OR same chr+loci, but different primary alt alleles.
			# if there are not any errors, we continue to adding it to the hash - which we really don't need to do
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
			push(@file_contents, $line);
		}
	
	}
	
	return %design_file if %design_file;
	return @file_contents if @file_contents;
	close $read_fh;
}

# write file subroutine, to create output files in the user-specified output directory
sub write_file{
	my ($data_to_write_ref, $dir, $prefix, $suffix) = @_;
	my @data_to_write = @$data_to_write_ref;

	my $full_path_to_output_file = $dir . "/" . $prefix . "-" . $suffix;
	print $full_path_to_output_file, "\n";
	open(my $write_fh, ">", $full_path_to_output_file) || die "Unable to open filehandle for writing $dir\\$prefix$suffix:$!\n";
	
	for my $line (@data_to_write){
		print $write_fh $line, "\n";
	}

	close $write_fh;
}

#TODO: need to add in new function for checking alt allele, and rearrange data if necessary.



#TODO: need to rewrite this filtering function. Will repurpose the name of this function for filtering for passing read depth.
#TODO: if corrected/uncorrected read depth total does not match -> throw out amplicon
#TODO: if amplicon is thrown out, discard 50nt entries (NOT VARIANT)

# Function for comparison of the design file & a given file, which is one of:
# allele_corrected, allele_uncorrected, corrected_50nt, or uncorrected_50nt.
# We need to see if the primary alternate allele that's in a given file matches 
# what we'd expect (e.g. what the design file contains). 

# renamed this sub from "check_primary_alt_allele" to "filter" as this gives a better indication of what it actually does.
# After speaking with Sonia, the script should not throw out lines with REF only called, should not throw out anything with
# the expected alt allele & should only run a check when the secondary alt allele is called.
sub filter{
	my ($design_file_ref, $file_to_compare_ref, $print_if_debug) = @_;

	# de-reference the references that were passed to the subroutine
	my %design_file = %$design_file_ref; 
	my @file_to_compare = @$file_to_compare_ref; 
	
	# examine the file of interest (allele or 50nt, and corrected or uncorrected)
	# see what matched up with the design file and what does not.
	my (@pass, @fail); 
	foreach my $line (@file_to_compare){
		my @temp = split(/\t/, $line);
		# ALLELE CORRECTED\UNCORRECTED OR 50nt CORRECTED\UNCORRECTED, columns of interest
		# (indexed starting at 0)
		# 0: chromosome
		# 1: loci
		# 3: reference allele
		# 4: alt allele
		my $chr = $temp[0]; my $loci = $temp[1]; my $ref = $temp[3]; my $alt = $temp[4];

		# check that we do have this chromosome + loci combo in the design file
		if ($design_file{$chr."_".$loci}){
			# check reference and alt.
			
			# going to extract the first character from the alt string to get the primary alt allele.
			# alt string could be only one character, such as A
			# it also could be multiple comma separated characters, such as T,C,G
			# and if the primary alt allele is actually the reference, it will start with <*>.

			# irregardless of what the string looks like, if the extracted first character does not match
			# the design file's alt, then we do not want to retain it.

			# update: this is not the intended behavior of the script, after speaking with Sonia.
			# here are the expected outcomes for each scenario:
			# * anything reference only --> pass
			# * anything with expected alt allele --> pass
			# * check should only be done when secondary alt allele is called.


			my $primary_alt = substr($alt, 0, 1); 
			#print length($alt), "\t", $alt, "\n";
			print "Design file ref and alt: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $design_file{$chr."_".$loci}->{'alt'}, "\n" if $print_if_debug;
			print "Input file ref and alt: ", $ref, "\t", $alt, "\n" if $print_if_debug; 

			# check first that reference only allele goes to pass.  
			if ( ($design_file{$chr."_".$loci}->{'ref'} eq $ref)  &&  $alt eq "<*>"){
				push(@pass, $line);
				print "PASS due to reference match: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\t", "and alt is $alt\n\n" if $print_if_debug;
			}
			# next check that expected alt allele goes to pass. 
			elsif ( ($design_file{$chr."_".$loci}->{'ref'} eq $ref ) && ( $design_file{$chr."_".$loci}->{'alt'} eq $primary_alt )){
				push(@pass, $line); 
				print "PASS due to primary alt match: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\t", "and\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$primary_alt\t $alt\n\n" if $print_if_debug; 
			}	
			# if reference/alt didn't match, we do not want to retain the line..
			# need both reference and alt to match the design file.
			else{
				push (@fail, $line);
				print "FAIL due to: ",  $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\tand\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$primary_alt\t $alt\n\n" if $print_if_debug; 
				next;  
			}	

		}else{
			# warn if we do not find a matching chr + loci combination in the design file. 
			print "Chromosome and loci combination $chr $loci not found in design file.\n"; 
		}			
	}

	return (\@pass, \@fail);
}


# for usage statement

__END__

=head1 NAME

Script to verify detected primary alternate allele matches that of the oligo target file.

=head1 SYNOPSIS

verify_primary_alt_allele.pl [options]

  Options: 
  	-h|help                 brief help message
  	-oligo_target_file      file containing list of oligos targeting which sites, tab-delimited
  	-allele_corrected       allele file (error corrected), tab-delimited
  	-allele_uncorrected     allele file (not error corrected), tab-delimited
        -corrected_50nt         50nt file (error corrected), tab-delimited
        -uncorrected_50nt       50nt file (not error corrected), tab-delimited
	-output_dir             folder to save output files
  	-sample_name            sample or chip name, which is used as a prefix in the output filenames
        -cutoff                 minimum read depth to take (default: 0)
  	-debug                  turn on debugging output
=cut
