#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;	# for script arguments for file input passing
use Pod::Usage;   	# for usage statement

my ($help, $oligo_target_file, $allele_corrected, $allele_uncorrected, $corrected_50nt, $uncorrected_50nt, $output_dir, $sample_name, $cutoff, $debug);

# command line arg handling
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

# debug flag will print out lots of info so you can verify the script is doing what we expect it to do
$debug = 0 unless ($debug);

# if the user running the script uses the help option or forgets to
# pass files, print the usage statement
pod2usage(-verbose =>1) if $help;
pod2usage("$0: Not enough files provided.\n") unless $oligo_target_file && $allele_corrected && $allele_uncorrected && $output_dir && $sample_name && $corrected_50nt && $uncorrected_50nt;

# unless you gave a cutoff, we'll use 0
$cutoff = 0 unless $cutoff;

# output directory name must be unique
if (-d $output_dir){
	die "Please select a new output directory name. $output_dir already exists.\n";
}else{
	mkdir $output_dir unless -d $output_dir;
}

# declaring hash and arrays for storing data from each of the files that will be read in
my (%oligo_target, @allele_corrected, @allele_uncorrected, @corrected_50nt, @uncorrected_50nt);

# Tell the user what files & parameters we are using
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
%oligo_target       = read_file($oligo_target_file,  1); # use hash as this data will be searched, everything else is in arrays.
@allele_corrected   = read_file($allele_corrected,   0);
@allele_uncorrected = read_file($allele_uncorrected, 0);
@corrected_50nt     = read_file($corrected_50nt,     0);
@uncorrected_50nt   = read_file($uncorrected_50nt,   0);

# compare each file (allele corrected/uncorrected & 50nt corrected/uncorrected)
# with the design file.
my (@allele_corrected_pass, @allele_corrected_fail, @allele_uncorrected_pass, @allele_uncorrected_fail);

# send allele_corrected file to alternate allele check function
print "Running alt allele check for allele corrected file.\n";
my ($allele_corrected_pass_ref, $allele_corrected_fail_ref)  = alt_allele_check(\%oligo_target, \@allele_corrected, $debug);
@allele_corrected_pass = @$allele_corrected_pass_ref; @allele_corrected_fail = @$allele_corrected_fail_ref;

# send allele uncorrected file to alternate allele check function
print "Running alt allele check for allele uncorrected file.\n";
my ($allele_uncorrected_pass_ref, $allele_uncorrected_fail_ref) = alt_allele_check(\%oligo_target, \@allele_uncorrected, $debug);
@allele_uncorrected_pass = @$allele_uncorrected_pass_ref; @allele_uncorrected_fail = @$allele_uncorrected_fail_ref;

# send allele corrected + uncorrected data from the alt_allele_check function to the filtering function
# also send the 50 nt files to the filtering function
my (@allele_corrected_pass_filtered,  @allele_uncorrected_pass_filtered, @corrected_50nt_filtered,  @uncorrected_50nt_filtered);
#my ($allele_corrected_pass_filtered_ref, @allele_uncorrected_pass_filtered_ref, @corrected_50nt_filtered_ref, @uncorrected_50nt_filtered_ref) = filter(\@allele_corrected_pass, \@allele_uncorrected_pass, \@corrected_50nt, \@uncorrected_50nt, $cutoff, $debug);
#TODO save returned >10k files from filter function
print "Running filter function to discard allele corrected, allele uncorrected, 50nt corrected, and 50nt uncorrected data below the cutoff.\n";
filter(\@allele_corrected_pass, \@allele_uncorrected_pass, \@corrected_50nt, \@uncorrected_50nt, $cutoff, $debug);
#@allele_corrected_pass_filtered = @$allele_corrected_pass_filtered_ref; @allele_uncorrected_pass_filtered = @$allele_uncorrected_pass_filtered_ref;
#@corrected_50nt_filtered = @$corrected_50nt_filtered_ref; @uncorrected_50nt_filtered =@$uncorrected_50nt_filtered_ref;


print "List of output files:\n";
# TODO: change this to the data returned from the filter sub
# write the files
# need output filename (using basename), array, and output directory
write_file(\@allele_corrected_pass,   $output_dir, $sample_name, "alleleErrC.pass.txt");
#write_file(\@allele_corrected_fail,   $output_dir, $sample_name, "alleleErrC.fail.txt");
write_file(\@allele_uncorrected_pass, $output_dir, $sample_name, "alleleNoC.pass.txt");
#write_file(\@allele_uncorrected_fail, $output_dir, $sample_name, "alleleNoC.fail.txt");
#TODO add in output files for 50nt corr/uncorr pass/fail

# read file subroutine to ingest input files
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

		# check that we do have this chromosome + loci combo in the design file
		if ($design_file{$chr."_".$loci}){
			# going to extract the first character from the alt string to get the primary alt allele.
			# alt column could be only one character, such as A
			# it also could be multiple comma separated characters, such as T,C,G
			# and if there is no alt allele - the column denotes the reference, represented as <*>.

			my $primary_alt = substr($alt, 0, 1);
			#print "Design file ref and alt: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $design_file{$chr."_".$loci}->{'alt'}, "\n" if $print_if_debug;
			#print "Input file ref and alt: ", $ref, "\t", $alt, "\n" if $print_if_debug;
			#print $line,"\n" if $print_if_debug;

			# verify reference in this line is the same as in the design file:
			if($design_file{$chr."_".$loci}->{'ref'} eq $ref){

				# ALT ALLELE CHECK CONDITIONS:
				# 1. ref allele only
				if ( $alt eq "<*>" ){
					push(@pass, $line);
					#print "PASS due to reference match: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\t", "and alt is $alt\n\n" if $print_if_debug;

				# 2. expected alt allele in primary allele position
				}elsif ( $design_file{$chr."_".$loci}->{'alt'} eq $primary_alt ){
					push(@pass, $line);
					#print "PASS due to primary alt match: ", $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\t", "and\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$primary_alt\t $alt\n\n" if $print_if_debug;

				# 3. expected alt allele in secondary/tertiary/etc position
				}elsif ( $alt =~m/($design_file{$chr."_".$loci}->{'alt'})/ ){
					#TODO: need to modify:
					# ALT (col 8)
					# TOTAL DP (col 9)  = ref + ALT
					# AAF (col 10)  = new ALT / new total DP = ((col 8)/ (col9))
					# % nt QC (col 11) = new col 9 / col 6
					# position of match in $alt
					my $alt_tmp = $alt;
					$alt_tmp =~ tr/,//d; # remove commas to get position
							     # used to pull out correct alt counts
					my $position = index($alt_tmp, $design_file{$chr."_".$loci}->{'alt'}) ;

					my $AD = $temp_line[6];
					$AD =~ tr/AD=//d;
					my @temp_for_alt_counts =split(",", $AD);

					# pull out correct alt counts based upon the position of the expected alt allele in col 8. So in the AD column, the correct alt counts will be position plus one (as reference is listed first in the AD string)
					my $correct_alt_counts = $temp_for_alt_counts[$position + 1];
					my $recalc_dp = $temp_for_alt_counts[0] + $correct_alt_counts; # ref + alt
					my $aaf = ($correct_alt_counts/ $recalc_dp);
					my $percent_nt_qc = ($recalc_dp / $temp_line[5]);

					my @end_of_line = @temp_line[11...38];
					my $end_of_line = join("\t", @end_of_line) ;
					my $modified_line = "$chr\t$loci\t$temp_line[2]\t$ref\t$alt\t$temp_line[5]\tAD=$AD\t$correct_alt_counts\t$recalc_dp\t$aaf\t$percent_nt_qc\t$end_of_line\n";

					#print $modified_line, "\n\n" if $print_if_debug;
					push(@pass, $modified_line);
					#print "PASS due to: ",  $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\tand\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$alt\n" if $print_if_debug;

				# check for alternative alleles that are not the expected one (e.g. not in the design file)
				# modify alt, total dp, aaf if so

				# 4. alternative alleles that are not the expected one (e.g. not in the design file)
				}elsif( $alt !~/$design_file{$chr."_".$loci}->{'alt'}/ ){
					#print "PASS (needs modifications) due to: ",  $design_file{$chr."_".$loci}->{'ref'}, "\t", $ref, "\tand\t", $design_file{$chr."_".$loci}->{'alt'}, "\t$alt\n" if $print_if_debug;
					# Alt goes to 0
					# Total dp is (original total dp – original alt) or AD=(probably easier)
					# AAF goes to 0
					# % reads that are good quality (%nt QC), put this to 0??? Not sure yet
					my $AD = $temp_line[6];
					my $ref_count_from_AD = $1 if $AD =~m/AD=(\d+),/;
					my $percent_nt_qc = $ref_count_from_AD / $temp_line[5];
					my @end_of_line = @temp_line[11...38];
					my $end_of_line = join("\t", @end_of_line);
					my $modified_line = "$chr\t$loci\t$temp_line[2]\t$ref\t$alt\t$temp_line[5]\t$AD\t0\t$ref_count_from_AD\t0\t$temp_line[10]\t$end_of_line\n";
					#print $modified_line,"\n\n" if $print_if_debug;

					# save modified line
					push(@pass, $modified_line);

				# catch for any other condition we were not expecting
				}else{
					print "Exception in logic found in alt_allele_check function, due to the following line!\n$line\n";
					push(@fail, $line);
					print "FAIL due to a condition that is not covered by the logic in the alt_allele_check function.\n" if $print_if_debug;
				}


			# reference in this line IS NOT the same as in the design file, yield an error:
			}else{
				print "Reference listed for $chr $loci in the following line does not match the ref in the design file.\n";
				push (@fail, $line);
				print "FAIL due to different ref in design file for $chr $loci. DF has  ", $design_file{$chr."_".$loci}->{'ref'}, "  and this line has $ref\n" if $print_if_debug;
			}


		}else{
			# warn if we do not find a matching chr + loci combination in the design file.
			print "Chromosome and loci combination $chr $loci not found in design file.\n";
			push(@fail, $line);
			print "FAIL due to: $chr & $loci missing in design file.\n" if $print_if_debug;
		}
	} #foreach close
	my $input_array_size = $#file_to_compare + 1; my $pass_array_size = $#pass +1; my $fail_array_size = $#fail +1;
	print "\tSize of input data to function:($input_array_size)\n\tSize of output \"passing\" data:($pass_array_size)\n\tSize of output \"failing\" data:($fail_array_size)\n\n";
	return (\@pass, \@fail);
}


#TODO: need to rewrite this filtering function. Will repurpose the name of this function for filtering for passing read depth.
#TODO: if corrected/uncorrected read depth total does not match -> throw out amplicon
#TODO: if amplicon is thrown out, discard 50nt entries (NOT VARIANT)

sub filter{
	my ($allele_corrected_pass_ref, $allele_uncorrected_pass_ref, $corrected_50nt_ref, $uncorrected_50nt_ref, $cutoff, $print_if_debug) = @_;

	# de-reference the references that were passed to the subroutine
	my @allele_corrected_pass = @$allele_corrected_pass_ref;
	my @allele_uncorrected_pass = @$allele_uncorrected_pass_ref;
	my @corrected_50nt = @$corrected_50nt_ref;
	my @uncorrected_50nt = @$uncorrected_50nt_ref;

	# process allele corrected "pass" data and segment data into >= DP cutoff, and < DP cutoff
	my %ac_above_cutoff; my %ac_below_cutoff;

	# $temp_line[8] - total DP
	# $temp_line[37] - primer ID (three digit code)
	foreach my $line (@allele_corrected_pass){
		my @temp_line = split("\t", $line);
		if ($temp_line[8] >= $cutoff){
			$ac_above_cutoff{$temp_line[37]}=$line;
		}elsif($temp_line[8] < $cutoff){
			$ac_below_cutoff{$temp_line[37]}=$line;
		}else{
			print "Error with allele corrected pass data!:\n$line\n";
		}
	}

	# go through allele uncorrected pass
	# for the discarded ones in allele corrected pass -> discard from allele uncorrected too
	# if any discarded from allele uncorrected (<cutoff) -> discard from allele corrected, too
	# then for the discarded ones, remove 50nt as well

	my %auc_above_cutoff; my %auc_below_cutoff; my %auc_missing_in_ac;
	my @remove_from_50nt_files; # this means the amplicon was < cutoff in at least one of allele corrected/uncorrected

	foreach my $line (@allele_uncorrected_pass){
		my @temp_line = split("\t", $line);
		#print "$line\t$temp_line[8]\n" if $print_if_debug;

		# if above cutoff, must check to see if allele_corrected is also above the cutoff for this amplicon
		if ( $temp_line[8]>=$cutoff  ){

			# allele uncorrected above cutoff, allele corrected below cutoff --> DISCARD
			if ($ac_below_cutoff{$temp_line[37]}){
				# add to remove_from_50nt_files
				push (@remove_from_50nt_files, $temp_line[37]);
				$auc_below_cutoff{$temp_line[37]}=$line;
			# both above cutoff ---> KEEP
			}elsif($ac_above_cutoff{$temp_line[37]}){
				$auc_above_cutoff{$temp_line[37]}=$line;
			# allele uncorrected primer ID is not in allele corrected---> PUT IN SEPARATE BUCKET
			}else{
				print "\tError in filter function. Missing primer $temp_line[37] for allele corrected\n" unless $ac_above_cutoff{$temp_line[37]} || $ac_below_cutoff{$temp_line[37]};
				print "Unspecified error in filter function" if $ac_above_cutoff{$temp_line[37]} || $ac_below_cutoff{$temp_line[37]};
				$auc_missing_in_ac{$temp_line[37]} = $line;
			}

		# allele uncorrected below cutoff:
		# check to see if allele_corrected is above the cutoff.
		}else{

			# allele uncorrected below cutoff, allele corrected below cutoff --> DISCARD
			if ($ac_below_cutoff{$temp_line[37]}){
				push(@remove_from_50nt_files, $temp_line[37]);
				$auc_below_cutoff{$temp_line[37]}=$line;
			# allele uncorrected below, allele corrected above cutoff ---> DISCARD
			}elsif($ac_above_cutoff{$temp_line[37]}){
				push(@remove_from_50nt_files, $temp_line[37]);
				# remove from above cutoff array for allele corrected, as the allele uncorrected is below
				my $to_keep = $ac_above_cutoff{$temp_line[37]};
				delete $ac_above_cutoff{$temp_line[37]};
				$ac_below_cutoff{$temp_line[37]} = $to_keep;
				$auc_below_cutoff{$temp_line[37]}=$line;

			# allele uncorrected primer ID is not in allele corrected---> PUT IN SEPARATE BUCKET
			}else{
				print "\tError in filter function. Missing primer $temp_line[37] for allele corrected\n" unless $ac_below_cutoff{$temp_line[37]} || $ac_above_cutoff{$temp_line[37]};
				print "Unspecified error in filter function" if $ac_below_cutoff{$temp_line[37]} || $ac_above_cutoff{$temp_line[37]};
				$auc_missing_in_ac{$temp_line[37]} = $line;
			}

		}

	}

	#TODO: check for primers that are in allele corrected data but not in allele uncorrected
	# shunt these to a bucket and print 'em out later
	my %ac_missing_in_auc;
	for my $primer (sort keys %ac_below_cutoff){
		if (! $auc_below_cutoff{$primer}){
			print "\tError in filter function. Missing primer $primer for allele uncorrected\n";
			$ac_missing_in_auc{$primer} = $ac_below_cutoff{$primer};
			delete $ac_below_cutoff{$primer};
		}
	}	 


	#TODO: remove the irrelevant lines from 50nt files
	my @corrected_50nt_above_cutoff; my @corrected_50nt_below_cutoff;
	my @uncorrected_50nt_above_cutoff; my @uncorrected_50nt_below_cutoff;

	# iterate through 50nt files for corrected and uncorrected
	# remove lines that we discarded from ac and auc
	# use @remove_from_50nt_files
	#TODO TODO TODO


	# table to check logic, make sure the numbers make sense
	my $allele_corrected_size = $#allele_corrected_pass + 1; my $allele_corrected_pass_dp_above_cutoff_size = keys(%ac_above_cutoff); my $allele_corrected_pass_dp_below_cutoff_size = keys(%ac_below_cutoff); my $allele_corrected_pass_not_in_allele_uncorrected = keys(%ac_missing_in_auc);
	print "\n\tSize of input allele corrected \"pass\" data to function:($allele_corrected_size)";
	print "\n\tSize of allele corrected \"pass\" data above cutoff:($allele_corrected_pass_dp_above_cutoff_size)";
	print "\n\tSize of allele corrected \"pass\" data below cutoff:($allele_corrected_pass_dp_below_cutoff_size)\n";
	print "\t\t";
	for my $ac_below (sort keys %ac_below_cutoff){
		print "$ac_below";
		print ", ";
	}
	print "\n\tSize of allele corrected \"pass\" data (primer IDs) that are not in allele uncorrected:($allele_corrected_pass_not_in_allele_uncorrected)\n";
	print "\t\t";
	for my $missing_primer (sort keys %ac_missing_in_auc){
		print "$missing_primer";
		print ", ";
	}
	print "\n\n";
	
	
	my $allele_uncorrected_size = $#allele_uncorrected_pass + 1; my $allele_uncorrected_pass_dp_above_cutoff_size = keys(%auc_above_cutoff); my $allele_uncorrected_pass_dp_below_cutoff_size = keys(%auc_below_cutoff); my $allele_uncorrected_pass_not_in_allele_corrected = keys(%auc_missing_in_ac);
	print "\tSize of input allele uncorrected \"pass\" data to function:($allele_uncorrected_size)\n";
	print "\tSize of allele uncorrected \"pass\" data above cutoff:($allele_uncorrected_pass_dp_above_cutoff_size)\n";
	print "\tSize of allele uncorrected \"pass\" data below cutoff:($allele_uncorrected_pass_dp_below_cutoff_size)\n";
	print "\t\t";
	for my $below_cutoff (sort keys %auc_below_cutoff){
		print "$below_cutoff";
		print ", ";
	}
	print "\n";
	print "\tSize of allele uncorrected \"pass\" data (primer IDs) that are not in allele corrected:($allele_uncorrected_pass_not_in_allele_corrected)\n";
	print "\t\t";
	for my $missing_primer (sort keys %auc_missing_in_ac){
		print "$missing_primer";
		print ", ";
	}
	print "\n\n";
	my $remove_from_50nt_files_size = $#remove_from_50nt_files +1;
	print "\tNumber of amplicons to remove:($remove_from_50nt_files_size)\n\n";
	my $corrected_50nt_size = $#corrected_50nt + 1;
	print "\tSize of input corrected 50nt data to function:($corrected_50nt_size)\n\n";

	my $uncorrected_50nt_size = $#uncorrected_50nt + 1;
	print "\tSize of input corrected 50nt data to function:($uncorrected_50nt_size)\n\n";

	# TODO:
	# return data structures of interest, so we can print them out!
	# @corrected_50nt_above_cutoff, @uncorrected_50nt_above_cutoff
	# %allele_uncorrected_pass_dp_above_cutoff, %allele_corrected_pass_dp_above_cutoff
	#return ();

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
