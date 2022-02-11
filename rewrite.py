#!/usr/bin/env python3


import argparse
import csv
import os
import pandas as pd
from pathlib import Path


def read_file(filename, hash=False):
	"""reads input file into structure (either csv object (or maybe pandas)
       or dictionary if design file
	"""

	temp_struct = ""
	if hash:
		temp_struct = {}
		with open(filename, 'r') as f:
			for line in f:
				fields = line.split('\t')
				# build key from Chromosome #, Loci, Primer
				key = f"{fields[2]}_{fields[3]}_{fields[0]}"
				print(f"{fields[2]}_{fields[3]}_{fields[0]}")
				temp_struct[key] = {}
				temp_struct[key]["ref"] = fields[5]
				temp_struct[key]["alt"] = fields[6]
				temp_struct[key]["primer_code"] = fields[0]
				temp_struct[key]["primer_name"] = fields[1]
	else:
		# need to remember that the input files should not have placeholders or have duplicate columns
		temp_struct = pd.read_csv(filename, delimiter='\t', header=None)
		#TODO: come to consensus on col names and verify them with Sonia
		temp_struct.columns = ['chromosome', 'site1', 'site2', 'ref', 'alt', 'total_reads', 'ad', 'altreads', 'dpreads', 'aaf', 'pass_read_fraction', 'chip_number', 'chip_number-sample_name', 'demultiplexed', 'is_it_corrected', 'downsampling', 'mutation_ID', 'fastq_name', 'primer_in_batch', 'sample_name']			

	return(temp_struct)


def alt_allele_check(design_file, allele_object_to_check):
	"""
	checks that don't result in modification of data:
		- ref allele was called
		- expected alt allele was called

	checks that result in modification of data:
		- expected alt allele called at non-primary position
		- non-expected allele was called
	
	fields to modify under the above illustrated conditions:
		- ALT, total DP, AAF, %nt QC
	"""

	for row in allele_object_to_check.index:
		# need chr, loci (site1), ref, alt, primer (primer_in_batch)

		# verify with design file
		# print(f"{allele_object_to_check['chromosome'][row]}_{allele_object_to_check['site1'][row]}_{allele_object_to_check['primer_in_batch'][row]}")

		# chromosome_site1_primer which is the unique key in design_file
		to_check =  f"{allele_object_to_check['chromosome'][row]}_{allele_object_to_check['site1'][row]}_{allele_object_to_check['primer_in_batch'][row]}"

		# in the design file
		if to_check in design_file:
			#print("The following data exists in the design file:", allele_object_to_check.loc[row,:])

			# need first character of 'alt' from allele_object_to_check. This'll be the primary alt: 
			primary_alt = allele_object_to_check['alt'][row][0]
			reference_in_allele_file = allele_object_to_check['ref'][row]
			allele_field = allele_object_to_check['alt'][row]

			# compare primary alt to reference
			reference_in_design_file =  design_file[to_check]['ref']
			expected_alt = design_file[to_check]['alt']

			# verify reference in this line is the same as in the design file:
			if reference_in_allele_file == reference_in_design_file:

				print(primary_alt, " ", expected_alt, " ", reference_in_design_file) #TODO: remove after debugging

				# string match to see if we have the expected alt allele
				# (relevant when identifying expected alt allele in non-primary position)
				pos_of_expected_alt_match = allele_field.find(expected_alt)

				# check if ref was called
				if allele_field == "<*>":
					print("reference match") #TODO: remove after debugging. Not tested yet
					#TODO: save to some data structure

				# was the expected alt allele called?
				elif primary_alt == expected_alt:
					print("primary alt match") #TODO: remove after debugging. Not tested yet
					#TODO: save to some data structure

				# was the expected allele called at the non-primary position?
				elif pos_of_expected_alt_match != -1:
					print("expected alt match in non-primary position") #TODO: remove after debugging. Not tested yet
					#TODO: need to modify altreads, dpreads, aaf, pass_read_fraction, and position of match in alt
					#TODO: save to some data structure
				
				# was a non-expected alt allele called?
				elif pos_of_expected_alt_match == -1:
					print("non-expected alt allele called") #TODO: remove after debugging. Not tested yet
					#TODO: need to modify altreads, dpreads, aaf, pass_read_fraction, and position of match in alt
					#TODO: save to some data structure

				else: 
					print("ERROR: Exception in logic found in alt_allele_check function due to", allele_object_to_check.loc[row,:])

			# reference is NOT the same as in the design file, which shouldn't ever happen.
			else:
				print("ERROR: Reference listed for chromosome, loci, primer ", to_check, "does not match the reference ", reference_in_design_file, "in the design file.")
				#TODO: need to verify this prints appropriately.

		# not in the design file, which shouldn't ever happen.	
		else:
			print("ERROR: The following data is missing in the design file:", allele_object_to_check.loc[row,:])
			#TODO: need to verify this prints appropriately.


if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()

	parser.add_argument("--oligo_target_file", 
							required=True, 
							help="file containing list of oligos targeting which sites, tab-delimited")
	parser.add_argument("--corrected_allele_file", 
							required=True, 
							help="allele file (error corrected), tab-delimited")
	parser.add_argument("--uncorrected_allele_file", 
							required=True, 
							help="allele file (not error corrected), tab-delimited")
	parser.add_argument("--corrected_50nt_file", 
							required=True, 
							help="50nt file (error corrected), tab-delimited")
	parser.add_argument("--uncorrected_50nt_file", 
							required=True, 
							help="50nt file (not error corrected), tab-delimited")
	parser.add_argument("--sample_name", 
							required=True, 
							help="sample or chip name, which is used as a prefix in the output filenames")
	parser.add_argument("--cutoff", 
							type=int, 
							default=0, 
							help="minimum read depth to take (default: 0)")
	parser.add_argument("--output_dir",
							required=True,
							help="folder to save output files (default: current working directory)")
	parser.add_argument("--debug")

	args = parser.parse_args()

	design_file_dict = read_file(args.oligo_target_file, hash=True)

	uncorrected_allele = read_file(args.uncorrected_allele_file)
	corrected_allele = read_file(args.uncorrected_allele_file)
	uncorrected_50nt = read_file(args.uncorrected_50nt_file)
	corrected_50nt = read_file(args.corrected_50nt_file)

	
	# start pre-processing/filtering step
	#TODO: alternative naming for "pass" and "fail". These don't accurately reflect the contents of these variables
	#TODO: actually implement the below functions, lol
	corrected_allele_pass, corrected_allele_fail     = alt_allele_check(design_file_dict, corrected_allele)
	uncorrected_allele_pass, uncorrected_allele_fail = alt_allele_check(design_file_dict, uncorrected_allele)
	
	collated_failed_primers = find_conflicting_primers(corrected_allele_fail, uncorrected_allele_fail)

	filtered_corrected_allele, filtered_uncorrected_allele, filtered_corrected_50nt, filtered_uncorrected_50nt, allele_corrected_missing_in_uncorrected, allele_uncorrected_missing_in_corrected  = filter(corrected_allele_pass, uncorrected_allele_pass, corrected_50nt, uncorrected_50nt, collated_failed_primers, cutoff)



	#TODO: handle output_dir creation/management
	if not Path(args.output_dir).is_dir():
		Path(args.output_dir).mkdir()
	#elif 
	
	


	
