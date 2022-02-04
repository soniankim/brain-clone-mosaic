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

	




if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()

	parser.add_argument("--oligo_target_file", 
							required=True, 
							help="file containing list of oligos targeting which sites, tab-delimited")
	parser.add_argument("--allele_corrected_file", 
							required=True, 
							help="allele file (error corrected), tab-delimited")
	parser.add_argument("--allele_uncorrected_file", 
							required=True, 
							help="allele file (not error corrected), tab-delimited")
	parser.add_argument("--50nt_corrected_file", 
							required=True, 
							help="50nt file (error corrected), tab-delimited")
	parser.add_argument("--50nt_uncorrected_file", 
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

	uncorrected_allele = read_file(args.allele_uncorrected_file)
	corrected_allele = read_file(args.allele_uncorrected_file)
	uncorrected_50nt = read_file(args.50nt_uncorrected_file)
	corrected_50nt = read_file(args.50nt_corrected_file)

	
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
	elif 
	
	


	
