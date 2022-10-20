1. Input
    - toggles (if specified, turns on that option): 
        - `-help|h`
            - help output can also be toggled if you forget some input files or parameters
        - `-debug`
            - will not print extra debug output unless specified
        - `-remove_ref_lt_alt`
            - data will not be removed if reference counts are fewer than the alt counts unless you switch this on
    - input files (requires filename with each option): 
        - `-oligo_target_file` - Design file
        - `-allele_corrected`
        - `-allele_uncorrected`
        - `-corrected_50nt`
        - `-uncorrected_50nt`
    - misc input parameters (requires specifying corresponding input for each of the below options):
        - `-output_dir`
        - `-sample_name` - used in constructing output file names
        - `-cutoff` 
            - will use `0` unless specified
        - `-qc_pass_filter`
            - default is `.7` unless specified
2. Output
    - output goes to directory specified with `output_dir`. A new directory must be used for each run of the tool.
        - Subdirectories:
            - `fail` - data that did not adhere to all specified conditions
            - `pass` - data that adhered to all specified conditions
    - allele check output
        - `fail` - files that did not pass the alt allele check function
            - `samplename-alleleErrC.failAltAlleleCheck.txt`
            - `samplename-alleleNoC.failAltAlleleCheck.txt`
        - `pass` - files that passed the alt allele check function
            - `samplename-alleleErrC.passAltAlleleCheck.txt`
            - `samplename-alleleNoC.failAltAlleleCheck.txt`
    - filter output
        - `fail` 
            - data below the DP cutoff
                - `samplename-alleleErrC.passAltAlleleCheck.belowCutoff.txt`
                - `samplename-alleleNoC.passAltAlleleCheck.belowCutoff.txt`
                - `samplename-50ntErrC.belowCutoff.txt`
                - `samplename-50ntNoC.belowCutoff.txt`
            - data missing in one of the allele files
                - `samplename-alleleErrC.passAltAlleleCheck.missingPrimersInAlleleNoC.txt`
                - `samplename-alleleNoC.passAltAlleleCheck.missingPrimersInAlleleErrC.txt`
        - `pass` 
            - data above DP cutoff
                - `samplename-alleleErrC.passAltAlleleCheck.aboveCutoff.txt`
                - `samplename-alleleNoC.passAltAlleleCheck.aboveCutoff.txt`
                - `samplename-50ntErrC.aboveCutoff.txt`
                - `samplename-50ntNoC.aboveCutoff.txt`
    - ref vs alt output
        - `pass`
            - files with ref > alt
                - `samplename-alleleErrC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refGTalt.txt`
                - `samplename-alleleNoC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refGTalt.txt`
                - `samplename-50ntErrC.aboveCutoff.aboveQCfilter.refGTalt.txt`
                - `samplename-50ntNoC.aboveCutoff.aboveQCfilter.refGTalt.txt`
        - `fail`
            - files with ref < alt
                - `samplename-alleleErrC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refLTalt.txt`
                - `samplename-alleleNoC.passAltAlleleCheck.aboveCutoff.aboveQCfilter.refLTalt.txt`
                - `samplename-50ntErrC.aboveCutoff.aboveQCfilter.refLTalt.txt`
                - `samplename-50ntNoC.aboveCutoff.aboveQCfilter.refLTalt.txt`

3. Functions
    - `read_file()`
        - **Input**: filename & `1` or `0` to designate if this is the design file
        - **Output**: hash for design file (so we can search it easily), or array for any other file
        - **Method**:
            - for design file:
                - split on tabs
                - columns of interest (indexing by 0):
                    - 0: primer
                    - 2: chromosome
                    - 3: loci
                    - 5: ref
                    - 6: alt
                - key for hash is `chr_loci_primer`
                - before adding the current line's information to the hash, we check if this key (`chr_loci_primer`) already exists.
                    - if it does exist:
                        - give warning if reference alleles for same chr, loci and primer do not match
                        - give warning if primary alt alleles for same chr, loci, primer do not match
                        - re-add to hash (useless step if ref and alt alleles match, shouldn't reach this if a warning is caught)
                    - if it doesn't exist: add to hash
            - for all other files: remove newlines, put each line as an element in an array
    - `write_file()`
       - **Input**: variable to write, directory, sample name for filename, suffix for filename, `1` or `0` to specify if we have a hash (1) or an array 
       - **Output**: File written to `directory/samplename-suffix.txt`
       - **Method**: For arrays, we write the full contents to file. For hashes, we write just the values.
    - `remove_data_if_alt_allele_check_failed()`
        - **Input**: design file hash, allele corrected array, allele uncorrected array, corrected 50nt array, corrected 50nt array, debug flag
        - **Output**: `pass` and `fail` arrays for both corrected and uncorrected
        - **Method**: Sends allele corrected and allele uncorrected data to `alt_allele_check()` (individually); generates a list of primers that are not consistent with the design file (`find_primers_that_conflict_with_design_file`); if any such primers, removes this data using `remove_lines_by_specified_data()`; returns all "passing" and "failing" data from the aforementioned steps.
    - `alt_allele_check()`
        - **Input**: design file hash, allele uncorrected OR allele corrected array, debug flag
        - **Output**: "pass" and "fail" arrays for allele corrected/uncorrected
        - **Method**: 
            - split lines on tabs, save fields of interest
            - check we have chr + loci + primer in design file
                - error if chr + loci + primer combo if not found (never seen this happen), fail
            - verify ref is the same as ref in design file
                - error if ref is not the same (never seen this happen), fail
            - check alt allele field:
                - pass if REF <*>
                - pass if expected alt allele in primary position
                - pass if expected alt allele in secondary/tertiary/etc. position
                    - modify columns 8-11
                - pass if alt allele that is not what we expect
                    - modify columns 8-11
                - error if condition does not match one of the above (never seen this happen), fail 
        - **Modifications of cols 8-11 (counting from 1)**
            - expected alt allele in non-primary position:
                - col 8 (ALT): correct alt count from AD
                    - pull out based upon position of the expected primary alt allele within original col 8, plus 1 (the 1 is to account for reference being listed first)
                - col 9 (TOTAL DP): REF + ALT
                    - REF = position 0 when AD (col 6) is split on commas, and "AD=" is removed
                    - ALT = recalculated col 8
                - col 10 (AAF): new col 8 / new col 9
                    - col 8 = correct alt counts
                    - col 9 = total DP
                - col 11 (%nt QC): new col 9 / col 6
                    - col 9 = total DP
                    - col 6 = total read count
            - expected alt allele not found:
                - col 8 (ALT): 0
                - col 9 (TOTAL DP): REF from AD
                    - REF from AD = position 0 when AD (col 6) is split on commas, and "AD=" is removed
                - col 10 (AAF): 0
                - col 11 (%nt QC): new col 9 / col 6
                    - col 9 = total DP
                    - col 6 = total read count
    - `find_primers_that_conflict_with_design_file()`
        - **Input**: array with allele corrected "fail" lines, array with allele error corrected "fail" lines, debug flag
        - **Output**: array with the combined list of "fail" primers to remove 
        - **Method**: 
            - called by `remove_data_if_alt_allele_check_failed`
            - arrays with primers that failed for allele corrected or allele uncorrected data are combined into one list.
            - primer fields are extracted (lines split on tabs, pull out field 38 when counting starting at 1), put into a primers to remove array.
    - `remove_lines_by_specified_data()`
        - **Input**: corrected OR uncorrected 50nt array, array containing list of "fail" primers, field number to check for a match, debug flag
        - **Output**: corrected 50nt "pass" array, uncorrected 50nt "pass" array
        - **Method**: 
            - called by `remove_data_if_alt_allele_check_failed` if there are any "failed" primers from the alt allele checking stage.
            - `remove_lines_by_specified_data()` is designed for eliminating lines from 50nt files for primers that failed the alt allele check.
            - this is called twice - once for corrected 50nt, and once for uncorrected 50nt.
            - the input parameter of "field number to check for a match" is specified as 37, which is the primer field when counting starting from 0.
            - the list of primers to remove is switched to a hash for easy searching. 
                - if data corresponding to such a primer is found in the input 50nt data, then it is discarded.
    - `filter()`
        - **Input**: "pass" arrays for allele corrected & uncorrected from alt_allele_check, 50 nt corrected & uncorrected "pass" arrays, cutoff, qc pass filter, remove if ref less than alt flag, debug flag
        - **Output**: below and above read depth cutoff for allele un/corrected and 50nt un/corrected; data in allele corrected but missing in allele uncorrected and vice versa; below and above qc filter for allele un/corrected and 50 nt un/corrected; allele un/corrected above qc pass with ref less than and greater than alt counts
        - **Method**: 
            - go through allele corrected pass. If above cutoff for DPREADS (col 9), save in "allele corrected above cutoff" hash (referenced below when looking through allele uncorrected)
                - else -> discard in allele corrected below cutoff" hash
            - go through allele uncorrected pass. 
                - both allele corrected + allele uncorrected must be >= cutoff to be retained.
                - If either allele corrected or allele uncorrected <= cutoff, put in relevant "below cutoff" hash & add to remove_from_50nt_files hash.
            - for all primer IDs in remove_from_50nt_files hash, remove relevant matching primer ID lines from 50 nt files
            - "missing" primer IDs that are in just one of allele corrected or allele uncorrected are identified & saved to variables
            - for qc pass filter:
                - go through each of these individually: (1) allele corrected pass, above DPREADS cutoff, (2) allele uncorrected pass, above DPREADS cutoff, then:
                    - segment each into "above qc filter" and "below qc filter hashes"
                    - remove discarded primers from the 50nt files that correspond to the appropriate allele collected file.
                        - also segment each into "above qc filter" and "below qc filter hashes"

    - `remove_data_if_ref_lt_alt()`
        - skipped if `-remove_ref_lt_alt` not toggled
        - **Input**: un/corrected 50nt above cutoff above qc filter, allele un/corrected above cutoff above qc filter, debug flag
        - **Output**: 
        - **Method**:
            - for each of the input allele data structures, segment based upon comparison of ref and alt counts.
                - if ref > alt counts, retain in "ref greater than alt" bucket
                - if ref < alt counts, discard into "ref less than alt" bucket
                - data discarded in allele un/corrected will be discarded in corresponding 50nt file using the appropriate primer.