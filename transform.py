from __future__ import division
import argparse
import ntpath
import os
import errno
import math
import sys

from scipy import stats
import numpy as np


# from https://stackoverflow.com/questions/8384737/extract-file-name-from-path-no-matter-what-the-os-path-format
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def read_from_file(filename, cutoff=None):
    """a wrapper for the input file reading execution block."""
    #print(filename,cutoff)
    obj = []
    with open(filename, "r") as f:
        for row in f:
            entry = row.rstrip().split('\t')
            if cutoff != None: # if applicable, only include rows that exceed the minimum DP threshold (default: 0)
                if int(entry[8]) >= cutoff: # filter on read depth (column 9 of the relevant inputs)
                    obj.append(entry)
            else:
                obj.append(entry)

    return obj


def write_to_file(filestem, filename, obj, target_dir):
    """a wrapper for the output file writing execution block."""

    with open(target_dir + '/' + filestem + '-' + filename, "w") as f:
        for row in obj:
            f.write('\t'.join(str(x) for x in row))
            if row[-1] != '\n':
                f.write('\n')


def simulate_mosaic_input(args, filestem):
    """This function serves to simulate the creation of the MOSAIC_INPUT spreadsheet.
    The primary difference is that everything is held in memory, compared to
    being manipulated within the excel document.
    """


    # load the designfile
    c_designfile_tab = read_from_file(args.oligo_target_file)

#    write_to_file("output-c-designfile_tab.txt", c_designfile_tab, args.output_dir)

    # load the corrected allele file (e.g. *-alleleErrC.txt.input.txt)
    allele_tab_corrected = read_from_file(args.allele_corrected, args.cutoff)

    num_corrected_allele_entries = len(allele_tab_corrected)
    #print(num_corrected_allele_entries)

    # load the uncorrected allele file (e.g. *-alleleNoC.txt.input.txt)
    allele_tab_uncorrected = read_from_file(args.allele_uncorrected, args.cutoff)

    # glue them together (create ALLELE tab)
    allele_tab = allele_tab_corrected + allele_tab_uncorrected

    num_total_allele_entries = len(allele_tab)
    #print(num_total_allele_entries)

    write_to_file(filestem, "output-allele_tab.txt", allele_tab, args.output_dir)

    # ALLELE-TAB tab creation
    """
    The ALLELE-TAB file ultimately takes the form of a tab-separated text file with the following columns:
    column# column_letter columm_name
    ------- ------------- -----------
       0          A        UniqueID
       1          B        Primer_ID
       2          C        MutationID
       3          D        Chip-MutationID
       4          E        Chrom
       5          F        Start
       6          G        End
       7          H        Ref
       8          I        Alts
       9          J        DP (original DP computed by samtools mpileup)
       10         K        INFO
       11         L        AD
       12         M        DP (modified with Ref/Alt alleles)
       13         N        AAF
       14         O        % pass
       15         P        Project ID
       16         Q        Version (demultiplexed vs UMI)
       17         R        Raw/Corrected?
       18         S        Downsampling
       19         T        ChipID
       20         U        Local Realignment?
       21         V        (blank) ???????
       22         W        MutationID
    """
    allele_tab_tab = []

    allele_tab_AJ = [] # used for construction of mutation_ID in allele_tab_tab (column W)
    allele_tab_AL = [] # used for construction of mutation_ID in allele_tab_tab (column W)
    allele_tab_AK = [] # needed for construction of sample_id in sensitivity_tab_tab (column M)
    allele_tab_AM = [] # needed for construction of sample_id in sensitivity_tab_tab (column M)
    for row in allele_tab:
        allele_tab_AJ.append(row[35])
        allele_tab_AL.append(row[37])
        allele_tab_AK.append(row[36])
        allele_tab_AM.append(row[38])

    cdf_mutation_ids = [] # primer code from design file, used for construction of primer_id in allele_tab_tab (column B)
    cdf_labels = []       # primer name from design file, used for construction of primer_id in allele_tab_tab (column B)
    for row in c_designfile_tab:
        cdf_mutation_ids.append(row[0])
        cdf_labels.append(row[1])

    for row in allele_tab:
        # prefill in cells from allele_tab ("less" dynamic cells)
        entry = ["","",row[30],row[29]+'-'+row[36],row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[23],row[24],row[25],row[26],row[22],"yes","",""]
        mutation_id = allele_tab_AL[allele_tab_AJ.index(entry[2])]
        primer_id = cdf_labels[cdf_mutation_ids.index(mutation_id)]
        unique_id = row[22] + '-' + primer_id

        entry[-1] = mutation_id
        entry[1] = primer_id
        entry[0] = unique_id

        allele_tab_tab.append(entry)

    write_to_file(filestem, "output-allele_tab_tab.txt", allele_tab_tab, args.output_dir)

    window_tab_corrected = read_from_file(args.corrected_50nt, args.cutoff)

    window_tab_uncorrected = read_from_file(args.uncorrected_50nt, args.cutoff)

    window_tab = window_tab_corrected + window_tab_uncorrected

    write_to_file(filestem, "output-window_tab.txt", window_tab, args.output_dir)

    # WINDOW-TAB tab creation
    """
    The WINDOW-TAB output file ultimately takes the form of a tab-separated text file with the following columns:
    column# column_letter column_name
    ------- ------------- -----------
       0          A        UniqueID
       1          B        Primer_ID
       2          C        MutationID
       3          D        Chip-MutationID
       4          E        Chrom
       5          F        Start
       6          G        End
       7          H        Ref
       8          I        Alts
       9          J        DP (original DP computed by samtools mpileup)
       10         K        INFO
       11         L        AD
       12         M        DP (modified with Ref/Alt alleles)
       13         N        AAF
       14         O        % pass
       15         P        Project ID
       16         Q        Version (demultiplexed vs UMI)
       17         R        Raw/Corrected?
       18         S        Downsampling
       19         T        ChipID
       20         U        Local Realignment?
       21         V        (blank) ???????
       22         W        MutationID
    """
    window_tab_tab = []

    window_tab_AJ = []
    window_tab_AL = []
    for row in window_tab:
        window_tab_AJ.append(row[35])
        window_tab_AL.append(row[37])

    for row in window_tab:
        entry = ["","",row[30],row[29]+'-'+row[36],row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[23],row[24],row[25],row[26],row[22],    "yes","",""]

        mutation_id = window_tab_AL[window_tab_AJ.index(entry[2])]
        primer_id = cdf_labels[cdf_mutation_ids.index(mutation_id)]
        unique_id = row[22] + '-' + primer_id

        entry[-1] = mutation_id
        entry[1] = primer_id
        entry[0] = unique_id

        window_tab_tab.append(entry)

    write_to_file(filestem, "output-window_tab_tab.txt", window_tab_tab, args.output_dir)

    # SENSITIVITY-TAB tab creation
    """
    The sensitivity-tab output file ultimately takes the form of a tab-separated text file with the following columns:
    column# column_letter column_name
    ------- ------------- -----------
       0          A        MutationID
       1          B        Chip-MutationID
       2          C        UniqueID
       3          D        PrimerID
       4          E        Primer#_in_Batch
       5          F        PrimerID
       6          G        Chromosome
       7          H        AlleleStart
       8          I        AlleleEnd
       9          J        Ref
       10         K        Alt
       11         L        Gene
       12         M        Sample ID
       13         N        Product Start
       14         O        Product end
       15         P        InsertStart
       16         Q        InsertEnd
       17         R        Forward
       18         S        Reverse
       19         T        Longer Barcode (15nt)
       20         U        Primer#
    """
    sensitivity_tab_tab = []

    allele_tab_tab_A = [] # used to construct primer#_in_batch (column E)
    allele_tab_tab_W = [] # used to construct primer#_in_batch (column E)

    # collect values
    for row in allele_tab_tab:
        allele_tab_tab_A.append(row[0])
        allele_tab_tab_W.append(row[22])

    # begin building the table row by row
    for i in range(num_corrected_allele_entries):
        primer_id = allele_tab_tab[i][1]
        entry = [allele_tab_tab[i][2][:-2],allele_tab_tab[i][3],allele_tab_tab[i][0],primer_id,"",primer_id,"","","","","","NA","","","","","","","","",primer_id[-1]]

        primer_no_in_batch = allele_tab_tab_W[allele_tab_tab_A.index(entry[2])]
        sample_id = allele_tab_AM[allele_tab_AK.index(entry[0])]

        entry[4] = primer_no_in_batch
        entry[12] = sample_id

        sensitivity_tab_tab.append(entry)

    write_to_file(filestem, "output-sensitivity_tab_tab.txt", sensitivity_tab_tab, args.output_dir)

    mosaic_input_output = [allele_tab, window_tab, allele_tab_tab, window_tab_tab, sensitivity_tab_tab, c_designfile_tab]

    return mosaic_input_output
    

def simulate_sensitivity_assay(mosaic_input, downsample):
    """This function simulates the creation of the SENSITIVITYASSAY spreadsheet.
    The primary difference is that everything is held in memory, compared to
    being manipulated within the excel document.
    """

    q20_allele_tab = mosaic_input[2]

    q20_allele_tab_A = []
    q20_allele_tab_D = []
    q20_allele_tab_L = []
    q20_allele_tab_M = []
    q20_allele_tab_N = []
    q20_allele_tab_O = []
    q20_allele_tab_R = []
    q20_allele_tab_S = []
    for row in q20_allele_tab:

        tmpCMID_components = row[3].split('-')
        tmpCMID_components[2] = tmpCMID_components[2][tmpCMID_components[2].find('N'):]
        row[3] = '-'.join(tmpCMID_components)

        q20_allele_tab_A.append(row[0])
        q20_allele_tab_D.append(row[3])
        q20_allele_tab_L.append(row[11])
        q20_allele_tab_M.append(row[12])
        q20_allele_tab_N.append(row[13])
        q20_allele_tab_O.append(row[14])
        q20_allele_tab_R.append(row[17])
        q20_allele_tab_S.append(row[18])
    
    q20_nts_tab = mosaic_input[3]

    q20_nts_tab_A = []
    q20_nts_tab_D = []
    q20_nts_tab_M = []
    q20_nts_tab_N = []
    q20_nts_tab_R = []
    q20_nts_tab_S = []
    for row in q20_nts_tab:
        tmpCMID_components = row[3].split('-')
        tmpCMID_components[2] = tmpCMID_components[2][tmpCMID_components[2].find('N'):]
        row[3] = '-'.join(tmpCMID_components)

        q20_nts_tab_A.append(row[0])
        q20_nts_tab_D.append(row[3])
        q20_nts_tab_M.append(row[12])
        q20_nts_tab_N.append(row[13])
        q20_nts_tab_R.append(row[17])
        q20_nts_tab_S.append(row[18])

    # SENSITIVITY-ANALYSIS tab creation
    """
    The sensitivity-analysis output file ultimately takes the form of a tab-separated text file with the following columns:
    column# column_letter column_name
    ------- ------------- -----------
       0          A        MutationID
       1          B        Chip-MutationID
       2          C        UniqueID
       3          D        PrimerID
       4          E        Primer#_in_Batch
       5          F        PrimerID
       6          L        Gene
       7          M        Sample ID
       8          U        Primer#
       9          X        UMI?
       10         Y        Pool#
       11         AA       Description of Pool
       12         AC       Input DNA (ng)
       13         AD       Method (UMI/NonUMI)
       14         AF       Dilution
       15         AG       Expected AAF (from original sequencing or dilution)
       16         AI       Raw: IT Read Depth
       17         AJ       Raw: Alt Allele Depth
       18         AK       Raw: AAF (pileup)
       19         AL       Raw: % NTs Passing QC
       20         AM       Raw: Background AAF (within 50nts)
       21         AN       Raw: Stdev Background AAF
       22         AO       Raw: Variance Background AAF
       23         AP       Raw: Stdev of average Background
       24         AQ       Raw: Confidence Interval of Background
       25         AR       IT Read Depth
       26         AS       Alt Allele Depth
       27         AT       AAF (pileup)
       28         AU       % NTs Passing QC
       29         AV       Background AAF (within 50nts)
       30         AW       Stdev Background AAF
       31         AX       Variance Background AAF
       32         AY       Stdev of Average Background
       33         AZ       Confidence Interval of Background
       34         BA       Average AAF
       35         BB       Confidence interval of Average AAF
       36         BC       P value (t-test) of allele frequencees to background error rate
       37         BD       Downsample
       38         BE       Delta AAF (Uncorrect vs correct)
       39         BF       Ratio AAF (Uncorrect vs correct)
       40         BG       Delta Error (Uncorrect vs correct)
       41         BH       Ratio Error (Uncorrect vs correct)
       42         BI       % Reduced Error
    """
    sensitivity_analysis_tab = []

    # begin building the table row by row
    for row in mosaic_input[4]:
        entry = row[:6] # A-F
        entry.extend(row[11:13]) # L, M
        entry.extend(row[20])    # U
        entry.extend(["","","","50ng","nonUMI","",""]) # X, Y, AA, AC, AD, AF, AG
        entry.extend([""]*21) # placeholders for AI-BC
        entry.append(downsample) #BD
        entry.extend([""]*5) # BE-BI
        chipmutationID_components = entry[1].split('-')
        chipmutationID_components[2] = chipmutationID_components[2][chipmutationID_components[2].find('N'):]
        entry[1] = '-'.join(chipmutationID_components)

        # match on PrimerID
        check1 = []
        for pos, row2 in enumerate(q20_allele_tab_A):
            if row2 == entry[2]:
                check1.append(pos)

        # match on UniqueID (single amplicon)
        check2 = []
        for pos, row2 in enumerate(q20_nts_tab_A):
            if row2 == entry[2]:
                check2.append(pos)

        # match on Chip-MutationID (all amplicons)
        check3 = []
        for pos, row2 in enumerate(q20_nts_tab_D):
            if row2 == entry[1]:
                check3.append(pos)

        for pos in check1:
            if q20_allele_tab_S[pos] == entry[37]:
                if q20_allele_tab_R[pos] == "uncorrected":
                    entry[16] = q20_allele_tab_M[pos]
                    entry[17] = q20_allele_tab_L[pos]
                    entry[18] = float(q20_allele_tab_N[pos])
                    entry[19] = q20_allele_tab_O[pos]
                if q20_allele_tab_R[pos] == "Corrected":
                    entry[25] = q20_allele_tab_M[pos]
                    entry[26] = q20_allele_tab_L[pos]
                    entry[27] = float(q20_allele_tab_N[pos])
                    entry[28] = q20_allele_tab_O[pos]

        vals_unc = []
        vals_c = []
        ns_unc = []
        ns_c = []
        for pos in check2:
            if q20_nts_tab_S[pos] == entry[37]:
                if q20_nts_tab_R[pos] == "uncorrected":
                    vals_unc.append(float(q20_nts_tab_N[pos]))
                    ns_unc.append(float(q20_nts_tab_N[pos]))
                                
                if q20_nts_tab_R[pos] == "Corrected":
                    vals_c.append(float(q20_nts_tab_N[pos]))
                    ns_c.append(float(q20_nts_tab_N[pos]))
                    
        vals_c_error = []
        vals_unc_error = []
        for pos in check3:
            if q20_nts_tab_S[pos] == entry[37]:
                if q20_nts_tab_R[pos] == "Corrected":
                    vals_c_error.append(float(q20_nts_tab_N[pos]))
                    # THESE VALUES GET PROCESSED AT THE VERY END IN ORDER TO AVOID HAVING TO REINDEX EVERYTHING IN BETWEEN
                    # THIS NEEDS TO BE FIXED LATER
                if q20_nts_tab_R[pos] == "uncorrected":
                    vals_unc_error.append(float(q20_nts_tab_N[pos]))

        if len(vals_unc) > 0:
            entry[20] = sum(vals_unc) / len(vals_unc)
            #entry[21] = math.sqrt(sum(pow(x-entry[20],2) for x in vals_unc) / (len(vals_unc)-1))
            entry[21] = np.std(vals_unc, ddof=1)

        if len(vals_c) > 0:
            entry[29] = sum(vals_c) / len(vals_c)
            entry[30] = math.sqrt(sum(pow(x-entry[29],2) for x in vals_c) / (len(vals_c)-1))

        if len(ns_unc) > 0:
            mean_ns_unc = sum(ns_unc) / len(ns_unc)
            #entry[22] = sum(pow(x-mean_ns_unc,2) for x in ns_unc) / (len(ns_unc) - 1)
            entry[22] = np.var(ns_unc, ddof=1)


        if len(ns_c) > 0:
            mean_ns_c = sum(ns_c) / len(ns_c)
            #entry[31] = sum(pow(x-mean_ns_c,2) for x in ns_c) / (len(ns_c) - 1)
            entry[31] = np.var(ns_c, ddof=1)

        # this is first appended to the end of the row due to indexing disgustingness, fixed below

        # average background error
        if len(vals_c_error) > 0:
            entry[32] = np.std(vals_c_error, ddof=1)
            mean_vals_error = np.mean(vals_c_error)
            #ci_raw = stats.norm.interval(0.95, loc=mean_vals_error, scale=entry[32]/math.sqrt(len(vals_c_error)))
            ci_raw = stats.t.interval(alpha=0.95, loc=mean_vals_error, df=len(vals_c_error)-1, scale=stats.sem(vals_c_error))
            if not np.isnan(ci_raw[1]):
                entry[33] = ci_raw[1] - mean_vals_error
            else:
                entry[33] = ""
            entry.append(mean_vals_error)
        else:
            entry[32] = ""
            entry[33] = ""
            #entry.append("")

        # raw average background error
        mean_vals_unc_error = ""
        if len(vals_unc_error) > 0:
            entry[23] = np.std(vals_unc_error, ddof=1)
            mean_vals_unc_error = np.mean(vals_unc_error)
            #ci_raw = stats.norm.interval(0.95, loc=mean_vals_unc_error, scale=entry[23]/math.sqrt(len(vals_unc_error)))
            ci_raw = stats.t.interval(alpha=0.95, loc=mean_vals_unc_error, df=len(vals_unc_error), scale=stats.sem(vals_unc_error))
            if not np.isnan(ci_raw[1]):
                entry[24] = ci_raw[1] - mean_vals_unc_error
            else:
                entry[24] = ""
            entry.append(mean_vals_unc_error)
        else:
            entry[23] = ""
            entry[23] = ""
            #entry.append("")

        sensitivity_analysis_tab.append(entry)
    
    sensitivity_B = []
    sensitivity_AK = []
    sensitivity_AM = []
    sensitivity_AO = []
    sensitivity_AT = []
    sensitivity_AV = []
    sensitivity_AX = []
    sensitivity_BD = []

    for row in sensitivity_analysis_tab:
        sensitivity_B.append(row[1])
        sensitivity_AK.append(row[18])
        sensitivity_AM.append(row[20])
        sensitivity_AO.append(row[22])
        sensitivity_AT.append(row[27])
        sensitivity_AV.append(row[29])
        sensitivity_AX.append(row[31])
        sensitivity_BD.append(row[37])

    for row in sensitivity_analysis_tab:
        #print("row:")
        #print(row)
        #print("row[1]:")
        #print(row[1])
        #print("row[37]: ")
        #print(row[37])
        #print("sensitivity_B:")
        #print(sensitivity_B)
        #print("sensitivity_BD:")
        #print(sensitivity_BD)
        #print("sensitivity_AO:")
        #print(sensitivity_AO)

        #AOs = [float(sensitivity_AO[pos]) for pos in range(len(sensitivity_AO)) if all([sensitivity_BD[pos] == row[37], sensitivity_B[pos] == row[1]])]
        #mean_AO = sum(AOs) / len(AOs)
        #row[23] = math.sqrt(mean_AO)
       
        try:
            AKs = [float(sensitivity_AK[pos]) for pos in range(len(sensitivity_AK)) if all([sensitivity_BD[pos] == row[37], sensitivity_B[pos] == row[1]])]
        except ValueError:
            AKs = []
 
 #       AMs = [float(sensitivity_AM[pos]) for pos in range(len(sensitivity_AM)) if all([sensitivity_BD[pos] == row[37], sensitivity_B[pos] == row[1]])]
        

#        if len(AMs) > 1 and len(set(AMs)) > 1:
#            mean_AMs = np.mean(AMs)
#            #print(mean_AMs)
#            ci_raw = stats.norm.interval(0.95, loc=mean_AMs, scale=np.std(AMs)/math.sqrt(len(AMs)))
#            if not np.isnan(ci_raw[1]):
#                row[24] = ci_raw[1] - mean_AMs
#            else:
#                row[24] = ""
#        else:
#            row[24] = ""

#        AVs = [float(sensitivity_AV[pos]) for pos in range(len(sensitivity_AV)) if all([sensitivity_BD[pos] == row[37], sensitivity_B[pos] == row[1]])]

#        AXs = [float(sensitivity_AX[pos]) for pos in range(len(sensitivity_AX)) if all([sensitivity_BD[pos] == row[37], sensitivity_B[pos] == row[1]])]
#        mean_AX = sum(AXs) / len(AXs)
#        row[32] = math.sqrt(mean_AX)

#        if len(AVs) > 1 and len(set(AVs)) > 1:
#            mean_AVs = np.mean(AVs)
#            ci = stats.norm.interval(0.95,loc=mean_AVs,scale=np.std(AVs)/math.sqrt(len(AVs)))
#            if not np.isnan(ci[1]):
#                row[33] = ci[1] - mean_AVs
#            else:
#                row[33] = ""
#        else:
#            row[33] = ""

        # because it might be the case that there isn't an associated raw AAF to pull due to thresholding, we need to check for any empty values.
        # if we encounter an empty value, throw out the calculation and use a blank instead
        try:
            raw_aafs = [float(q20_allele_tab_N[pos]) for pos in range(len(q20_allele_tab_N)) if all([q20_allele_tab_D[pos] == row[1], q20_allele_tab_R[pos] == "uncorrected", q20_allele_tab_S[pos] == row[37]])]
            row.append(sum(raw_aafs) / len(raw_aafs)) # raw average AAF
        except ValueError:
            row.append("")

        aafs = [float(q20_allele_tab_N[pos]) for pos in range(len(q20_allele_tab_N)) if all([q20_allele_tab_D[pos] == row[1], q20_allele_tab_R[pos] == "Corrected", q20_allele_tab_S[pos] == row[37]])]
        row[34] = sum(aafs) / len(aafs)

        ATs = [float(sensitivity_AT[pos]) for pos in range(len(sensitivity_AT)) if all([sensitivity_BD[pos] == row[37], sensitivity_B[pos] == row[1]])]
        
        if len(ATs) > 1 and len(set(ATs)) > 1:
            mean_ATs = np.mean(ATs)
            #ci_average = stats.norm.interval(0.95, loc=mean_ATs, scale=np.std(ATs)/math.sqrt(len(ATs)))
            ci_average = stats.t.interval(alpha=0.95, loc=mean_ATs, df=len(ATs)-1, scale=stats.sem(ATs))
            if not np.isnan(ci_average[1]):
                row[35] = ci_average[1] - mean_ATs
            else:
                row[35] = ""
        else:
            row[35] = ""

        #print(ATs)
        #print(var_ATs)
        #print(row[34])
        #print(row)
        
        if row[18]:
            if row[18] > 0:
                row[38] = row[18] - row[27]
                if row[27] > 0:
                    row[39] = row[18] / row[27]
            # this currently does not work as intended at the moment
            # due to it being nontrivial to test equality of floats
            # but if this elif evaluates, that simply means that there are no meaningful
            # values in either AAF (pileup) column.
            elif np.isclose([row[18]], [0], atol=1e-8, equal_nan=False):
                row[38] = 0

        if row[20] and row[20] > 0:
            row[40] = row[20] - row[29]
            if row[29] > 0:
                row[41] = row[20] / row[29]
            
            row[42] = row[40] / row[20] * 100

        if len(AKs) > 1 and len(set(AKs)) > 1:
            mean_AKs = np.mean(AKs)
            #print(meanAKs)
            #ci_raw_average_aaf = stats.norm.interval(0.95, loc=mean_AKs, scale=np.std(AKs, ddof=1)/math.sqrt(len(AKs)))
            ci_raw_average_aaf = stats.t.interval(alpha=0.95, loc=mean_AKs, df=len(AKs)-1, scale=stats.sem(AKs))
            if not np.isnan(ci_raw_average_aaf[1]):
                row.append(ci_raw_average_aaf[1] - mean_AKs)
            else:
                row.append("")
        else:
            #print(-1)
            row.append("")

        # re-inserts the previously appended row values at the desired position
        row.insert(25, row.pop())   # raw confidence interval on average AAF
        row.insert(25, row.pop())   # raw average AAF
        row.insert(21, row.pop())   # raw: average background AAF
        row.insert(33, row.pop())   # average background AAF
        row.insert(28, "")            # raw: p-value (t-test)

        row.pop(-5)   # removes Delta AAF (we no longer need it, but is currently impractical to untangle it upstream at the moment)

    sensitivity_analysis_tab.insert(0, [
                                        "MutationID",
                                        "Chip-MutationID",
                                        "UniqueID",
                                        "PrimerID",
                                        "Primer#_in_Batch",
                                        "PrimerID",
                                        "Gene",
                                        "Sample ID",
                                        "Primer#",
                                        "UMI?",
                                        "Pool#",
                                        "Description of Pool",
                                        "Input DNA (ng)",
                                        "Method (UMI/NonUMI)",
                                        "Dilution",
                                        "Expected AAF (from original sequencing or dilution)",
                                        "Raw: IT Read Depth",
                                        "Raw: Alt Allele Depth",
                                        "Raw: AAF (pileup)",
                                        "Raw: % NTs Passing QC",
                                        "Raw: Background AAF (within 50nts)",
                                        "Raw: Average Background AAF",
                                        "Raw: Stdev Background AAF",
                                        "Raw: Variance Background AAF",
                                        "Raw: Stdev of Average Background AAF",
                                        "Raw: Confidence interval of Background AAF",
                                        "Raw: Average AAF",
                                        "Raw: Confidence Interval of Average AAF",
                                        "Raw: p-value (t-test)",
                                        "IT Read Depth",
                                        "Alt Allele Depth",
                                        "AAF (pileup)",
                                        "% NTs Passing QC",
                                        "Background AAF (within 50nts)",
                                        "Average Background AAF",
                                        "Stdev Background AAF",
                                        "Variance Background AAF",
                                        "Stdev of Average Background AAF",
                                        "Confidence interval of Background AAF",
                                        "Average AAF",
                                        "Confidence interval of Average AAF",
                                        "P value (t-test) -- of allele frequencies to background error rate.",
                                        "Downsample",
#                                        "Delta AAF (Uncorrect vs correct)",
                                        "Ratio AAF (Uncorrect vs correct)",
                                        "Delta Error (Uncorrect vs correct)",
                                        "Ratio Error (Uncorrect vs correct)",
                                        "% Reduced Error"
                                       ])
        
    return sensitivity_analysis_tab


if __name__ == "__main__":
    # create script arguments as a way to standardize input file passing
    parser = argparse.ArgumentParser()
    parser.add_argument("--allele_corrected", help="allele file (error corrected), tab-delimited")
    parser.add_argument("--allele_uncorrected", help="allele file (not error corrected), tab-delimited")
    parser.add_argument("--corrected_50nt", help="50nt file (error corrected), tab-delimited")
    parser.add_argument("--uncorrected_50nt", help="50nt file (not error corrected), tab-delimited")
    parser.add_argument("--oligo_target_file", help="file containing list of oligos targeting which sites, tab-delimited")
    parser.add_argument("--cutoff", help="minimum read depth to take (default: 0)", type=int, default=0)
    parser.add_argument("--output_dir", help="folder to save outputs")
    parser.add_argument("--remove_duplicates", action="store_true", default=False, help="eliminate duplicate output rows before writing to file")

    args = parser.parse_args()

    try:
        os.makedirs(args.output_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            print("output directory %s already found, skipping creation" % args.output_dir)

    basename = path_leaf(args.allele_corrected)
    filestem = '-'.join(basename.split('-')[:3])

    # generate first spreadsheet-representative object
    mosaic_input = simulate_mosaic_input(args, filestem)

	# check for consistent downsampling values
    downsampling_values = []
    for thing in mosaic_input[:-2]:
        try:
            downsampling_values.extend([row[26] for row in thing])
        except IndexError:
            downsampling_values.extend([row[18] for row in thing])
    downsampling_values = set(downsampling_values)
    if len(downsampling_values) > 1:
        sys.exit(1, "Downsampling values not consistent; aborting. Please ensure all Downsampling columns in input files are consistent.")

    # generate second spreadsheet-representative object 
    sensitivity_assay = simulate_sensitivity_assay(mosaic_input, list(downsampling_values)[0])

    filtered_assay = []
    if args.remove_duplicates:
        for row in sensitivity_assay:
            if row not in filtered_assay:
                filtered_assay.append(row)
    else:
        filtered_assay = sensitivity_assay

    write_to_file(filestem, "output-sensitivity_assay.txt", filtered_assay, args.output_dir)

