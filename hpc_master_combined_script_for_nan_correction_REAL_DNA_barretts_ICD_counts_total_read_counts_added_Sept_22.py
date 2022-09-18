#!/data/home/hfx381/presto/bin/python

import pandas as pd
import numpy as np
import random
from itertools import combinations
from itertools import product
import os
import fnmatch
import glob
import fileinput
from datetime import datetime

print("Packages imported, commencing processing now", "", sep="\n")

starttime = datetime.now()

## EDIT -> CHOOSE DIRECTORY OF .txt files to process:
txt_file_directory = "/data/BCI-ESCS/richard/ngs_data/output_meth_reports/methseq_barretts/"

## EDIT -> and a NEW SUBDIRECTORY to save into:
databatch = "barretts_with_intercrypt/"

# EDIT -> choose which axis (0 = CpG Locus {column}; 1 = Read {row}) to determine the median value to replace NaNs
# If template is real world DNA then assign CpG locus medians (axis=0)
# If template target is Control DNA then assign read medians (axis=1)
cpglocus_or_read = 1 # keep as one - code re-written to transpose dataframe so rows are analysed

## EDIT -> The fraction of missing values within a cpg locus (column) or read (row) that results in filtering a particular cpg locus/read:
nan_fraction_cpg_locus_exclusion_threshold = 0.2 # i.e. Methylation status successfully called in 80% of the cpg locus/read

## EDIT -> Decide on threshold values to rename a NaN within a certin CpG Locus, i.e. probability based or random 0 / 1 chance
mean_1 = 0.2 # this is the mean value of the methylation calls per CpG Locus i.e 0.2 = 20% calls are methylated; 2x SD = 0.8
mean_2 = 0.8 # 0.8 = 80% of reported calls are methylated.

##EDIT -> change the value of the CpG loci number in the gene target variables depending which dataset is being analysed
## i.e. whether the primers have been trimmed by seqk or not in the data pre-processing (DEFAULT = 0)
cpg_loci_change = 0

master_df_output_directory = "/data/BCI-ESCS/richard/dataframes/"

try:
    os.mkdir(master_df_output_directory + databatch)
except FileExistsError:
    print("output folder already exists")

output_subdirectory = master_df_output_directory + databatch

#Move to the PATH where the txt files are:
os.chdir(txt_file_directory)

# Select the Files for Analysis within the Folder listed above -> os.chdir() -> Assign to a common variable
selected_txt_files = os.listdir()

# create a variable of the file names to add to final dataframes / dictionaries later
file_names = []
for i in range(len(selected_txt_files)):
    file_names.append(selected_txt_files[i][:-4])


open_files = []
for file in selected_txt_files:
    with open(file) as data:
        open_files.append(data.readlines())

print(file_names)

print("", "Files loaded into memory", "Proceeding to file clean-up", "", sep="\n")

#Multiplex Pool Number 5:
targets_in_dataset = ["ANKRD2", "CAMK2B", "CSRP3", "LOC-L", "1-MYOD1", "2-MYOD1", "1-NKX2-5", "2-NKX2-5", "NPPB", "PXDNL", "SBK2", "1-SBK3", "SCN5A", "1-TNNI3", "2-TNNI3"]

# split tab delimited text file into separate lists for each row - input can be any number of open files,
# each file is saved as its own list within the full list - all_files_datasplit
single_file_datasplit = []
def split_data(file):
    for data in file:
        single_file_datasplit.append(data.split("\t"))


all_files_datasplit = []
for i in range(len(open_files)):
    if single_file_datasplit == []:
        split_data(open_files[i])
        all_files_datasplit.append(single_file_datasplit)
    else:
        single_file_datasplit = []
        split_data(open_files[i])
        all_files_datasplit.append(single_file_datasplit)
    single_file_datasplit = []

# Test Code:
#filetest = open_files[2]
#datatest = all_files_datasplit[2]
#print(filetest[0:5])
#print(datatest[0:5])
#print(len(all_files_datasplit))


# Remove "Bismark methylation extractor version..." header from each dataset
no_headers = []
def header_removal(file):
    no_headers.append(file[1:])

for i in range(len(all_files_datasplit)):
    header_removal(all_files_datasplit[i])

# Test code:
#headtest = no_headers[2]
#print(headtest[0:5])
#print(len(headtest))


# Further clean-up - Remove newlines "\n"
newlines = []
def newline_removal(file):
    for data in file:
        for element in data:
            if "\n" in element:
                newlines.append(element.replace("\n", ""))
            else:
                newlines.append(element)


no_newlines = []
for i in range(len(no_headers)):
    if newlines == []:
        newline_removal(no_headers[i])
        no_newlines.append(newlines)
    else:
        newlines = []
        newline_removal(no_headers[i])
        no_newlines.append(newlines)
    newlines = []

# Test code:
#newlinestest = no_newlines[2]
#print(newlinestest[0:10])
#print(len(newlinestest))
#print(len(no_newlines))


# clean up gene data - remove all coordinates of genes

coordinates = []
def no_coord(file):
    for element in file:
        if "ANKRD2_10:97578093-97578518" in element:
            coordinates.append(element.replace("ANKRD2_10:97578093-97578518", "ANKRD2"))
        elif "CAMK2B_7:44324809-44325372" in element:
            coordinates.append(element.replace("CAMK2B_7:44324809-44325372", "CAMK2B"))
        elif "CSRP3_11:19210424-19210874" in element:
            coordinates.append(element.replace("CSRP3_11:19210424-19210874", "CSRP3"))
        elif "LOC_X:153485781-153486312" in element:
            coordinates.append(element.replace("LOC_X:153485781-153486312", "LOC_L"))
        elif "MYOD1_11:17719658-17720173" in element:
            coordinates.append(element.replace("MYOD1_11:17719658-17720173", "2_MYOD1"))
        elif "MYOD1_11:17720208-17720574" in element:
            coordinates.append(element.replace("MYOD1_11:17720208-17720574", "1_MYOD1"))
        elif "NKX2-5_5:173232483-173232952" in element:
            coordinates.append(element.replace("NKX2-5_5:173232483-173232952", "2_NKX2-5"))
        elif "NKX2-5_5:173234608-173235170" in element:
            coordinates.append(element.replace("NKX2-5_5:173234608-173235170", "1_NKX2-5"))
        elif "NPPB_1:11858993-11859542" in element:
            coordinates.append(element.replace("NPPB_1:11858993-11859542", "NPPB"))
        elif "PXDNL_8:51809235-51809721" in element:
            coordinates.append(element.replace("PXDNL_8:51809235-51809721", "PXDNL"))
        elif "SBK2_19:55535946-55536434" in element:
            coordinates.append(element.replace("SBK2_19:55535946-55536434", "SBK2"))
        elif "SBK3_19:55549846-55550306" in element:
            coordinates.append(element.replace("SBK3_19:55549846-55550306", "1_SBK3"))
        elif "SCN5A_3:38649483-38650055" in element:
            coordinates.append(element.replace("SCN5A_3:38649483-38650055", "SCN5A"))
        elif "TNNI3_19:55156291-55156830" in element:
            coordinates.append(element.replace("TNNI3_19:55156291-55156830", "1_TNNI3"))
        elif "TNNI3_19:55156839-55157273" in element:
            coordinates.append(element.replace("TNNI3_19:55156839-55157273", "2_TNNI3"))
        else:
            coordinates.append(element)

no_coordinates = []
for i in range(len(no_newlines)):
    if coordinates == []:
        no_coord(no_newlines[i])
        no_coordinates.append(coordinates)
    else:
        coordinates = []
        no_coord(no_newlines[i])
        no_coordinates.append(coordinates)
    coordinates = []

# Test code:
#genetest = no_coordinates[2]
#print(genetest[:5])
#print(len(genetest))
#print(len(no_coordinates))


# clean up barcode_id and barcode_counts data - needs further splitting of the string

bc_sep = []
def barcode(file):
    for element in file:
        if "CONSCOUNT" in element:
            barcode_split=[]
            barcode_split.append(element.split("|"))
            barcode_list_access=barcode_split[0]
            bc_sep.append(barcode_list_access[0])
            count_split=[]
            count_split.append(barcode_list_access[1].split("="))
            count_split_access=count_split[0]
            bc_sep.append(int(count_split_access[1]))
        else:
            bc_sep.append(element)

bc_tidied = []
for i in range(len(no_coordinates)):
    if bc_sep == []:
        barcode(no_coordinates[i])
        bc_tidied.append(bc_sep)
    else:
        bc_sep = []
        barcode(no_coordinates[i])
        bc_tidied.append(bc_sep)
    bc_sep = []

# Test code:
#bc_test = bc_tidied[2]
#print(bc_test[0:10])
#print(len(colontest))
#print(len(bc_tidied))

#print(len(no_coordinates))

print("", "Files cleaned-up", "Proceeding to dataframe construction", "", sep="\n")

# Data is now cleaned-up -> need to sort into dataframes now
clean_data = bc_tidied

# assign data values to new column headers, convert CpG Loci to integers (less 2 [added by Bismark]) and Meth/UnMeth Calls to Base 2

grouped_data = []
def group_data_columns(file):
    barcode_id = []
    barcode_count = []
    symbol_meth = []
    gene_target = []
    cpg_locus = []
    methZ_unmethz = []
    for i in list(range(0, len(file), 6)):
        barcode_id.append(file[i])
    for i in list(range(1, len(file), 6)):
        barcode_count.append(file[i])
    for i in list(range(2, len(file), 6)):
        symbol_meth.append(file[i])
    for i in list(range(3, len(file), 6)):
        gene_target.append(file[i])
    for i in list(range(4, len(file), 6)):
        cpg_locus.append(file[i])
    for i in list(range(5, len(file), 6)):
        methZ_unmethz.append(file[i])
    grouped_data.append(barcode_id)
    grouped_data.append(barcode_count)
    grouped_data.append(symbol_meth)
    grouped_data.append(gene_target)
    cpg_locus_int = []
    for num in cpg_locus:
        cpg_locus_int.append(int(num))
    cpg_locus_int_less2 = [num - 2 for num in cpg_locus_int]
    grouped_data.append(cpg_locus_int_less2)
    meth1_unmeth0 = []
    for data in methZ_unmethz:
        if data == "z":
            meth1_unmeth0.append(0)
        elif data == "Z":
            meth1_unmeth0.append(1)
        else:
            meth1_unmeth0.append(data)
    grouped_data.append(meth1_unmeth0)

reorg_data = []
for i in range(len(clean_data)):
    if grouped_data == []:
        group_data_columns(clean_data[i])
        reorg_data.append(grouped_data)
    else:
        grouped_data = []
        group_data_columns(clean_data[i])
        reorg_data.append(grouped_data)
    grouped_data = []

# Test code:
#reorgtest = reorg_data[2]
#si = reorgtest[0]
#gt = reorgtest[2]
#loci = reorgtest[3]
#m_um = reorgtest[4]
#print(si[0:10])
#print(gt)
#print(loci[0:10])
#print(m_um[0:10])
#print(reorgtest)
#print(len(reorgtest))
#print(len(si))
#print(len(m_um))
#print(reorg_data)
#print(len(reorg_data))


# Convert lists/files into individual Dataframes (saved within the list "df_list")

df_list = []
def dataframe_constructor(file):
    for data in file:
        df = pd.DataFrame(
        {"barcode_id": data[0], "barcode_count": data[1], "symbol_meth": data[2], "gene_target": data[3], "cpg_locus": data[4],
        "meth1_unmeth0": data[5]})
        df_list.append(df)

dataframe_constructor(reorg_data)

#print(df_list[0])
#print(len(df_list))


# Re-organise each df into dfs grouped by gene target

dfs_by_target = []
def group_by_target(list):
    for df in list:
        ankrd2 = df[df.gene_target.isin(["ANKRD2"])].reset_index(drop=True)
        camk2b = df[df.gene_target.isin(["CAMK2B"])].reset_index(drop=True)
        csrp3 = df[df.gene_target.isin(["CSRP3"])].reset_index(drop=True)
        loc_l = df[df.gene_target.isin(["LOC_L"])].reset_index(drop=True)
        myod1_1 = df[df.gene_target.isin(["1_MYOD1"])].reset_index(drop=True)
        myod1_2 = df[df.gene_target.isin(["2_MYOD1"])].reset_index(drop=True)
        nkx2_1 = df[df.gene_target.isin(["1_NKX2-5"])].reset_index(drop=True)
        nkx2_2 = df[df.gene_target.isin(["2_NKX2-5"])].reset_index(drop=True)
        nppb = df[df.gene_target.isin(["NPPB"])].reset_index(drop=True)
        pxdnl = df[df.gene_target.isin(["PXDNL"])].reset_index(drop=True)
        sbk2 = df[df.gene_target.isin(["SBK2"])].reset_index(drop=True)
        sbk3_1 = df[df.gene_target.isin(["1_SBK3"])].reset_index(drop=True)
        scn5a = df[df.gene_target.isin(["SCN5A"])].reset_index(drop=True)
        tnni3_1 = df[df.gene_target.isin(["1_TNNI3"])].reset_index(drop=True)
        tnni3_2 = df[df.gene_target.isin(["2_TNNI3"])].reset_index(drop=True)
        targets = [ankrd2, camk2b, csrp3, loc_l, myod1_1, myod1_2, nkx2_1, nkx2_2, nppb, pxdnl, sbk2, sbk3_1, scn5a, tnni3_1, tnni3_2]
        dfs_by_target.append(targets)

group_by_target(df_list)

#print(len(dfs_by_target))
#print(dfs_by_target)
#print(dfs_by_target[1])


# Select columns of interest before pivot
pre_pivot_dfs = []
def pre_pivot_prep(list):
    for data in list:
        pre_pivot = []
        for df in data:
            df_pre_pivot = df[["barcode_id", "barcode_count", "cpg_locus", "meth1_unmeth0"]]
            pre_pivot.append(df_pre_pivot)
        pre_pivot_dfs.append(pre_pivot)

pre_pivot_prep(dfs_by_target)

#print(len(pre_pivot_dfs))
#print(pre_pivot_dfs[0])
#print(pre_pivot_dfs[1])


# Pivot all dataframes into individual strand methylation tags

pivotted_dfs = []
def pivot(list):
    for data in list:
        pivot = []
        for df in data:
            df_pivot = df.pivot(columns="cpg_locus", index="barcode_id", values="meth1_unmeth0").reset_index()
            pivot.append(df_pivot)
        pivotted_dfs.append(pivot)

pivot(pre_pivot_dfs)

# Test code:
#print(len(pivotted_dfs))
#sample1 = pivotted_dfs[2]
#camk2b = sample1[10]
#print(camk2b)
#print(list(camk2b.columns.values))
#print(pivotted_dfs[0])
#print(pivotted_dfs[1])


# Need to create another dataframe with the barcode counts per barcode and then merge this with the pivotted dataframe above to retain this data

bc_count_dfs = []
def bc_count_pre_grouping(list):
    for data in list:
        bc_count_pre_grouping = []
        for df in data:
            bc_count_df_pre_grouping = df[["barcode_id", "barcode_count"]]
            bc_count_pre_grouping.append(bc_count_df_pre_grouping)
        bc_count_dfs.append(bc_count_pre_grouping)

bc_count_pre_grouping(dfs_by_target)


bc_count_grouped_dfs = []
def bc_count_singles(list):
    for data in list:
        bc_count_groups = []
        for df in data:
            tar_bc_count_groups = df.groupby(["barcode_id"]).mean().reset_index()
            bc_count_groups.append(tar_bc_count_groups)
        bc_count_grouped_dfs.append(bc_count_groups)

bc_count_singles(bc_count_dfs)


merged_dfs = []
def merge_dfs(df1, df2):
    for i in range(len(df1)):
        dataset1=df1[i]
        dataset2=df2[i]
        tar_dfs_merged = []
        for k in range(len(dataset1)):
            tar_df1=dataset1[k]
            tar_df2=dataset2[k]
            tar_df_merge = pd.merge(tar_df1, tar_df2, on="barcode_id")
            tar_dfs_merged.append(tar_df_merge)
        merged_dfs.append(tar_dfs_merged)

merge_dfs(pivotted_dfs, bc_count_grouped_dfs)

print("", "DataFrames made!", "Saving as CSV files...", "", sep="\n")

#print(merged_dfs[0])

## FINAL OUTPUT --> Now script to save dataframes into the output subdirectory as .CSV files

os.chdir(output_subdirectory)

bc_info_subdir = "barcode_information/"
try:
    os.mkdir(bc_info_subdir)
except FileExistsError:
    print("output folder for barcode information already exists")

os.chdir(output_subdirectory + bc_info_subdir)

bc_info_folder = "_barcode_info_dataframes"
bc_info_file = "_barcode_info_for_"
for i in range(len(file_names)):
    try:
        os.mkdir(file_names[i] + bc_info_folder)
    except FileExistsError:
        print("Barcode info output folder for " + file_names[i] + " already exists")

def make_barcode_csv_files (list):
    for i in range(len(list)):
        df_list=list[i]
        subsubdir = file_names[i] + bc_info_folder
        os.chdir(output_subdirectory + bc_info_subdir + subsubdir)
        for k in range(len(df_list)):
            df_list[k].to_csv(file_names[i] + bc_info_file + targets_in_dataset[k] + ".csv", index=False)

make_barcode_csv_files(merged_dfs)

# Make Dataframe of Raw consensus Reads successfully aligned for each sample and gene target
consensus_reads_list=[]
for i in range(len(merged_dfs)):
    sample=merged_dfs[i]
    all_targets=[]
    for k in range(len(sample)):
        target_length=len(sample[k])
        all_targets.append(target_length)
    consensus_reads_list.append(all_targets)

consensus_aligned_reads_df = pd.DataFrame(consensus_reads_list,
                                          columns=targets_in_dataset, index=file_names).reset_index()
consensus_aligned_reads_df.rename(columns={"index": "Sample Name"}, inplace=True)

#print(consensus_aligned_reads_df)

# Make CSV of Raw Consensus Reads for the DataBatch

os.chdir(output_subdirectory)

consensus_reads_subdir = "raw_consensus_aligned_reads_per_target"
try:
    os.mkdir(consensus_reads_subdir)
except FileExistsError:
    print("output folder for raw consensus reads frequency tables already exists")

os.chdir(output_subdirectory + consensus_reads_subdir)

print("", "Raw Consensus Reads Per Target", "", consensus_aligned_reads_df, "", sep="\n")
consensus_aligned_reads_df.to_csv(databatch[:-1] + "_raw_consensus_reads_per_target.csv", index=True)

print("", "CSV files saved", "Proceeding to Total Reads and Standard Deviation task", "", sep="\n")

# Make dataframe of total read for each sample and each gene target and dataframe for standard deviation of the total reads

total_reads_per_target_list = []
stdev_total_reads_per_target_list = []
for i in range(len(merged_dfs)):
    sample=merged_dfs[i]
    all_targets_tr=[]
    all_targets_std=[]
    for k in range(len(sample)):
        total_read_count=sample[k]["barcode_count"].sum()
        all_targets_tr.append(total_read_count)
        stdev_total_count=sample[k]["barcode_count"].std()
        all_targets_std.append(stdev_total_count)
    total_reads_per_target_list.append(all_targets_tr)
    stdev_total_reads_per_target_list.append(all_targets_std)

total_reads_count_df = pd.DataFrame(total_reads_per_target_list,
                                          columns=targets_in_dataset, index=file_names).reset_index()
total_reads_count_df.rename(columns={"index": "Sample Name"}, inplace=True)

stdev_total_reads_df = pd.DataFrame(stdev_total_reads_per_target_list,
                                          columns=targets_in_dataset, index=file_names).reset_index()
stdev_total_reads_df.rename(columns={"index": "Sample Name"}, inplace=True)

# Make CSV of Total Reads Per Target for the DataBatch

os.chdir(output_subdirectory)

total_reads_subdir = "total_unfiltered_reads_per_target"
try:
    os.mkdir(total_reads_subdir)
except FileExistsError:
    print("output folder for total reads tables already exists")

os.chdir(output_subdirectory + total_reads_subdir)

print("", "Total Reads Per Target", "", total_reads_count_df, "", sep="\n")
total_reads_count_df.to_csv(databatch[:-1] + "_total_unfiltered_reads_per_target.csv", index=True)

# Make CSV of Standard Deviation of Total Reads Per Target for the DataBatch

os.chdir(output_subdirectory)

stdev_total_reads_subdir = "stdev_total_unfiltered_reads_per_target"
try:
    os.mkdir(stdev_total_reads_subdir)
except FileExistsError:
    print("output folder for stdev total reads tables already exists")

os.chdir(output_subdirectory + stdev_total_reads_subdir)

print("", "Stdev for Total Reads Per Target", "", stdev_total_reads_df, "", sep="\n")
stdev_total_reads_df.to_csv(databatch[:-1] + "_stdev_total_unfiltered_reads_per_target.csv", index=True)

print("", "CSV files saved", "Proceeding to CpG loci exclusion task", "", sep="\n")


# Specific CpG Loci information within targeted amplicons -> For 300PE Sequencing -> Note, CpGs in primers removed

ankrd2_r1r2_300pe = [x - cpg_loci_change for x in [115, 141, 143, 157, 159, 163, 170, 184, 187, 198, 207, 225, 232, 241, 256,
                     261, 270, 275, 284, 297, 306, 329, 331, 339, 352, 362, 375, 380, 391]]
# ANRKD2 Nil exclusions

camk2b_r1r2_300pe = [x - cpg_loci_change for x in [27, 30, 32, 34, 36, 46, 63, 72, 82, 84, 87, 90, 100, 107, 113, 115, 117,
                     120, 123, 139, 149, 153, 159, 161, 163, 169, 174, 176, 184, 190, 213, 215,
                     217, 219, 242, 244, 251, 256, 264, 270, 272, 280, 282, 285, 403, 409, 423, 429,
                     432, 434, 440, 459, 463, 470, 473, 482, 495, 507, 509, 521, 523, 535, 538]]
# CAMK2B CpGs excluded in the 3' tail -> 291, 307, 312, 316, 323, 328, 337, 340, 362, 364, 366, 379, 390, 393

csrp3_r1r2_300pe = [x - cpg_loci_change for x in [90, 94, 98, 100, 104, 138, 143, 147, 154, 168, 188, 199, 219, 223, 228, 231,
                    255, 279, 294, 311, 317, 319, 341, 357, 416]]
# CSRP3 Nil exclusions

loc_l_r1r2_300pe = [x - cpg_loci_change for x in [44, 57, 66, 169, 181, 199, 272, 359, 417, 419,
                    423, 425, 439, 443, 452, 462, 469, 481, 483, 493, 500, 511]]
# LOC-L CpGs excluded in the 3' tail -> 274, 315, 321, 326, 328, 347

myod1_1_r1r2_300pe = [x - cpg_loci_change for x in [26, 31, 38, 46, 68, 70, 79, 82, 84, 94, 97, 103, 114, 117, 123, 132,
                      137, 139, 142, 145, 157, 160, 166, 169, 171, 180, 182, 196, 199, 213,
                      242, 245, 249, 263, 265, 271, 277, 293, 319]]
# 1-MYOD1 Nil exclusions

myod1_2_r1r2_300pe = [x - cpg_loci_change for x in [39, 48, 52, 65, 79, 92, 104, 109, 111, 117, 139, 145, 150, 152, 155, 166,
                      173, 176, 202, 206, 218, 223, 230, 238, 246, 254, 268, 270, 281, 287, 289,
                      302, 313, 323, 325, 328, 337, 343, 347, 351, 359, 369, 371, 373, 380, 394,
                      399, 424, 432, 446, 449, 453, 456, 464, 474, 476, 480, 483, 486]]
# 2-MYOD1 Nil exclusions

nkx2_1_r1r2_300pe = [x - cpg_loci_change for x in [31, 40, 43, 53, 64, 78, 110, 116, 160, 174, 177, 193, 203, 212, 234, 255,
                     264, 269, 275, 278,  461, 467, 503, 505, 516, 530, 533]]
# 1-NKX2-5 CpGs excluded in the 3' tail -> 295, 297, 310, 312, 314, 323, 335, 348, 351, 357, 370, 381, 387, 411

nkx2_2_r1r2_300pe = [x - cpg_loci_change for x in [35, 38, 41, 44, 47, 55, 63, 65, 82, 84, 110, 113, 120, 122, 132, 135, 153,
                     162, 171, 179, 186, 189, 192, 194, 225, 234, 237, 251, 257, 267, 270, 285,
                     294, 297, 300, 314, 327, 335, 347, 359, 373, 395, 397, 402, 404, 414, 443]]
# 2-NKX2-5 Nil exclusions

nppb_r1r2_300pe = [x - cpg_loci_change for x in [45, 62, 69, 78, 84, 92, 114, 117, 130, 135, 139, 153, 155, 157, 172, 175,
                   177, 184, 196, 208, 223, 232, 234, 240, 243, 250, 258, 285, 354, 384,
                   386, 399, 404, 409, 430, 440, 456, 459, 465, 467, 483, 504]]
# NPPB CpGs excluded in the 3' tail -> 312

pxdnl_r1r2_300pe = [x - cpg_loci_change for x in [26, 44, 96, 135, 141, 258, 422, 456]]
# PXDNL CpGs excluded in the 3' tail -> 305, 341, 353, 355, 367, 376


sbk2_r1r2_300pe = [x - cpg_loci_change for x in [29, 31, 134, 146, 169, 172, 187, 197, 336, 363, 365, 384, 397, 437, 444]]
# SBK2 CpG exclusions in the 3' tail -> 215, 258, 260, 262, 280, 303, 308, 317, 326

sbk3_1_r1r2_300pe = [x - cpg_loci_change for x in [26, 45, 64, 87, 157, 162, 173, 183, 190, 193, 223, 234, 238, 261, 266, 287,
                     304, 352, 378, 383, 388, 402, 410, 419]]
# 1-SBK3 CpGs excluded in the 3' tail -> 324, 335, 338

scn5a_r1r2_300pe = [x - cpg_loci_change for x in [34, 36, 40, 42, 66, 69, 72, 77, 80, 94, 112, 114, 118, 120, 126, 150, 170,
                    175, 178, 182, 201, 205, 209, 219, 223, 235, 237, 240, 250, 261, 281, 291,
                    298, 387, 390, 393, 396, 400, 416, 424, 430, 441, 453, 460, 462, 465, 470,
                    476, 478, 492, 501, 508, 512, 531, 541, 543, 547]]
# SCN5A CpGs excluded in the 3' tail -> 303, 309, 315, 321, 332, 338, 340, 343, 356, 364, 374, 376, 382

tnni3_1_r1r2_300pe = [x - cpg_loci_change for x in [40, 50, 74, 93, 132, 135, 141, 204, 209, 236, 240, 275,  430, 440, 464, 475]]
# 1-TNNI3 CpGs excluded in the 3' tail -> 301, 314, 324, 335, 339, 344, 355, 369, 379

tnni3_2_r1r2_300pe = [x - cpg_loci_change for x in [29, 51, 79, 84, 153, 177, 180, 195, 197, 208, 214, 218, 231, 258, 261,
                                                    324, 326, 343, 367, 388, 392, 400, 409]]
# 2-TNNI3 CpGs excluded in the 3' tail -> 306, 315

cpg_loci_300_all_targets = [ankrd2_r1r2_300pe,
                            camk2b_r1r2_300pe,
                            csrp3_r1r2_300pe,
                            loc_l_r1r2_300pe,
                            myod1_1_r1r2_300pe,
                            myod1_2_r1r2_300pe,
                            nkx2_1_r1r2_300pe,
                            nkx2_2_r1r2_300pe,
                            nppb_r1r2_300pe,
                            pxdnl_r1r2_300pe,
                            sbk2_r1r2_300pe,
                            sbk3_1_r1r2_300pe,
                            scn5a_r1r2_300pe,
                            tnni3_1_r1r2_300pe,
                            tnni3_2_r1r2_300pe]


# Exclude CpG Loci not formally covered by 300bp PE sequencing;
# if dataframe is empty, this is added to maintain the df number in each sample (can be fully excluded later). if
# if CpG loci missing - i.e. incorrect mapping of target within sample - then df is replaced by "empty_df" and Keyerror thrown
empty_df = pd.DataFrame(columns=["barcode_id"])
error_df_list = []
specific_300pe_dfs = []
def exclude_cpg_loci(list):
    index = 0
    error_msgs = []
    for data in list:
        pe300_dfs = []
        index += 1

        ankrd2_pivot = data[0]
        if ankrd2_pivot.empty == False:
            try:
                ankrd2_300pe_data = ankrd2_pivot[ankrd2_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(ankrd2_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in ANKRD2 Dataset Sample " + str(index),
                      data[0], "ANKRD2 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index-1],
                                                "Gene Target": "ANKRD2",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(ankrd2_pivot)

        camk2b_pivot = data[1]
        if camk2b_pivot.empty == False:
            try:
                camk2b_300pe_data = camk2b_pivot[camk2b_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(camk2b_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in CAMK2B Dataset Sample " + str(index),
                      data[1], "CAMK2B dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "CAMK2B",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(camk2b_pivot)

        csrp3_pivot = data[2]
        if csrp3_pivot.empty == False:
            try:
                csrp3_300pe_data = csrp3_pivot[csrp3_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(csrp3_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in CSRP3 Dataset Sample " + str(index),
                      data[2], "CSRP3 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "CSRP3",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(csrp3_pivot)

        loc_l_pivot = data[3]
        if loc_l_pivot.empty == False:
            try:
                loc_l_300pe_data = loc_l_pivot[loc_l_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(loc_l_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in LOC-L Dataset Sample " + str(index),
                      data[3], "LOC-L dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "LOC-L",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(loc_l_pivot)

        myod1_1_pivot = data[4]
        if myod1_1_pivot.empty == False:
            try:
                myod1_1_300pe_data = myod1_1_pivot[myod1_1_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(myod1_1_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in 1-MYOD1 Dataset Sample " + str(index),
                      data[4], "1-MYOD1 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "1-MYOD1",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(myod1_1_pivot)

        myod1_2_pivot = data[5]
        if myod1_2_pivot.empty == False:
            try:
                myod1_2_300pe_data = myod1_2_pivot[myod1_2_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(myod1_2_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in 2-MYOD1 Dataset Sample " + str(index),
                      data[5], "2-MYOD1 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "2-MYOD1",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(myod1_2_pivot)

        nkx2_1_pivot = data[6]
        if nkx2_1_pivot.empty == False:
            try:
                nkx2_1_300pe_data = nkx2_1_pivot[nkx2_1_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(nkx2_1_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in 1-NKX2-5 Dataset Sample " + str(index),
                      data[6], "1-NKX2-5 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "1-NKX2-5",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(nkx2_1_pivot)

        nkx2_2_pivot = data[7]
        if nkx2_2_pivot.empty == False:
            try:
                nkx2_2_300pe_data = nkx2_2_pivot[nkx2_2_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(nkx2_2_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in 2-NKX2-5 Dataset Sample " + str(index),
                      data[7], "2-NKX2-5 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "2-NKX2-5",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(nkx2_2_pivot)

        nppb_pivot = data[8]
        if nppb_pivot.empty == False:
            try:
                nppb_300pe_data = nppb_pivot[nppb_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(nppb_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in NPPB Dataset Sample " + str(index),
                      data[8], "NPPB dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "NPPB",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(nppb_pivot)

        pxdnl_pivot = data[9]
        if pxdnl_pivot.empty == False:
            try:
                pxdnl_300pe_data = pxdnl_pivot[pxdnl_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(pxdnl_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in PXDNL Dataset Sample " + str(index),
                      data[9], "PXDNL dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "PXDNL",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(pxdnl_pivot)

        sbk2_pivot = data[10]
        if sbk2_pivot.empty == False:
            try:
                sbk2_300pe_data = sbk2_pivot[sbk2_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(sbk2_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in SBK2 Dataset Sample " + str(index),
                      data[10], "SBK2 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "SBK2",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(sbk2_pivot)

        sbk3_1_pivot = data[11]
        if sbk3_1_pivot.empty == False:
            try:
                sbk3_1_300pe_data = sbk3_1_pivot[sbk3_1_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(sbk3_1_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in 1-SBK3 Dataset Sample " + str(index),
                      data[11], "1-SBK3 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "1-SBK3",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(sbk3_1_pivot)

        scn5a_pivot = data[12]
        if scn5a_pivot.empty == False:
            try:
                scn5a_300pe_data = scn5a_pivot[scn5a_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(scn5a_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in SCN5A Dataset Sample " + str(index),
                      data[12], "SCN5A dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "SCN5A",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(scn5a_pivot)

        tnni3_1_pivot = data[13]
        if tnni3_1_pivot.empty == False:
            try:
                tnni3_1_300pe_data = tnni3_1_pivot[tnni3_1_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(tnni3_1_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in 1-TNNI3 Dataset Sample " + str(index),
                      data[13], "1-TNNI3 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "1-TNNI3",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(tnni3_1_pivot)

        tnni3_2_pivot = data[14]
        if tnni3_2_pivot.empty == False:
            try:
                tnni3_2_300pe_data = tnni3_2_pivot[tnni3_2_r1r2_300pe].reset_index(drop=True)
                pe300_dfs.append(tnni3_2_300pe_data)
            except KeyError:
                print("", "KeyError - CpG loci missing in 2-TNNI3 Dataset Sample " + str(index),
                      data[14], "2-TNNI3 dataset replaced with empty dataframe", "", sep="\n")
                pe300_dfs.append(empty_df)
                error_msgs.append(pd.DataFrame({"Error": "Excl CpG Loci",
                                                "Sample Number": index,
                                                "Sample Name": file_names[index - 1],
                                                "Gene Target": "2-TNNI3",
                                                "Output": "Empty df"}, index=[0]))
        else:
            pe300_dfs.append(tnni3_2_pivot)

        specific_300pe_dfs.append(pe300_dfs)
    try:
        concat = pd.concat(error_msgs, ignore_index=True)
        error_df_list.append(concat)
    except ValueError:
        print("", "Excellent! - No error messages for 'Exclude CpG Loci'", "", sep="\n")

exclude_cpg_loci(pivotted_dfs)

# Test code:
#print(len(specific_300pe_dfs))
#test = specific_300pe_dfs[1]
#print(test[1])
#print(specific_300pe_dfs)

# Time to pre-process and clean-up data into new dataframes:
time_check1 = datetime.now()
df_made_time = time_check1 - starttime
print("Time taken for data processing into dataframes for analysis: " + str(df_made_time) + " seconds", "", sep="\n")

print("", "Done!", "Proceeding to correct and filter any missing values", "", sep="\n")

### PROCESSING FOR:
###         UNIQUE METHYLATION SEQUENCES DATAFRAMES
###         ABSOLUTE AND MEAN PAIRWISE DISTANCES TABLE
###         METHYLATION DENSITIES TABLE

# Keep sequences with NaN Values to prevent excessive filtering that can potentially skew the data
# use of the .groupby() function with flag "dropna=False") - doesnt seem to work
# NaN is of "Missing Completely at Random" (MCAR) - status
# -> this is certainly the case for NaN's within the main body of the amplicon, except the values quality trimmed at 3' area of R1/R2
# For methylation gradient / control DNA experiments where template should be either 100% Meth or 0% Unmeth -> take median of the sequence (row) to fill NaN
# For experiment proper samples i.e. barrett's, cell lines etc use median of the CpG locus
# Thus different Master_combined_python_scripts are required for each data cohort

# This next function ultimately:
#   i) removes/replaces "NaN" data points within and under certain thresholds set above
#   ii) converts and ensures data is astype integers to prep for UMS, Pairwise distance analysis and Methylation density analysis
#   iii) Fills Nan's will the consensus/majority value of a CpG Locus (median), assuming this is >80% of the calls (equates to a 2x SD of 0.8)


analysis_prepped_dfs = []
def analysis_prep(list):
    counter = -1
    for data in list:
        all_nan_excl = []
        counter += 1
        for i in range(len(data)):
            target = data[i]
            cpg_loci = cpg_loci_300_all_targets[i]
            df_t = target.T.reset_index()
            try:
                df_t.drop("cpg_locus", axis=1, inplace=True)
            except KeyError:
                print("cpg_locus header not found in " + file_names[counter] + targets_in_dataset[i])
            nans_num = df_t.isnull().sum(axis=1)
            row_or_col_num = df_t.shape[1]
            partial_nan_excl = []
            for j in range(len(nans_num)):
                if nans_num[j] / row_or_col_num <= nan_fraction_cpg_locus_exclusion_threshold and row_or_col_num >= 3:
                    partial_nan_excl.append(df_t.loc[j])
                else:
                    pass
            partial_nan_excl_df = pd.DataFrame(partial_nan_excl)
            nan_replaced = partial_nan_excl_df.apply(lambda row_col: row_col.fillna(
                round(random.random()) if mean_1 < row_col.mean() < mean_2 else row_col.median()),
                                                               axis=1)
            nan_replaced["cpg_locus"] = [cpg_loci[l] for l in partial_nan_excl_df.index]
            nan_replaced_pivot = nan_replaced.T
            integer_df = nan_replaced_pivot.dropna().astype(int).reset_index(drop=True)
            integer_df.columns = integer_df.iloc[-1]
            final_reorganised_df = integer_df.iloc[:-1, :].reset_index(drop=True)
            all_nan_excl.append(final_reorganised_df)
        analysis_prepped_dfs.append(all_nan_excl)

analysis_prep(specific_300pe_dfs)


print("", "Done!", "Proceeding to generate unique methylation sequences DataFrames", "", sep="\n")

# Add "id" column to analysis_prepped_dfs to allow a count of the final UMS

analysis_prepped_dfs_with_id = []
for data in analysis_prepped_dfs:
    new_id = []
    for target in data:
        target["id"] = list(range(1, len(target) + 1))
        new_id.append(target)
    analysis_prepped_dfs_with_id.append(new_id)

# Generate unique sequence dataframes and abundance of each sequence (methylation tag)
ums_dfs = []
def unique_methylation_sequences(list):
    index = 0
    error_msgs = []
    for data in list:

        index += 1
        target_ums = []

        try:
            ankrd2_300pe_ums_0nan = data[0].groupby([x for x in data[0].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(ankrd2_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group ANKRD2 dataset into UMS for Sample " + str(index) + ":",
                  data[0], sep="\n")
            target_ums.append(data[0])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "ANKRD2",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            camk2b_300pe_ums_0nan = data[1].groupby([x for x in data[1].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(camk2b_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group CAMK2B dataset into UMS for Sample " + str(index) + ":",
                  data[1], sep="\n")
            target_ums.append(data[1])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "CAMK2B",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            csrp3_300pe_ums_0nan = data[2].groupby([x for x in data[2].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(csrp3_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group CSRP3 dataset into UMS for Sample " + str(index) + ":",
                  data[2], sep="\n")
            target_ums.append(data[2])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "CSRP3",
                                            "Output": "Unfiltered Data"}, index=[0]))
        try:
            loc_l_300pe_ums_0nan = data[3].groupby([x for x in data[3].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(loc_l_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group LOC-L dataset into UMS for Sample " + str(index) + ":",
                  data[3], sep="\n")
            target_ums.append(data[3])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "LOC-L",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            myod1_1_300pe_ums_0nan = data[4].groupby([x for x in data[4].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(myod1_1_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group 1-MYOD1 dataset into UMS for Sample " + str(index) + ":",
                  data[4], sep="\n")
            target_ums.append(data[4])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "1-MYOD1",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            myod1_2_300pe_ums_0nan = data[5].groupby([x for x in data[5].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(myod1_2_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group 2-MYOD1 dataset into UMS for Sample " + str(index) + ":",
                  data[5], sep="\n")
            target_ums.append(data[5])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "2-MYOD1",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            nkx2_1_300pe_ums_0nan = data[6].groupby([x for x in data[6].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(nkx2_1_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group 1-NKX2-5 dataset into UMS for Sample " + str(index) + ":",
                  data[6], sep="\n")
            target_ums.append(data[6])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "1-NKX2-5",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            nkx2_2_300pe_ums_0nan = data[7].groupby([x for x in data[7].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(nkx2_2_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group 2-NKX2-5 dataset into UMS for Sample " + str(index) + ":",
                  data[7], sep="\n")
            target_ums.append(data[7])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "2-NKX2-5",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            nppb_300pe_ums_0nan = data[8].groupby([x for x in data[8].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(nppb_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group NPPB dataset into UMS for Sample " + str(index) + ":",
                  data[8], sep="\n")
            target_ums.append(data[8])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "NPPB",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            pxdnl_300pe_ums_0nan = data[9].groupby([x for x in data[9].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(pxdnl_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group PXDNL dataset into UMS for Sample " + str(index) + ":",
                  data[9], sep="\n")
            target_ums.append(data[9])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "PXDNL",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            sbk2_300pe_ums_0nan = data[10].groupby([x for x in data[10].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(sbk2_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group SBK2 dataset into UMS for Sample " + str(index) + ":",
                  data[10], sep="\n")
            target_ums.append(data[10])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "SBK2",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            sbk3_1_300pe_ums_0nan = data[11].groupby([x for x in data[11].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(sbk3_1_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group 1-SBK3 dataset into UMS for Sample " + str(index) + ":",
                  data[11], sep="\n")
            target_ums.append(data[11])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "1-SBK3",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            scn5a_300pe_ums_0nan = data[12].groupby([x for x in data[12].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(scn5a_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group SCN5A dataset into UMS for Sample " + str(index) + ":",
                  data[12], sep="\n")
            target_ums.append(data[12])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "SCN5A",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            tnni3_1_300pe_ums_0nan = data[13].groupby([x for x in data[13].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(tnni3_1_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group 1-TNNI3 dataset into UMS for Sample " + str(index) + ":",
                  data[13], sep="\n")
            target_ums.append(data[13])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "1-TNNI3",
                                            "Output": "Unfiltered Data"}, index=[0]))

        try:
            tnni3_2_300pe_ums_0nan = data[14].groupby([x for x in data[14].columns[:-1]]).id.count().reset_index().astype(int)
            target_ums.append(tnni3_2_300pe_ums_0nan)
        except KeyError:
            print("", "KeyError - Unable to group 2-TNNI3 dataset into UMS for Sample " + str(index) + ":",
                  data[14], sep="\n")
            target_ums.append(data[14])
            error_msgs.append(pd.DataFrame({"Error": "Unique Meth Seq",
                                            "Sample Number": index,
                                            "Sample Name": file_names[index - 1],
                                            "Gene Target": "2-TNNI3",
                                            "Output": "Unfiltered Data"}, index=[0]))

        ums_dfs.append(target_ums)
        print("Completed UMS for Sample " + str(index))
    try:
        concat = pd.concat(error_msgs, ignore_index=True)
        error_df_list.append(concat)
    except ValueError:
        print("", "Excellent! - No error messages for 'unique methylation sequences' generator", "", sep="\n")

unique_methylation_sequences(analysis_prepped_dfs_with_id)

#print(ums_dfs)
#print(len(ums_dfs))

print("", "Done!", "Saving as CSV Files...", "", sep="\n")

# Make the unique methylation sequences individual csv files for each gene target

os.chdir(output_subdirectory)

ums_subdir = "unique_meth_sequences_dataframes/"
try:
    os.mkdir(ums_subdir)
except FileExistsError:
    print("output folder for unique methylation sequences already exists")

os.chdir(output_subdirectory + ums_subdir)

ums_folder = "_ums_dataframes"
ums_file = "_unique_meth_sequences_for_"
for i in range(len(file_names)):
    try:
        os.mkdir(file_names[i] + ums_folder)
    except FileExistsError:
        print("Unique methylation sequences output folder for " + file_names[i] + " already exists")

def make_ums_csv_files (list):
    for i in range(len(list)):
        df_list=list[i]
        subsubdir = file_names[i] + ums_folder
        os.chdir(output_subdirectory + ums_subdir + subsubdir)
        for k in range(len(df_list)):
            df_list[k].to_csv(file_names[i] + ums_file + targets_in_dataset[k] + ".csv", index=False)

make_ums_csv_files(ums_dfs)


# Make Dataframe of UMI consensus Reads successfully aligned for each sample and gene target
ums_consensus_reads_list=[]
for i in range(len(ums_dfs)):
    ums_sample=ums_dfs[i]
    ums_all_targets=[]
    for k in range(len(ums_sample)):
        ums_target_length=len(ums_sample[k])
        ums_all_targets.append(ums_target_length)
    ums_consensus_reads_list.append(ums_all_targets)

ums_consensus_aligned_reads_df = pd.DataFrame(ums_consensus_reads_list,
                                          columns=targets_in_dataset, index=file_names)
ums_consensus_aligned_reads_df.index.names = ["Sample Name"]

# Make CSV of UMS Consensus Reads Frequency Per Target for the DataBatch

os.chdir(output_subdirectory)

ums_consensus_reads_subdir = "ums_consensus_aligned_reads_per_target"
try:
    os.mkdir(ums_consensus_reads_subdir)
except FileExistsError:
    print("output folder for UMS consensus reads frequency tables already exists")

os.chdir(output_subdirectory + ums_consensus_reads_subdir)

print("", "UMS Consensus Reads Per Target", "", ums_consensus_aligned_reads_df, "", sep="\n")
ums_consensus_aligned_reads_df.to_csv(databatch[:-1] + "_ums_consensus_reads_per_target.csv", index=True)

# Time to make UMS dataframes:
time_check2 = datetime.now()
ums_made_time = time_check2 - time_check1
print("Time taken to complete UMS dataframes: " + str(ums_made_time) + " seconds")


print("", "Files saved", "Proceeding to pairwise analysis", "", sep="\n")

# Pairwise (Hamming) distance Functions - output of a variable is delivered as a list to "all_hd" - these functions are for Dataframe inputs only
# Need to remove the "id" column again

pd_prepped_dfs = []
def pd_prep(list):
    for data in list:
        excl_nan = []
        for target in data:
            excl_nan.append(target.drop(["id"], axis=1).astype(int).reset_index(drop=True))
        pd_prepped_dfs.append(excl_nan)

pd_prep(analysis_prepped_dfs)

# Pairwise distance functions for INTRA-CRYPT pairwise distance

all_hd = []
def hammingdistance(single_pair):
    seq1 = single_pair[0]
    seq2 = single_pair[1]
    diff = 0
    length = len(seq1)
    for i in range(length):
        if seq1[i] != seq2[i]:
            diff += 1
    all_hd.append(diff)
    #individual_hd_lists(all_hd)


def get_pairs(reads_all_combinations):
    index = 0
    while index < len(reads_all_combinations):
        single_pair = reads_all_combinations[index]
        index += 1
        hammingdistance(single_pair)
        #print(single_pair)


def get_combinations(reads_list):
    reads_all_combinations = []
    for combs in combinations(reads_list, 2):
        reads_all_combinations.append(combs)
    get_pairs(reads_all_combinations)
    #print(reads_all_combinations)


def return_reads(df):
    reads_list = []
    for i in range(len(df)):
        reads_list.append(list(df.iloc[i]))
    get_combinations(reads_list)
    #print(reads_list)


def all_hamming_distances(df):
    return_reads(df)

# Loop to run a set of variables through the Pairwise Functions above - Output is to the list "hd_list" with each variable's
# total pairwise combinations as a separate list object (a list of lists)

pd_raw_outputs = []
pd_count = 0
for data in pd_prepped_dfs:
    sample_hd = []
    pd_count += 1
    for i in range(len(data)):
        if all_hd == []:
            all_hamming_distances(data[i])
            sample_hd.append(all_hd)
        else:
            all_hd = []
            all_hamming_distances(data[i])
            sample_hd.append(all_hd)
        all_hd = []
    pd_raw_outputs.append(sample_hd)
    print("INTRA-CRYPT raw pairwise distance combinations completed for sample " + str(pd_count))

#print(pd_raw_outputs)
#print(len(pd_raw_outputs))


# Time to make collate raw pd data:
time_check3 = datetime.now()
pd_made_time = time_check3 - time_check2
print("Time taken to complete pairwise distance combinations for all samples: " + str(pd_made_time) + " seconds")

# Pairwise distance processing of raw data
mean_target_pd = []
def get_mean_target_pw(list):
    for data in list:
        sample_pds = []
        for target in data:
            sample_pds.append(np.mean(target))
        mean_target_pd.append(sample_pds)

get_mean_target_pw(pd_raw_outputs)

#print(len(mean_target_pd))
#print(len(mean_target_pd[0]), len(mean_target_pd[1]), len(mean_target_pd[2]))
#print("", mean_target_pd, "", sep="\n")

print("", "Done!", "Saving results as CSV files...", "", sep="\n")

## Make the mean pairwise distances table csv file for the data batch being processed - THIS IS FOR INTRA-CRYPT Pairwise Distance
os.chdir(output_subdirectory)

pd_subdir = "mean_pairwise_distance_tables"
try:
    os.mkdir(pd_subdir)
except FileExistsError:
    print("output folder for mean pairwise distance already exists")

os.chdir(output_subdirectory + pd_subdir)

mean_target_pd_df = pd.DataFrame(mean_target_pd, columns=targets_in_dataset, index=file_names).reset_index()
mean_target_pd_df.index.names = ["Sample Name"]
print("", "Mean Pairwise Distances", "", mean_target_pd_df, "", sep="\n")
mean_target_pd_df.to_csv(databatch[:-1] + "_mean_intracrypt_pairwise_distances.csv", index=False)


# These are the pairwise distance functions for INTER-CRYPT pairwise comparisons - a "CARTESIAN PRODUCT"
# Note that the target e.g. myod or scn5a within each list of dfs now have uneven lengths of CpG loci as if a
# particular CpG locus did not meet the Nan replacement criteria then it is thrown out (usually CpG loci at the 3' end)
# So for the INTERCRYPT PWD Functions to work must standardise the dfs again to the minimum number of CpG loci per
# target pair.  e.g. B013 has 47 CpG Loci (columns) in SCN5A whereas B017 has 57 CpG Loci

all_hd_intercrypt = []
def hammingdistance_intercrypt(single_cartesian_pair):
    seq1 = single_cartesian_pair[0]
    seq2 = single_cartesian_pair[1]
    if len(seq1) == len(seq2):
        diff = 0
        length = len(seq1)
        for i in range(length):
            if seq1[i] != seq2[i]:
                diff += 1
        all_hd_intercrypt.append(diff)
    else:
        pass

def get_single_pair(reads_all_pairs):
    index = 0
    while index < len(reads_all_pairs):
        single_cartesian_pair = reads_all_pairs[index]
        index += 1
        #print(single_cartesian_pair)
        hammingdistance_intercrypt(single_cartesian_pair)


def get_cartesian_product(reads_list_df1, reads_list_df2):
    reads_all_pairs = []
    for pairs in product(reads_list_df1, reads_list_df2):
        reads_all_pairs.append(pairs)
    get_single_pair(reads_all_pairs)

def return_reads_intercrypt(df1, df2):
    reads_list_df1 = []
    reads_list_df2 = []
    for i in range(len(df1)):
        reads_list_df1.append(list(df1.iloc[i]))
    for i in range(len(df2)):
        reads_list_df2.append(list(df2.iloc[i]))
    get_cartesian_product(reads_list_df1, reads_list_df2)

def all_hamming_distances_intercrypt(df1, df2):
    #print("df1 and df2 being printed", df1, df2)
    return_reads_intercrypt(df1, df2)

all_intercrypt_combs = []
for inter_crypt_combs in combinations(pd_prepped_dfs, 2):
    all_intercrypt_combs.append(inter_crypt_combs)

print("", "Number of INTERcrypt combinations to analyse is: ", len(all_intercrypt_combs))

pd_raw_outputs_intercrypt = []
pd_icd_count = 0
for data in all_intercrypt_combs:
    pd_icd_count += 1
    sample_hd_intercrypt = []
    sample1 = data[0]
    sample2 = data[1]
    for i in range(len(sample1)):
        cpg_standardised_list = []
        for cpg in list(sample1[i].columns.values):
            try:
                fake_df = sample2[i][cpg]
                cpg_standardised_list.append(cpg)
            except KeyError:
                pass
        if all_hd_intercrypt == []:
            all_hamming_distances_intercrypt(sample1[i][cpg_standardised_list], sample2[i][cpg_standardised_list])
            sample_hd_intercrypt.append(all_hd_intercrypt)
        else:
            all_hd_intercrypt = []
            all_hamming_distances_intercrypt(sample1[i][cpg_standardised_list], sample2[i][cpg_standardised_list])
            sample_hd_intercrypt.append(all_hd_intercrypt)
        all_hd_intercrypt = []
    pd_raw_outputs_intercrypt.append(sample_hd_intercrypt)
    print("INTER-CRYPT raw pairwise distance combinations completed for sample " + str(pd_icd_count))

# Need to Label each INTER-CRYPT pair in the final CSV file, so generate new labels:

all_merged_intercrypt_labels = []
def merge_labels(single_label_pair):
    label1 = single_label_pair[0]
    label2 = single_label_pair[1]
    merged_labels = label1[:-11] + "vs" + label2[:-11]
    all_merged_intercrypt_labels.append(merged_labels)

def get_label_pair(all_intercrypt_labels):
    index = 0
    while index < len(all_intercrypt_labels):
        single_label_pair = all_intercrypt_labels[index]
        index += 1
        merge_labels(single_label_pair)

all_intercrypt_combs_labels = []
for label_combs in combinations(file_names, 2):
    all_intercrypt_combs_labels.append(label_combs)
get_label_pair(all_intercrypt_combs_labels)

print("", "INTERCRYPT PWD WORKED!", sep="\n")

print("", all_merged_intercrypt_labels, "", sep="\n")

# Pairwise distance processing of raw data - INTER-CRYPT MEANS

mean_icd_target_pd = []
def get_mean_icd_target_pw(list):
    for data in list:
        sample_pds = []
        for target in data:
            sample_pds.append(np.mean(target))
        mean_icd_target_pd.append(sample_pds)

get_mean_icd_target_pw(pd_raw_outputs_intercrypt)

## Make the mean pairwise distances table csv file for the data batch being processed - THIS IS FOR INTER-CRYPT Pairwise Distance
os.chdir(output_subdirectory)

pd_subdir = "mean_pairwise_distance_tables"
try:
    os.mkdir(pd_subdir)
except FileExistsError:
    print("output folder for mean pairwise distance already exists")

os.chdir(output_subdirectory + pd_subdir)

mean_icd_target_pd_df = pd.DataFrame(mean_icd_target_pd, columns=targets_in_dataset, index=all_merged_intercrypt_labels).reset_index()
mean_target_pd_df.index.names = ["Sample Name"]
print("", "Mean Pairwise Distances", "", mean_target_pd_df, "", sep="\n")
mean_icd_target_pd_df.to_csv(databatch[:-1] + "_mean_intercrypt_pairwise_distances.csv", index=False)


# Time to make pd_mean dataframe:
time_check4 = datetime.now()
pd_mean_time = time_check4 - time_check3
print("Time taken to calculate mean pairwise distance for targets for all samples: " + str(pd_mean_time) + " seconds")

print("", "Files saved", "Proceding to methylation density analysis", "", sep="\n")

#print(pd_prepped_dfs)
#test = pd_prepped_dfs[0]
#test2 = test[2]
#print(len(test2))

#Methylation density processing of the raw data
md_cpg_locus = []
md_amplicon = []
md_target = []
def get_meth_density(list):
    for data in list:
        cpg_locus = []
        amplicon = []
        target = []
        for i in range(len(data)):
            if len(data[i]) > 0:
                cpg_locus.append(data[i].mean())
                amplicon.append(data[i].mean(1))
                target.append(data[i].stack().mean())
            else:
                cpg_locus.append(pd.Series(np.nan, index=cpg_loci_300_all_targets[i]))
                amplicon.append(np.nan)
                target.append(np.nan)
        md_cpg_locus.append(cpg_locus)
        md_amplicon.append(amplicon)
        md_target.append(target)

get_meth_density(pd_prepped_dfs)

print("", "Done!", "Saving results as CSV files...", "", sep="\n")

#print("CpG Locus", md_cpg_locus, "", sep="\n")
#print("Amplicons", md_amplicon, "", sep="\n")
#print("Gene Targets", md_target, "", sep="\n")

# Time to make calculate meth density dataframe:
time_check5 = datetime.now()
md_time = time_check5 - time_check4
print("Time taken to calculate mean methylation density for targets for all samples: " + str(md_time) + " seconds")

os.chdir(output_subdirectory)

meth_den_subdir = "methylation_density_tables"
try:
    os.mkdir(meth_den_subdir)
except FileExistsError:
    print("output folder for methylation density tables already exists")

os.chdir(output_subdirectory + meth_den_subdir)

# Methylation density CSV for the entire target - per sample, per gene target
md_target_df = pd.DataFrame(md_target, columns=targets_in_dataset, index=file_names).reset_index()
md_target_df.index.names = ["Sample Name"]
print("", "Methylation Densities", "", md_target_df, "", sep="\n")
md_target_df.to_csv(databatch[:-1] + "_methylation_densities.csv", index=False)


# Make CSV files for Methylation density for individual CpG Loci on Individual Targets
flat_list_cpg_md = [target for sample in md_cpg_locus for target in sample]

for i in range(len(targets_in_dataset)):
    target = [flat_list_cpg_md[k] for k in range(i, len(flat_list_cpg_md), 15)]
    #list_of_lists = [a for a in target]
    target_df = pd.DataFrame(target, index=file_names).reset_index()
    target_df.to_csv(databatch[:-1] + "cpg_loci_methylation_densities_for_" + targets_in_dataset[i] + ".csv", index=False)

# Time to make calculate meth density dataframe:
time_check6 = datetime.now()
md_df_time = time_check6 - time_check5
print("Time taken to create mean methylation density dataframe: " + str(md_df_time) + " seconds")

print("", "Files saved", "Creating error messages DataFrame", "", sep="\n")

# Create a concatenated dataframe of all error messages for review later if necessary
os.chdir(output_subdirectory)

error_msgs_subdir = "error_messages"
try:
    os.mkdir(error_msgs_subdir)
except FileExistsError:
    print("output folder for error messages already exists")

os.chdir(output_subdirectory + error_msgs_subdir)

try:
    error_msgs_df = pd.concat(error_df_list, ignore_index=True)
    print("", "Error Messages", "", error_msgs_df, "", sep="\n")
    error_msgs_df.to_csv(databatch[:-1] + "_error_messages.csv", index=False)
except ValueError:
    print("", "Excellent! No error messages dataframe from these samples", "", sep="\n")

print("Processing complete!", "", sep="\n")

# Total time to run script:
total_script_time = datetime.now() - starttime
print("Total time taken to process " + str(len(file_names)) + " samples: " + str(total_script_time) + " seconds")

exit()