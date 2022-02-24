import pandas as pd
import os
import sys
from os.path import join as join
from itertools import islice

nm_np = {}
def get_allele_freq(label, value):
    if label == "GT:AB:AD:DP:GQ:PL":
            #0/1:0.450:9,11:20:99:361,0,283
            get_ad = value.split(':')[2].split(',')[1] # get the alternate allele for ex 9, 11(alt_allele_)
            get_dp = value.split(':')[3]
            allele_freq = float(get_ad)/float(get_dp) #AD/DP
            return allele_freq, get_ad, get_dp

    elif label == "GT:AD:DP:GQ:PL":
            #1/1:0,37:37:99:1399,108,0
            get_ad1 = value.split(':')[1].split(',')[1]
            get_dp1 = value.split(':')[2]
            allele_freq1 = float(get_ad1)/float(get_dp1)
            return allele_freq1, get_ad1, get_dp1

def proteomics_header(header_p, sep_field):
    header = ""
    for i in open(header_p):
        header = header + sep_field + i.strip()
    return header.strip()

def genomics_header(header_g, sep_field):
    header = ""
    for i in open(header_g):
        header = header + sep_field + i.rstrip()
    return header.strip()

def seperator(sep_type):
    if sep_type == "tsv":
        return "\t"
    elif sep_type == "csv":
        return ","

def load_nm_np(db_path):
    for i in open(db_path):
        split_i = i.split('\t')
        nm_np[split_i[1].rstrip()] = split_i[0]

load_nm_np(join('meta_info', 'NP_and_NM.txt'))

def vcf_header_idx(infile, sep_field):
    df = pd.read_csv(infile, sep=sep_field)
    return df.columns.get_loc('AAChange.refGene'), df.columns.get_loc('Genotype'), df.columns.get_loc('GT_value')
    
def read_vcf_file(infile, sep_field):
    dicts_vcfs = {}
    aa_change = vcf_header_idx(infile, sep_field)
    with open(infile) as f:
        for i in islice(f, 1, None):
            split_i = i.split(sep_field)
            if split_i[aa_change[0]] != '.' and split_i[aa_change[0]] != "UNKNOWN":
                for split_vars in split_i[aa_change[0]].split(','):
                    split_var_info = split_vars.split(':')
                    dicts_vcfs[split_var_info[1] + "_" + split_var_info[-1][2:]] = i.rstrip()
    return dicts_vcfs

def pep_header_idx(pep_file, sep_field):
    df = pd.read_csv(pep_file, sep=sep_field)
    return df.columns.get_loc('Master Protein Descriptions')

def read_peptide_file(input_file_p, input_file_g, output_file, sep_field, header_p, header_g):
    #create a function to read selected columns, peptide and genomics
    vcf_header = vcf_header_idx(input_file_g, sep_field)
    dicts_vcf = read_vcf_file(input_file_g, sep_field)
    #p_header = proteomics_header(header_p, sep_field)
    prot_acc_idx = pep_header_idx(input_file_p, sep_field)
    p_header = proteomics_header(header_p, sep_field)
    g_header = genomics_header(header_g, sep_field)
    write_file = open(output_file, 'w')
    write_file.write(p_header + sep_field + g_header + sep_field + 'Allele_freq(AD/DP)' + sep_field + 'Allele_depth(AD)' + sep_field + 'Total_depth(DP)' + sep_field + 'Gene' + '\n')
    with open(input_file_p) as f:
        for j in islice(f, 1, None):
            split_j = j.split('\t')
            split_var_info = split_j[prot_acc_idx].split(';')[0].split('#')
            acc = split_var_info[0].split('_')[0] + "_" + split_var_info[0].split('_')[1]
            #print (split_var_info[0], len(split_var_info[0].split('_')))
            if len(split_var_info[0].split('_')) > 2:
                mut = split_var_info[0].split('_')[2]
                gene = split_var_info[1]
                if acc.replace('"', '') in nm_np:
                    if nm_np[acc.replace('"', '')].split('.')[0] + "_" + mut in dicts_vcf:
                        #print (nm_np[acc.replace('"', '')].split('.')[0] + "_" + mut)
                        matched_vcf_info = dicts_vcf[nm_np[acc.replace('"', '')].split('.')[0] + "_" + mut]
                        a_freq = get_allele_freq(matched_vcf_info.split('\t')[vcf_header[1]], matched_vcf_info.split('\t')[vcf_header[2]])
                        write_file.write(j.rstrip() + sep_field + matched_vcf_info + sep_field + str(a_freq[0]) + sep_field + str(a_freq[1]) + sep_field + str(a_freq[2]) + sep_field + gene + '\n')

    write_file.close()

if __name__ == "__main__":
    input_file_p= sys.argv[1]
    input_file_g = sys.argv[2]
    output_file = sys.argv[3]
    sep_field = sys.argv[4]
    header_p = sys.argv[5]
    header_g = sys.argv[6]
    #db = sys.argv[7]
    sep_type = seperator(sep_field)
    #load_dbs = load_db(db)
    #for lst_file in os.listdir(input_dir):
    read_peptide_file(input_file_p, input_file_g, output_file, sep_type, header_p, header_g)
