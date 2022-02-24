import os
import pandas as pd
import sys
from os.path import join as join
from itertools import islice
import read_file
import re

dicts_seq = {}

def read_prot_db(infile):
    rds = read_file.read_fasta(infile)
    for rows in rds:
        dicts_seq[rows[0].split('|')[3]] = rows[1].rstrip()

def header(infile, sep_field):
    header = ""
    with open(infile) as f:
        for i in islice(f, 0, 1):
            split_i = i.split(sep_field)
            for j in split_i:
                header = header + sep_field + j
        
    return header.strip() + sep_field + "Variant_peptide" + sep_field + "Normal_peptide"

def fetch_index(infile, sep_field):
    df = pd.read_csv(infile, sep="\t")
    return df.columns.get_loc('Sequence'), df.columns.get_loc('Master Protein Descriptions'), df.columns.get_loc('Sample_name')

def fetch_index1(infile, sep_field):
    df = pd.read_csv(infile, sep="\t")
    return df.columns.get_loc('Sequence'), df.columns.get_loc('Master Protein Descriptions'), df.columns.get_loc('Sample_name'), df.columns.get_loc('Variant_peptide'), df.columns.get_loc('Normal_peptide')

def seperator(sep_type):
    if sep_type == "tsv":
        return "\t"
    elif sep_type == "csv":
        return ","

def process(infile, sep_field, output_dir_v, output_dir_n):
    #print (infile, sep_field, output_dir)
    dicts_seq_v = {}
    dicts_seq_n = {}
    dicts_prot_info = {}
    header_idx = fetch_index1(infile, sep_field)
    with open(infile) as f:
        for i in islice(f, 1, None):
            split_i = i.split(sep_field)
            if split_i[header_idx[2]].rstrip() not in dicts_seq_v:
                dicts_seq_v[split_i[header_idx[2]].rstrip()] = [split_i[header_idx[3]]]
                dicts_prot_info[split_i[header_idx[2]].rstrip()] = [split_i[header_idx[1]]]
                dicts_seq_n[split_i[header_idx[2]].rstrip()] = [split_i[header_idx[4]]]
            else:
                dicts_seq_v[split_i[header_idx[2]].rstrip()].append(split_i[header_idx[3]])
                dicts_prot_info[split_i[header_idx[2]].rstrip()].append(split_i[header_idx[1]])
                dicts_seq_n[split_i[header_idx[2]].rstrip()].append(split_i[header_idx[4]])

                #store sample name as key and seq as value

    cnt = 1 
    write_idx = open(join(output_dir_v, 'file_idx.idx'), 'w')
    for k, v in dicts_seq_v.items():
        write_file = open(join(output_dir_v, k + '.fasta'), 'w')
        for peps_idx in range(len(v)):
            write_file.write('v' + str(cnt) + '\n' + v[peps_idx].rstrip() + '\n')
            write_idx.write('v' + str(cnt) + '\t' + dicts_prot_info[k][peps_idx] + '\n')
            cnt+=1
        write_file.close()
    write_idx.close()

    cnt = 1
    write_idx = open(join(output_dir_n, 'file_idx.idx'), 'w')
    for k, v in dicts_seq_n.items():
        write_file = open(join(output_dir_n, k + '.fasta'), 'w')
        for peps_idx in range(len(v)):
            write_file.write('n' + str(cnt) + '\n' + v[peps_idx].rstrip() + '\n')
            write_idx.write('n' + str(cnt) + '\t' + dicts_prot_info[k][peps_idx] + '\n')
            cnt+=1
        write_file.close()
    write_idx.close()

# +-15 peptide for neoantigen_prediction
def create_pep_for_neoantigen(infile, sep_field, db_path):
    write_file = open(infile.split('.')[0] + "_res.txt", 'w')
    header_idx = fetch_index(infile, sep_field)
    header_pg = header(infile, sep_field)
    write_file.write(header_pg + '\n')
    read_prot_db(db_path)
    with open(infile) as f:
        for j in islice(f, 1, None):
            split_j = j.split(sep_field)
            np = split_j[header_idx[1]].split('_')[0].replace('"', '') + "_" + split_j[header_idx[1]].split('_')[1]
            mut = split_j[header_idx[1]].split('_')[2].split('#')[0]
            ref = mut[0]
            alt = mut[-1]
            loc = re.findall('\d+', mut)
            #write_file.write(mut + '\t' + ref + '\t' +loc[0] + '\t' + alt + '\n')
            seq = dicts_seq[np]
            #write_file.write(seq + '\n')
            mut_seq = seq[:int(loc[0]) -1] + alt + seq[int(loc[0]):]
            if int(loc[0]) >= 15 and int(loc[0]) + 15 <= len(mut_seq):
                write_file.write(j.rstrip() + '\t' + mut_seq[(int(loc[0])) - 15: int(loc[0]) + 15] + "\t" + seq[(int(loc[0])) - 15: int(loc[0]) + 15] + '\n')
                #pass
            elif int(loc[0]) >= 15 and int(loc[0]) + 15 > len(mut_seq):
                write_file.write(j.rstrip() + '\t' + mut_seq[(int(loc[0])) - 15:] + "\t" + seq[(int(loc[0])) - 15:] + '\n')
                #pass
            elif int(loc[0]) < 15 and int(loc[0]) + 15 <= len(mut_seq):
                write_file.write(j.rstrip() + '\t' + mut_seq[:int(loc[0]) + 15] + "\t" + seq[:int(loc[0]) + 15] +'\n')
    write_file.close()

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_dir_v = sys.argv[2]
    output_dir_n = sys.argv[3]
    sep_field = sys.argv[4]
    db = sys.argv[5]
    sep_type = seperator(sep_field)
    #for lst_file in os.listdir(input_dir):
    create_pep_for_neoantigen(input_file, sep_type, db)
    process(input_file.split('.')[0] + "_res.txt", sep_type, output_dir_v, output_dir_n)

