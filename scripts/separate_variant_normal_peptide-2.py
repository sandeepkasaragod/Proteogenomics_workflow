import pandas as pd
import os
import sys
from os.path import join as join
from itertools import islice

db_seq = ""

def load_db(db_path):
    for i in open(db_path):
        db_seq = i.rstrip()
    return db_seq

def seperator(sep_type):
    if sep_type == "tsv":
        return "\t"
    elif sep_type == "csv":
        return ","

def read_header(header_file):
    lst_header = []
    for i in open(header_file):
        lst_header.append(i.rstrip())
    return lst_header

def write_file(input_file, output_file_variant, output_file_normal, sep_field, header_lst, db_seqs):
    df = pd.read_csv(input_file, sep="\t")
    seq_idx = df.columns.get_loc('Sequence')
    header_raw = df.columns
    header = ""
    for i in header_raw:
        header = header + sep_field + i.strip()
    write_variant = open(output_file_variant, 'w')
    write_variant.write(header.strip() + '\n')
    
    write_normal = open(output_file_normal, 'w')
    write_normal.write(header.strip() + '\n')

    with open(join(input_file)) as f:
        for i in islice(f, 1, None):
            split_i = i.split('\t')
            if split_i[seq_idx].strip() in db_seqs:
                write_normal.write(i.rstrip() + '\n')
            else:
                write_variant.write(i.rstrip() + '\n')
    write_normal.close()
    write_variant.close()
    
if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file_variant = sys.argv[2]
    output_file_normal = sys.argv[3]
    sep_field = sys.argv[4]
    header_file = sys.argv[5]
    db = sys.argv[6]
    sep_type = seperator(sep_field)
    load_dbs = load_db(db)
    #for lst_file in os.listdir(input_dir):
    write_file(input_file, output_file_variant, output_file_normal, sep_type, header_file, load_dbs)
