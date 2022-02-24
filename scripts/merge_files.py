import os
from itertools import islice
from os.path import join as join
import pandas as pd
import sys

def header(input_dir, sep_field):
    strs = ""
    lst = os.listdir(input_dir)
    df = pd.read_csv(join(input_dir, lst[0]), sep=sep_field)
    for i in df.columns:
        strs = strs + '"' +  i + '"' + sep_field
    return strs.strip()

def separator(sep_type):
    if sep_type == "tsv":
        return "\t"
    elif sep_type == "csv":
        return ","

def variant_merge(input_dir, output_file, sep_field):
    header_info = header(input_dir, sep_field)
    write_file = open(output_file, 'w')
    write_file.write('"ID"' + sep_field + header_info + sep_field + '"Sample_name"' + '\n')
    uniq_id = 1
    for list_file in os.listdir(input_dir):
        with open(join(input_dir, list_file)) as f:
            for i in islice(f, 1, None):
                write_file.write(str(uniq_id) + sep_field + i.rstrip() + sep_field + list_file.split('.')[0] + '\n')
                uniq_id+=1

    write_file.close()

def genes(infile, out_file, sep_field):
    df = pd.read_csv(infile, sep=sep_field)
    uniq_genes = []
    df = df['Gene'].drop_duplicates()
    df.to_csv(out_file, sep=sep_field, index=False)

if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    sep_type = separator(sys.argv[3])
    gene_file = sys.argv[4]
    variant_merge(input_dir, output_file, sep_type)
    genes(output_file, gene_file, sep_type)

