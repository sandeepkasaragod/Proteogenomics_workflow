import pandas as pd
import os
import sys
from os.path import join as join

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

def write_file(input_file, output_file, sep_field, header_file):
    df = pd.read_csv(join(input_file), sep=sep_field)
    #header = read_header(header_file, columns = read_header(header_file))
    df.to_csv(join(output_file), sep=sep_field, columns = read_header(header_file), index=False)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sep_field = sys.argv[3]
    header_file = sys.argv[4]
    sep_type = seperator(sep_field)
    #for lst_file in os.listdir(input_dir):
    write_file(input_file, output_file, sep_type, header_file)

#python filter_cols-1.py ./../Objective-2_article/November_03_2020/Raw_data/PeptideGroups_combined test tsv proteomics_data_header.txt
