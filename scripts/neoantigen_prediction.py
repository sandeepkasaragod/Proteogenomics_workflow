import pandas as pd
import os
import sys
from os.path import join as join

mhc = {}
def mhc_func(infile):
        for i in open(infile):
                split_i = i.split('\t')
                mhc[split_i[0]] = split_i[1]
        return mhc


def generate_command(input_dir, netmhc_path):
    for lst_file in os.listdir(input_dir):
        print (netmhc_path, join(input_dir, lst_file) -BA 8,9,10,11,12,13,14,15 -a 

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sep_field = sys.argv[3]
    header_file = sys.argv[4]
    sep_type = seperator(sep_field)
    #for lst_file in os.listdir(input_dir):
    write_file(input_file, output_file, sep_type, header_file)


