import pandas as pd
import sys


def col_index(infile, sep_type):
    df = pd.read_csv(infile, sep=sep_type)
    return (df['cosmic70'], df['CLINSIG'])

def read_file(infile, sep_type):
    idx = col_index(infile, sep_type)
    for i in open(infile):
        split_i = i.split(sep_type)
        if split_i[idx[0]] !='.' or split_i[idx[1]] !='.':
            write_file.write(i.rstrip() + '\n')
    write_file.close()



