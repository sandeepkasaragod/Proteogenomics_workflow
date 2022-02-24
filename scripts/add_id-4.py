import os
import sys
from itertools import islice
from os.path import join as join

cnt=1
def p_header(sep_field):
    header = ""
    for i in open('proteomics_data_header.txt'):
        header = header + sep_field + i.rstrip()
    return header 

def add_id(indir, sep_field, output_dir):
    prot_header = p_header(sep_field)
    for lst_file in os.listdir(indir):
        write_file = open(join(output_dir, lst_file), 'w')
        for i in islice(join(indir + lst_file), 0, None):
            split_i = i.split(sep_field)
            if split_i[0] == prot_header:
                write_file.write("Id" + sep_field + i.rstrip())
            else:
                cnt+=1
                write_file.write(str(cnt) + sep_field + i.rstrip())
    write_file.close()

if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    sep_field = sys.argv[3]
    create_id = add_id(input_dir, sep_field, output_dir)
    


