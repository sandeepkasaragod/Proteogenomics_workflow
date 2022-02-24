import os
import sys
from itertools import islice
from os.path import join as join
import pandas as pd

#step1: extracting the high binding affinity neoantigens
def variant_with_high_affinity(input_dir, output_file):
    write_file = open(output_file, "w") #variant_with_less_than_500.txt
    for lsts in os.listdir(input_dir): #List of neoantigen output variant
        for j in open(join(input_dir,lsts)):
            if " SB" in j.rstrip():
                ab = ""
                a = [x.rstrip() for x in j.replace('<=', ' ').split(' ') if len(x) != 0]
                for ls in a:
                    ab = ab + '\t' + ls
                if float(ab.split('\t')[-2]) <= 500.00:
                     write_file.write(ab.strip() + '\t' + lsts.split('.')[0].split('_')[0] + '\n')

    write_file.close()

#step2: extracting the entire wildtype peptides
def extract_wt_peps(input_dir, output_file):
    write_file = open(output_file, 'w') #WT_peptides.txt
    for lst_file in os.listdir(input_dir):
        for i in open(input_dir + "/" + lst_file):
            if i.startswith(' ') and not i.startswith(' Pos'):
                a = [x.rstrip() for x in i.replace('<=', ' ').split(' ') if len(x) != 0]
                ab = ""
                for ls in a:
                    ab = ab + '\t' + ls                
                write_file.write(ab.strip() + '\t' + lst_file.split('.')[0].split('_')[0] + '\n')

    write_file.close()


#step3: mapping the variant peptides corresponding wildtype peptides
def distc(a, b):
    f = (lambda x, y: 0 if x == y else 1)
    return sum(map(f, a, b))

def map_vars_wt(input_variant, input_wt, output_file):
    dicts_wt_pep = {}
    dicts_matched_pep = {}
    for i in open(input_wt):  #WT peptide input step2
        split_i = i.split('\t')
        if split_i[-1].rstrip() + "_" + split_i[1] + "_" + split_i[10] not in dicts_wt_pep:
            dicts_wt_pep[split_i[-1].rstrip() + "_" + split_i[1] + "_" + split_i[10]] = [i.rstrip()]
        else:
            dicts_wt_pep[split_i[-1].rstrip() + "_" + split_i[1] + "_" + split_i[10]].append(i.rstrip())
            #print (split_i[-1].rstrip() + "_" + split_i[1] + "_" + split_i[10])
            
    for j in open(input_variant): #variant peptide input step1
        split_j = j.split('\t')
        #print (split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('v', 'n'))
        if split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('v', 'n') in dicts_wt_pep:
            #if split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('V', 'W') == "BT20" + "_" + "HLA-B*15:01" + "_" + "W4":
            for k in dicts_wt_pep[split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('v', 'n')]:
                split_k = k.split('\t')
                if len(split_j[9].strip()) == len(split_k[9].strip()):
                    if  distc(split_j[9], split_k[9]) == 1:
                        if split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('v', 'n') not in dicts_matched_pep:
                            dicts_matched_pep[split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('v', 'n')] = j.rstrip() + '\t' + "WT" + '\t' + k.rstrip()
                        else:
                            aff = dicts_matched_pep[split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('v', 'n')].split('\t')[15]
                            if float(aff) > float(split_k[15]):
                                dicts_matched_pep[split_j[-1].rstrip() + "_" + split_j[1] + "_" + split_j[10].replace('v', 'n')] = j.rstrip() + '\t' + "WT" + '\t' + k.rstrip()

    write_file = open(output_file, 'w')
    for x, y in dicts_matched_pep.items():
        write_file.write(y + '\n')
    write_file.close()


#step4: replace SB to - in a file
def replace_sb(infile):
    write_file = open(infile.split('.')[0] + "_sb_replaced.txt", 'w')  #Neoantigens_with_wt-1.txt
    for i in open(infile): #Neoantigens_with_wt-1.txt'
        split_i = i.split('\t')
        if len(split_i) == 36:
            write_file.write(i.rstrip().replace(split_i[-1].rstrip(), "-") + '\t' + split_i[-1].rstrip() + '\n')
        else:
            write_file.write(i.rstrip() + '\n')

    write_file.close()

#Step5: Compare variant affinity with WT affinity
def compare_var_aff(infile, outfile):
    write_file = open(outfile, 'w')
    for i in open(infile):
        split_i = i.split('\t')
        if float(split_i[15]) > float(split_i[34]):
            write_file.write(i.rstrip() + '\n')
    write_file.close()

#step6: add header gene information to file
def gene_info(index_file):
    df_gene = []
    df_idx = []
    df = pd.read_csv(index_file, sep='\t', names=["ID", "Description"])
    return df

#gene_info("Input_peptides/variant/file_idx.idx")

def header_lbl():
    return ["Pos","MHC","Peptide","Core","Of","Gp","Gl","Ip","Il","Icore","ID","Score_EL","%Rank_EL","Score_BA","%Rank_BA",\
            "Aff(nM)","BindLevel","Sample","Pep_type","wt_Pos","wt_MHC","wt_Peptide","wt_Core","wt_Of","wt_Gp","wt_Gl","wt_Ip","wt_Il","wt_Icore","wt_Identity",\
            "wt_Score_EL","wt_%Rank_EL","wt_Score_BA","wt_%Rank_BA","wt_Aff(nM)","wt_BindLevel","Sample_wt"]

def prep_neoantigen_table(infile, index_file, outfile):
    header = header_lbl()
    gene_desc = gene_info(index_file)
    df = pd.read_csv(infile, sep='\t', names=header)
    merge = pd.merge(df, gene_desc, on="ID") 
    #df['Gene'] = 
    merge.to_csv(outfile, sep='\t')

if __name__ = "__main__":
    variant_dir = sys.argv[1]
    wt_dir = sys.argv[2]
    output_dir = sys.argv[3]
    header_idx = sys.argv[7]
    
    variant_with_high_affinity(variant_dir, join(output_dir, "neoantigens_with_less_than_500_aff.txt"))
    extract_wt_peps(wt_dir, join(output_dir, "wt_neoantigens.txt"))
    map_vars_wt(join(variant_dir, "neoantigens_with_less_than_500_aff.txt"), join(wt_dir, "wt_neoantigens.txt"), join(output_dir, "Variant_WT_combined.txt"))
    replace_sb(join(output_dir, "variant_wt_combined.txt")
    compare_var_aff(join(output_dir, "variant_wt_combined_sb_replaced.txt", output_dir, join(output_dir, "neoantigen_peptides.txt"))

    prep_neoantigen_table(join(output_dir, "neoantigen_peptides.txt"), header_idx, join(output_dir, "candidate_neoantigens.txt"))

#variant_with_high_affinity("Variant_results", "Variants_with_less_than_500_aff.txt")
#extract_wt_peps("Wildtype_results", "WT_peptides.txt")
#map_vars_wt("Variants_with_less_than_500_aff.txt", "WT_peptides.txt", "Variant_WT_combined.txt")
#replace_sb("Variant_WT_combined.txt")
#compare_var_aff("Variant_WT_combined_sb_replaced.txt", "Neoantigen_peptides.txt")
#prep_neoantigen_table("Neoantigen_peptides.txt", "Input_peptides/variant/file_idx.idx", "Candidate_neoantigens.txt")

