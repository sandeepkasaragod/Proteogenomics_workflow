import sys
from os.path import join as join
import os
import pandas as pd

configfile: 'config.yml'

input_dir_p = config['input_dir_proteomics']
input_dir_g = config['input_dir_genomics']
separator = config['separator']
database = config['db']

output_dir = config['output_dir']
filter_col = config['output_dir'] + config['filter_col'] #step1
variant_lst = config['output_dir'] + config['separate_variant_normal']  + "variant" #step2
normal_lst = config['output_dir'] + config['separate_variant_normal']  + "normal"
variant_info = config['output_dir'] + config['pep_with_variant_info'] #step3
merge_files = config['output_dir'] + config['merge_files'] #step4
#novel_variants = config['output_dir'] + config['variants_not_in_cosmic_and_clinvar']
enricher = config['output_dir'] + config['enricher'] 
prep_for_neoantigen_prediction = config['output_dir'] + config['prep_for_neoantigen_prediction']
raw_neoantigens = config['output_dir'] + config['raw_neoantigen']
pep_len_netmhc = config['pep_len_netmhc']
netmhc_path = config['netmhc_path']
neoantigen = config['neoantigen']

os.makedirs(variant_lst, exist_ok=True)
os.makedirs(normal_lst, exist_ok=True)
os.makedirs(variant_info, exist_ok=True)
os.makedirs(join(prep_for_neoantigen_prediction, "variant"), exist_ok=True)
os.makedirs(join(prep_for_neoantigen_prediction, "wildtype"), exist_ok=True)
os.makedirs(join(raw_neoantigens, "variant"), exist_ok=True)
os.makedirs(join(raw_neoantigens, "wildtype"), exist_ok=True)

#https://stackoverflow.com/questions/65728953/conditional-execution-of-snakemake-rules-based-on-column-in-metatable

hlatype = pd.read_csv('neoantigen_lbl.txt', sep='\t').set_index("samples", drop = False)
print (hlatype.hla_type)
#data list
(proteomics, ) = glob_wildcards(join(config['input_dir_proteomics'], "{sample}.txt"))
(genomics, ) = glob_wildcards(join(config['input_dir_proteomics'], "{sample}.txt"))

rule all:
	input:
		expand(join(filter_col, "{sample}.txt"), sample=proteomics),
		expand(join(variant_lst, "{sample}.txt"), sample=proteomics),
		expand(join(variant_info, "{sample}.txt"), sample=proteomics),
		join(merge_files, "variant_merged.txt"),
		#join(enricher, "enrichment.txt"),
		join(prep_for_neoantigen_prediction, "variant", "file_idx.idx"),
		expand(join(raw_neoantigens, "variant", "{sample}.txt"), sample=proteomics),
		expand(join(neoantigen, "neoantigens_with_less_than_500_aff.txt"), sample=proteomics),

#Filter columns
rule step1:
	params:
		sep_field = separator,
		script_path = join("scripts", "filter_cols-1.py")
	input:
		txt = join(input_dir_p, "{sample}.txt"),
		header = join("meta_info", "proteomics_data_header.txt")
	output:
		txt = join(filter_col, "{sample}.txt")
	shell:
		"python {params.script_path} {input.txt} {output.txt} {params.sep_field} {input.header}"


#variant and normal peptide seperation	
rule step2:
	params:
		sep_field = separator,
		script_path = join("scripts", "separate_variant_normal_peptide-2.py"),
	input:
		txt = join(filter_col, "{sample}.txt"),
		header = join("meta_info", "proteomics_data_header.txt"),
		db = join("meta_info", "singleline_db.fasta")
	output:
		txt_variant = join(variant_lst, "{sample}.txt"),
		txt_normal = join(normal_lst, "{sample}.txt")
	shell:
		"python {params.script_path} {input.txt} {output.txt_variant} {output.txt_normal} {params.sep_field} {input.header} {input.db}"


#adding variant information
rule step3:
	params:
		sep_field = separator,
		script_path = join("scripts", "add_variant_information-3.py")
	input:
		header_p = join("meta_info", "proteomics_data_header.txt"),
		header_g = join("meta_info", "genomics_data_header.txt"),
		input_p = join(variant_lst, "{sample}.txt"),
		input_g = join(input_dir_g, "{sample}.txt"),
		db = join("meta_info", "singleline_db.fasta"),
	output:
		txt = join(variant_info, "{sample}.txt")
	shell:
		"python {params.script_path} {input.input_p} {input.input_g} {output.txt} {params.sep_field} {input.header_p} {input.header_g}"
	
#merge files
rule step4:
	params:
		sep_field = separator,
		script_path = join("scripts", "merge_files.py")
	input:
		txt_dir = variant_info
	output:
		txt = join(merge_files, "variant_merged.txt"),
		txt_gene = join(merge_files, "genes.txt")
	shell:
		"python {params.script_path} {input.txt_dir} {output.txt} {params.sep_field} {output.txt_gene}"

'''
#gene enrichment
rule step5:
	params:
		script_path = join("scripts", "enricher.r")
	input:
		txt = join(merge_files, "genes.txt")
	output:
		png = join(enricher, "enrichment.txt")
	shell:
		"Rscript {params.script_path}" 

'''
#prepare sequence file for neoantigen prediction
rule prepare_seq_for_neoantigen_pred:
	params:
		script_path = join("scripts", "prep_for_neoantigen.py"),
		sep_field = separator,
		output_dir_v = join(prep_for_neoantigen_prediction, "variant"),
		output_dir_n = join(prep_for_neoantigen_prediction, "wildtype"),
	input:
		txt = join(merge_files, "variant_merged.txt"),
		prot_db = database
	output:
		txt_v = join(prep_for_neoantigen_prediction, "variant", "file_idx.idx"),
		txt_n = join(prep_for_neoantigen_prediction, "wildtype", "file_idx.idx"),
	shell:
		"python {params.script_path} {input.txt} {params.output_dir_v} {params.output_dir_n} {params.sep_field} {input.prot_db}"


#neoantigen prediction
rule neoantigen_prediction:
        params:
                sep_field = separator,
                netmhc_dir = netmhc_path,
		mhc_lbl = lambda wildcards: hlatype.loc[wildcards.sample],
		pep_len = pep_len_netmhc
        input:
                fasta_v = join(prep_for_neoantigen_prediction, "variant", "{sample}.fasta"),
		fasta_n = join(prep_for_neoantigen_prediction, "wildtype", "{sample}.fasta")
	output:
		txt_v = join(raw_neoantigens, "variant", "{sample}.txt"),
		txt_n = join(raw_neoantigens, "wildtype", "{sample}.txt")
	shell:
		"{params.netmhc_dir} {input.fasta_v} -BA -l {params.pep_len} -a {params.mhc_lbl.hla_type} > {output.txt_v} && {netmhc_path} {input.fasta_n} -BA -l {params.pep_len} -a {params.mhc_lbl.hla_type} > {output.txt_n}"

rule candidate_neoantigens:
	params:
		script_path = "scripts/neoantigen_anlysis.py",
		output_dir = join(neoantigen),
		header_index = join(raw_neoantigens, "variant", "file_idx.idx"),
	input:
		variant_neoantigen = join(raw_neoantigens, "variant"),
		wt_neoantigen = join(raw_neoantigens, "wildtype"),
		#header_index = join(raw_neoantigens, "variant", "file_idx.idx")
	output:
		join(neoantigen, "neoantigens_with_less_than_500_aff.txt"),
		join(neoantigen, "wt_neoantigens.txt"),
		join(neoantigen, "variants_wt_combined.txt"),
		join(neoantigen, "variants_wt_combined_sb_replaced.txt"),
		join(neoantigen, "neoantigen_peptides.txt"),
		join(neoantigen, "candidate_neoantigens.txt"),
	shell:
		"{params.script_path} {input.variant_neoantigen} {input.wt_neoantigen} {params.output_dir} {params.header_index}"
