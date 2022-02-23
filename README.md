This proteogenomics workflow utilized the results from genomics and proteomics to identify the cancer neoantigens. The workflow includes two main steps i) pre-processing step, which performs the data cleaning, maps the genomics results to its corresponding proteomics data, and merges the results to perform a ii) post-processing steps. The second step includes gene enrichment analysis and neoantigen prediction (netMHCPan). The workflow is well suited to run on Linux and Bash on Linux (Windows).

## Requirements
	- conda
	- netMHCPan 4.1
All the necessary packages are included in the environment.yml. 
In case of any problem in running netMHCPan4.1 (neoantigen prediction), the following packages and update can be installed manually. 
	-sudo apt-get update -y
	-sudo apt-get install csh
	-sudo apt-get install -y tcsh

### Conda channels
	- conda config --add channels defaults
	- conda config --add channels bioconda
	- conda config --add channels conda-forge

## Input files
### Genomics
The workflow requires Annotated variant results from the ANNOVAR tool. In case the user has raw NGS data, then the user can run CusVarDB tool to generate the variants (this includes a complete genomics pipeline to generate annotated VCF results)

### Proteomics
PeptigeGroups from the Proteome Discoverer tools will be the input from the proteomics end. 

## Setup
	- conda env create --file environment.yml
	- conda activate proteogenomics_workflow


