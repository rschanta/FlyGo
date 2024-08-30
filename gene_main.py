from gene_module import prep_gene_association_file, parse_obo_file,process_file
import glob
import os 
#%% Paths to the gene association and obo files
gene_association_path = 'gene_association.csv'
obo_path = 'go-basic.obo'

# Get the data for the gene association and obo files
fb_data = prep_gene_association_file(gene_association_path)
obo_data = parse_obo_file('go-basic.obo')

###
#EDIT THESE 2 VARIABLES AS NEEDED
###
hyperoxia_dir = './Hyperoxia'
mutant_dir = './Mutant'


## Gets all the files
hyperoxia_files = glob.glob(os.path.join(hyperoxia_dir, '*.csv'))
mutant_files = glob.glob(os.path.join(mutant_dir, '*.csv'))
files = hyperoxia_files + mutant_files



# Loop through both
for file in files:
    process_file(file,fb_data,obo_data,)
