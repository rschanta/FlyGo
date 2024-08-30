import pandas as pd
import requests
import os
import warnings
# Suppress specific warnings- The escape sequence is fine, and I'm using regular expressions :)
warnings.filterwarnings("ignore", category=SyntaxWarning, message="invalid escape sequence")
warnings.filterwarnings("ignore", category=UserWarning, message="This pattern is interpreted as a regular expression")

def prep_gene_association_file(file_path):
    '''
    INPUTS
        - file_path (string): Path to the gene association file

    OUTPUTS
        - filtered_data (pandas DataFrame): Data extracted from the gene association
            file, filtered to only include important info with the goal of 
            ultimately mapping gene names to GO_IDs
    PURPOSE
        - Finds biological process GO IDs, excludes any with 'NOT', and 
        reformats the alternate names column for easy searching later.

    '''
    data = pd.read_csv(file_path, sep='\t', header=None,skiprows=5)
    
    # Rename important columns
    data = data.rename(columns={1: 'FB_ID', 4: 'GO_ID',2:'name',3: 'keyword',5:'FB_ID2?',9:'full_name',10:'alt_names',8: 'type'})
    # Index out just those columns
    data = data[['FB_ID','GO_ID', 'name','alt_names', 'keyword', 'FB_ID2?', 'full_name','type']]
    
    # Only keep the biological Processes (P)
    filtered_data = data[data['type'] == 'P']
    
    # Filter out nots
    filtered_data = filtered_data[~filtered_data['keyword'].str.contains('NOT', case=False, na=False)]

    ## Reformat alternate name column for convenienve
    # Ensure string
    filtered_data['alt_names'] = filtered_data['alt_names'].astype(str)
    # Make sure begins and ends with |
    filtered_data['alt_names'] = filtered_data['alt_names'].apply(lambda x: f'|{x.strip("|")}|')
    
    return filtered_data


def get_info_name(filtered_data,name):
    '''
    INPUTS
        - filtered_data (pandas DataFrame): output of `prep_gene_association_file`
        - name (string): gene name

    OUTPUTS
        - name_data (pandas DataFrame): the filtered_data df filtered even further,
        to include just the specific gene name indicated.
    PURPOSE
        - Search to GO file for a specific gene

    '''
    # Find rows that have the input ID as its name, or in its list of alternate names
    alt_names_mask = filtered_data['alt_names'].str.contains(f'\|{name}\|', na=False)
    name_mask = filtered_data['name'].str.contains(name, na=False)
    combined_mask = alt_names_mask | name_mask
    name_data = filtered_data[combined_mask]
    
    return name_data


# Function to parse OBO file
def parse_obo_file(filename):
    '''
    INPUTS
        - filename (string): path to the .obo file from the DO website

    OUTPUTS
        - df (pandas DataFrame): all the information from the .obo file
            in an easier to read DataFrame
    PURPOSE
        - Parses out the info from the .obo file- note that some weird
        formatting in this one makes some of the more standard libraries 
        I tried failed, hence the manual parsing of parameters.

    '''
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    terms = []
    term = {}
    
    for line in lines:
        line = line.strip()
        if line.startswith('[Term]'):
            if term:
                terms.append(term)
            term = {}
        elif line.startswith('id:'):
            term['id'] = line.split('id: ')[1].strip()
        elif line.startswith('name:'):
            term['name'] = line.split('name: ')[1].strip()
        elif line.startswith('namespace:'):
            term['namespace'] = line.split('namespace: ')[1].strip()
        elif line.startswith('def:'):
            definition = line.split('def: ')[1].strip().strip('"')
            term['def'] = definition
        elif line.startswith('is_a:'):
            if 'is_a' not in term:
                term['is_a'] = []
            term['is_a'].append(line.split('is_a: ')[1].strip())
    
    # Append the last term
    if term:
        terms.append(term)
    
    # Convert to dataframe, rename some of the columns
    df = pd.DataFrame(terms)
    df = df.rename(columns={'id':'GO_ID','name':'name_obo'})
    return df

def get_name_info(fb_data,obo_data,name):
    '''
    INPUTS
        - fb_data (pandas DataFrame): output of `prep_gene_association_file`
        - obo_data (pandas DataFrame): output of `parse_obo_file`
        - name (string): gene name

    OUTPUTS
        - final_table (pandas DataFrame): dataframe merging together 
        relevant info from the obo file and the gene association file
    PURPOSE
        - This synthesizes the information from flybase and the larger obo
            file with the description of what things actuall are/do.

    '''
    # Find info for name in flybase data
    name_data = get_info_name(fb_data,name)
    # Merge with info in the obo file
    final_table = pd.merge(name_data, obo_data, on='GO_ID', how='left')
    # Index out just the columns 
    final_table = final_table[['FB_ID','GO_ID','name','name_obo']] 
    return final_table



def lambda_fun(name,fb_data,obo_data):
    '''
    INPUTS
        - fb_data (pandas DataFrame): output of `prep_gene_association_file`
        - obo_data (pandas DataFrame): output of `parse_obo_file`
        - name (string): gene name

    OUTPUTS
        - OBO_names (string): all of the names from the OBO file, concatenated as a string
            separated by |
        - FB_IDs (string): all of the FB_IDs, concatenated as a string
            separated by | . There really should only be one, but some have multiple.
        - GO_IDs (string): all of the relevant GO_IDs, concatenated as a string
            separated by |
        relevant info from the obo file and the gene association file
    PURPOSE
        - This lambda function is needed to add the new columns we want
        in Pandas without cumbersome loops.

    '''
    # Get data from a name
    try:
        parsed_data = get_name_info(fb_data,obo_data,name)
        # Get long strings, name, and FB_IDS
        OBO_names = parsed_data['name_obo'].str.cat(sep='|')
        FB_IDs = '|'.join(parsed_data['FB_ID'].unique())
        GO_IDs = parsed_data['GO_ID'].str.cat(sep='|')
        return OBO_names, FB_IDs,GO_IDs
    except:
        print(f'Problem with {name}! Data from FlyBase not found.')
        return None,None,None

#%% Read in files, name sure first column is named "gene"

def process_file(filepath,fb_data,obo_data,save=True):
    '''
    INPUTS
        - filepath (string): path to the .csv file
        - fb_data (pandas DataFrame): output of `prep_gene_association_file`
        - obo_data (pandas DataFrame): output of `parse_obo_file`

    OUTPUTS
        - down_regulated (pandas DataFrame)
        - up_regulated (pandas DataFrame)
        - both_regulated (pandas DataFrame)
    PURPOSE
        - This does the processing to check for significant p-values, determine
        whether a gene is down or up regulated, and then save these out to 
        an "Annotated" excel spreadsheet with 3 tabs, one for down, up, and both.
        Note that the file names/folders mirror the original.

    '''
    # Load in the csv, rename first column to gene
    file1 = pd.read_csv(filepath)
    file1 = file1.rename(columns={file1.columns[0]: 'gene'})
    
    # Remove NaNs
    file1 = file1.dropna(subset=['padj'])
    # Filter for significant p-value
    file1 = file1[file1['padj'] < 0.05]
    
    # Identify down-regulated entries
    down_regulated = file1[file1['log2FoldChange'] < -1]
    # Identify up-regulated entries
    up_regulated = file1[file1['log2FoldChange'] > 1]
    # Identified combined entries
    both_regulated = file1[(file1['log2FoldChange'] < -1) | (file1['log2FoldChange'] > 1)]
    
    # Make copies of dataframes since we need to modify them
    down_regulated = down_regulated.copy()
    up_regulated = up_regulated.copy()
    both_regulated = both_regulated.copy()
    
    # Apply mapping functions for the annotations
    down_regulated[['OBO_names', 'FB_IDs', 'GO_IDs']] = down_regulated['gene'].apply(lambda x: pd.Series(lambda_fun(x, fb_data,obo_data)))
    up_regulated[['OBO_names', 'FB_IDs', 'GO_IDs']] = up_regulated['gene'].apply(lambda x: pd.Series(lambda_fun(x, fb_data,obo_data)))
    both_regulated[['OBO_names', 'FB_IDs', 'GO_IDs']] = both_regulated['gene'].apply(lambda x: pd.Series(lambda_fun(x, fb_data,obo_data)))
    
    ## If save is true, save 
    if save == True:
        print('Save set to True (default). Saving Data')
        # Get directory, file name, and file name without extension
        directory_name = os.path.dirname(filepath)
        file_name_with_extension = os.path.basename(filepath)
        file_name, _ = os.path.splitext(file_name_with_extension)
        
        # Create new directory if it doesn't exist
        anno_dir = f'{directory_name}_ANNO'
        os.makedirs(anno_dir, exist_ok=True)
        
        # Create name of the file to mirror original name
        new_file_name = f'{anno_dir}/{file_name}_ANNO.xlsx'
        with pd.ExcelWriter(new_file_name) as writer:
            down_regulated.to_excel(writer, sheet_name='down_regulated', index=False)
            up_regulated.to_excel(writer, sheet_name='up_regulated', index=False)
            both_regulated.to_excel(writer, sheet_name='both', index=False)
        print(f'Data successfully processed and saved to {new_file_name}')
    return down_regulated, up_regulated,both_regulated

