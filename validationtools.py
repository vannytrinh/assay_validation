import pandas as pd
import numpy as np
from ete3 import NCBITaxa

ncbi = NCBITaxa()

# --- FUNCTIONS TO PARSE DATA ---

# date: date value of datasets file 
def datasets_get_year(date):
    '''
    Retrieve year from datasets file date elements 
    ex. '2000-01-09T00:00:00Z' --> '2000'
    '''
    if (date is not np.NaN):
        return int(str(date).split('-')[0])
    else:
        return np.NaN

# filepath: filepath to simulate_PCR file 
def get_accessions(filepath): 
    '''
    Returns list of hit accessions in simulate_PCR file 
    '''
    # read sim_pcr data
    data = pd.read_csv(filepath, sep='\t')

    # get list of pcr hits
    acc = pd.DataFrame()
    acc['Accession'] = data['Full_Hit_ID'].apply(lambda r: str(r).split(' ')[0])
    acc.drop_duplicates(subset=["Accession"], keep='first')

    return acc["Accession"].tolist()

# return df with accession number, taxid, and if the strain was hit by assay for all potential assay targets
# file: datatsets file
# pcr_acc: list of hit accessions
def load_data(file, pcr_acc):
    '''
    Parse dataset file and indicate if accession was hit by assay 
    '''
    # read in datasets file
    data = pd.read_csv(file, sep='\t')

    # filter for only human hosts
    data = data[data['Host Taxonomic ID'] == 9606]
    # filter out missing collection dates
    data = data[~data['Isolate Collection date'].isna()]
    # filter for complete sequences
    data = data[data['Completeness'] == 'COMPLETE']

    # assign boolean if accession was hit
    data['Hit'] = data['Accession'].apply(lambda x: x in pcr_acc)

    # retreive collection year
    data['Collection year'] = data['Isolate Collection date'].apply(datasets_get_year)
    data['Release year'] = data['Release date'].apply(datasets_get_year)

    # return df w/ relevant columns
    return data[['Accession', 'Virus Taxonomic ID', 'Collection year', 'Release year', 'Hit']].copy()

# --- FUNCTIONS TO ANALYZE DATA ---

# return df summarizing the hits and totals for each taxid from input df
# df: df from load_data, with rows of accession, taxid, and hit boolean
def count_data(df):
    '''
    Summarize total and hit data by taxid 
    '''
    counts = df.groupby(['Virus Taxonomic ID'])['Hit'] \
                  .agg(['count', lambda x: (x == True).sum()]).reset_index()
    counts.columns = ['Virus Taxonomic ID', 'Total','Hits']
    return counts

# returns true if it is a valid taxid for ncbi database
def valid_taxid(taxid):
    '''
    Returns true if taxid is valid in NCBI database     
    '''
    try:
        lineage = ncbi.get_lineage(taxid)
        return True  # TaxID exists
    except ValueError:
        print(f'Error: TaxID {taxid} not found')
        return False  # TaxID does not exist

# lineage: list of taxids 
# tagret: taxid of assay target 
def find_target(lineage, target):
    '''
    Returns the index of assay target taxid in the given linage. 
    Returns -1 if target is not in lineage 
    '''
    try:
        index = lineage.index(target)
        return index
    except ValueError:
        return -1

# taxid: taxid to find the ancestor of 
# target: target taxid of assay 
def find_collapse_to(taxid, target):
    '''
    Return the taxid of the ancestor of input taxid that is the direct child to target 
    '''
    # check if valid taxid was inputted 
    if not valid_taxid(taxid):
        return np.nan
    # get taxids lineage of input taxid 
    lineage = ncbi.get_lineage(taxid)
    # find where the target is in the lineage 
    index = find_target(lineage, target)
    # if taxid is not in lineage, do not collapse 
    if index == -1:
        return taxid
    # if taxid is in lineage but is last in lineage
    # return original taxid 
    if index + 1 >= len(lineage):
        return taxid
    # return taxid in lineage after target 
    return lineage[index + 1]


# df: outout from count_data or expand_data
# target_taxid: target taxid of assay 
def collapse_data(df, target_taxid):
    '''Returns new dataframe collapsing at species level'''

    # define collpase_to function for target_taxid 
    def collapse_to(taxid):
        return find_collapse_to(taxid, target_taxid)

    # define taxid to collapse data to  
    df["Collapsed TaxID"] = df["Virus Taxonomic ID"].apply(collapse_to)

    # rows to remove bc they have been added to a different level to collapse to
    rows_to_delete = []

    # parse through all rows in df
    for index, row in df.iterrows():
        # if curr taxid =/= taxid we want to collapse to 
        if row['Collapsed TaxID'] != row['Virus Taxonomic ID']:
            # find row of taxid we want to collapse to 
            matching_row = df[df['Virus Taxonomic ID'] == row['Collapsed TaxID']]
            # add curr row to be deleted later
            rows_to_delete.append(index)
            if not matching_row.empty:
                # if row exists already, add total&hits from curr taxid to species
                matching_index = matching_row.index[0]
                df.at[matching_index, 'Total'] += row['Total']
                df.at[matching_index, 'Hits'] += row['Hits']
            else:
                # if no row does not exist yet, create row and add to df
                new_row = row.copy()
                new_row['Virus Taxonomic ID'] = row['Collapsed TaxID']
                new_row = pd.DataFrame([new_row], columns=df.columns)
                df = pd.concat([df, new_row], ignore_index=True)

    # Delete the non-matching rows and np.nan
    df.drop(rows_to_delete, inplace=True)
    df.dropna(inplace=True)

    # make numbers ints 
    df['Virus Taxonomic ID'] = df['Virus Taxonomic ID'].astype(int)
    df['Total'] = df['Total'].astype(int)
    df['Hits'] = df['Hits'].astype(int)

    # add virus name 
    df['Virus'] = df["Virus Taxonomic ID"].apply(lambda x: ncbi.get_taxid_translator([x])[x])
    df = df.reset_index(drop=True)
    return df[["Virus", "Virus Taxonomic ID", "Hits", "Total"]]

# data: dataframe from count_data
# years: list of years
# target_taxid: taxid of assay target 
def expand_data(data, years, target_taxid):
    '''
    Expands given data so that it contains data from each given year 
    Data for a given year includes data collected up to that year 
    '''
    # collect data for each time point
    tps_data = []

    # iter through each input year 
    for yr in years:
        # count data for each taxid collected before or in given year
        counts = count_data(data[data["Collection year"] <= yr])
        # collapse data by species level
        counts = collapse_data(counts, target_taxid)
        counts.rename(columns={'Hits': 'Hits '+str(yr), 'Total': 'Total '+str(yr)}, inplace=True)
    
        tps_data.append(counts)

    # combine data for each timepoint into a combined dataframe
    combined_df = tps_data[0]
    for tp in tps_data[1:]:
        combined_df = pd.merge(combined_df, tp, on=['Virus', 'Virus Taxonomic ID'], how='outer')
    combined_df = combined_df.fillna(0)

    return combined_df

# hits: list of number of hits
# totals: list of number of totals
def get_ratio(hits, totals):
    '''
    Returns list of ratios
    Each hit is divided by corresponding total 
    '''
    # Make sure that both lists are equal lengths 
    if len(hits) != len(totals):
        print('Error: Unequal list lengths')
        return None

    # store ratios here 
    ratios = []

    # iter through values in lists 
    for elem1, elem2 in zip(hits, totals):
        if elem2 != 0:  # To avoid division by zero
            ratios.append(elem1 / elem2)
        else:
            # if denominator 0 
            ratios.append(np.nan)
    return ratios