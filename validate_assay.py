import os
import argparse
import pandas as pd
import numpy as np

from ete3 import NCBITaxa, Tree, TreeStyle, faces, AttrFace, TextFace, NodeStyle

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', type=str,
                        help='datasets file')
    parser.add_argument('-p', '--pcr', type=str,
                        help='list of hit accessions')
    parse.add_argument('-t', '--target', type=int,
                        help= 'target taxid')
    parser.add_argument('-n', '--name', type=str,
                        help='assay name')

    args = parser.parse_args()
    data_file = args.data
    pcr_file = args.pcr
    target_taxid = args.target
    assay_name = args.name

    # run 

    # image title
    imagetitle = assay_name + " Visualization"
    # image filename
    imagefile = assay_name + "_Validation.png"
    # tsv filename
    csvfile = assay_name + "Validation_Table.tsv"

    # parse pcr file 
    pcr_acc = get_pcr_acc(pcr_file)
    # parse datasets file
    data = load_data(data_file, pcr_acc)
    # count data for taxids
    counts = count_data(data)
    # collapse data by species level 
    counts = collapse_data(counts)
    
    # create tree
    tree = create_tree(counts)

    # save outputs
    save_tree(tree, show, imagetitle, imagefile)
    save_tsv(counts, csvfile)



if __name__ == '__main__':
    main()

# returns list of accession numbers of hits from simPCR file 
def get_pcr_acc(file):
    # read sim_pcr data
    data = pd.read_csv(file, sep='\t')
    
    # get list of pcr hits
    acc = pd.DataFrame()
    acc['Accession'] = data['Full_Hit_ID'].apply(lambda r: str(r).split(' ')[0])
    acc.drop_duplicates(subset=["Accession"], keep='first')
    
    return acc["Accession"].tolist()

# return df with accession number, taxid, and if the strain was hit by assay for all potential assay targets 
# file: datatsets file 
# pcr_acc: list of hit accessions 
def load_data(file, pcr_acc):
    # read in datasets file
    data = pd.read_csv(file, sep='\t')
    
    # filter for only human hosts
    data = data[data['Host Taxonomic ID'] == 9606]
    # filter out missing collection dates 
    data = data[~data['Isolate Collection date'].isna()]
    
    # assign boolean if accession was hit 
    data['Hit'] = data['Accession'].apply(lambda x: x in pcr_acc)

    # return df w/ relevant columns
    return data[['Accession', 'Virus Taxonomic ID', 'Hit']].copy()

# return df summarizing the hits and totals for each taxid from input df 
# df: df from load_data, with rows of accession, taxid, and hit boolean 
def count_data(df):
    counts = df.groupby(['Virus Taxonomic ID'])['Hit'] \
                  .agg(['count', lambda x: (x == True).sum()]).reset_index()
    counts.columns = ['Virus Taxonomic ID', 'Total','Hits']
    return counts

def find_collapse_to(taxid, target):
    if not valid_taxid(taxid):
        return np.nan
    lineage = ncbi.get_lineage(taxid)
    index = find_target(lineage, target)
    if index == -1:
        return taxid 
    if index + 1 >= len(lineage):
        return taxid 
    return lineage[index + 1]

def collapse_data(df):
    '''Returns new dataframe collapsing at species level'''
    # add species taxid 
    df["Collapsed TaxID"] = df["Virus Taxonomic ID"].apply(collapse_to)
    
    # rows to remove
    rows_to_delete = []
    
    # parse through all rows in df 
    for index, row in df.iterrows():
        # if curr taxid =/= species
        if row['Collapsed TaxID'] != row['Virus Taxonomic ID']:
            # find row taxid == species
            matching_row = df[df['Virus Taxonomic ID'] == row['Collapsed TaxID']]
            # add curr row to be deleted later
            rows_to_delete.append(index)
            if not matching_row.empty:
                # if row exists already, add total&hits from curr taxid to species
                matching_index = matching_row.index[0]
                df.at[matching_index, 'Total'] += row['Total']
                df.at[matching_index, 'Hits'] += row['Hits']
            else:
                # if no row exists where taxid == species, create row and add to df 
                new_row = row.copy()
                new_row['Virus Taxonomic ID'] = row['Collapsed TaxID']
                new_row = pd.DataFrame([new_row], columns=df.columns)
                df = pd.concat([df, new_row], ignore_index=True)

    # Delete the non-matching rows and np.nan
    df.drop(rows_to_delete, inplace=True) 
    df.dropna(inplace=True)
    
    df['Virus Taxonomic ID'] = df['Virus Taxonomic ID'].astype(int)
    df['Total'] = df['Total'].astype(int)
    df['Hits'] = df['Hits'].astype(int)
    
    df['Virus'] = df["Virus Taxonomic ID"].apply(lambda x: ncbi.get_taxid_translator([x])[x])
    df = df.reset_index(drop=True)
    return df[["Virus", "Virus Taxonomic ID", "Total", "Hits"]]

# saves df to file with filename
def save_tsv(df, filename):
    df.to_csv(filename, sep='\t', index=False)

# Return (Tree) of all hits + direct descendants of target parent
# df: df with taxids, total strains, hit strains
def create_tree(df):
    
    # Create an empty dictionary for 'Total' and 'Hits'
    total_dict = {}
    hits_dict = {}
    
    # Iterate over the DataFrame and populate the dictionaries
    for index, row in df.iterrows():
        taxid = row['Virus Taxonomic ID']
        total = row['Total']
        hits = row['Hits']
        total_dict[taxid] = total
        hits_dict[taxid] = hits
        
    taxids = df['Virus Taxonomic ID'].values.tolist()

    tree = ncbi.get_topology(taxids, intermediate_nodes=True)
    
    # add hit, total counts to node
    add_counts(tree, hits_dict, "hits")
    add_counts(tree, total_dict, "total")
    
    return tree

# Adds information about number of strains hit to each node
# Returns None
# tree: (Tree) tree 
# counts: (dictionary) key: (int) taxID; value: (int) count
# name: (string) name of count
def add_counts(tree, counts, name):
    # traverse through nodes of tree
    for node in tree.traverse("postorder"):
        # if taxid has count info
        if node.taxid in counts:
            # add counts
            node.add_feature(name, counts[node.taxid])
        else: 
            # add field as nan
            node.add_feature(name, np.nan)

# show_tree style input 
# shows name, shows target in black, highlights hits in green(intended target)/yellow 
def show(node):
    # node styles
    target_nstyle = NodeStyle()
    target_nstyle["shape"] = "circle"
    target_nstyle["size"] = 7
    target_nstyle["fgcolor"] = "Black"
    
    # show name 
    face = AttrFace("sci_name", fsize = 12, fgcolor = 'Black')
    node.set_style(target_nstyle)
    
    # highlight hits
    if (not np.isnan(node.hits) and int(node.hits) > 0):
        face.background.color = "LightPink"
        
    faces.add_face_to_node(face, node, column=0, position="branch-right")
        
    # add total/hit counts if it exists 
    if (not np.isnan(node.total) and int(node.total) > 0):
        hits_face = TextFace('hits:' + str(node.hits), fsize = 12, fgcolor = 'Black')
        total_face = TextFace('total:' + str(node.total), fsize = 12, fgcolor = 'Black')
        total_face.margin_bottom = 10
        faces.add_face_to_node(hits_face, node, column = 0)
        faces.add_face_to_node(total_face, node, column = 0)

def save_tree(tree, style, title, filename):
    ts = TreeStyle()
    # hide default name of taxid 
    ts.show_leaf_name = False
    # use style 
    ts.layout_fn = style
    # add title
    ts.title.add_face(TextFace(title, fsize=20), column=0)
    # save as file 
    tree.render(filename, tree_style=ts)