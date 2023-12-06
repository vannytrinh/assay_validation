import pandas as pd
import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
                        help='simulate_pcr filename')
    parser.add_argument('-o', '--output', type=str,
                        help='output filename')

    # parse arguments 
    args = parser.parse_args()
    filepath = args.input
    output = args.output

    # run
    accs = parse_file(filepath)
    save_file(accs, output)


def parse_file(filepath):
    '''
    Returns list of accessions numbers that hit in simulate_pcr output file 
    '''
    # Check if the input file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"The file '{filepath}' does not exist.")

    # read sim_pcr data
    data = pd.read_csv(filepath, sep='\t')
    acc = pd.DataFrame()
    # extract the accession from the Full_Hit_ID column
    acc['Accession'] = data['Full_Hit_ID'].apply(lambda r: str(r).split(' ')[0])
    # remove duplicate accessions
    acc.drop_duplicates(subset=["Accession"], keep='first')
        
    return acc["Accession"].tolist()

def save_file(accessions, output):
    '''
    Saves a .txt file of accession numbers 
    '''
    # Open the file in write mode
    with open(output, 'w') as file:
        # Write each string to the file
        for acc in accessions:
            file.write("%s\n" % acc)

if __name__ == '__main__':
    main()