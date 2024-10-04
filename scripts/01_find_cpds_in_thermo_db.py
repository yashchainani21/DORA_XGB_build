"""
Read in unique reactions from Brenda, KEGG, & MetaCyc then compile compounds that are already stored in eQuilibrator db
All compounds within these reactions would previously have had their SMILES canonicalized and sterochemistry removed
Compounds that already exist in eQuilibrator's compound db will be returned with a valid accession ID
Compounds that do not exist will have None for their accession ID
A list of dictionaries is finally obtained of the form [{'SMILES': 'cpd A smiles', 'accession ID': 'cpd A accession'}]
This final thermodynamics lookup table will be stored in data/interim/
"""

import time
import pickle
from thermo import thermo

with open('./postgres_conn_str.txt','r') as f:
    URI_EQ = f.read()

# select type of processed reactions to read in from "all_BKM_rxns", "all_AdH_mined_rxns", or "all_Monox_rxns"
dataset = "all_MetaCyc_V27_rxns"

processed_rxns_filepath = None
outfile_path = None

# note the different filepaths to load in data and store output data
if dataset == "all_BKM_rxns":
    processed_rxns_filepath = '../data/processed/all_processed_rxns_no_thermo.pkl'
    outfile_path = '../data/interim/all_processed_cpds_lookup_for_thermo.pkl'

if dataset == "all_AdH_rxns":
    processed_rxns_filepath = '../data/processed/all_AdH_mined_rxns_no_thermo.pkl'
    outfile_path = '../data/interim/all_AdH_cpds_lookup_for_thermo.pkl'

if dataset == "all_Monox_rxns":
    processed_rxns_filepath = '../data/processed/all_Monox_mined_rxns_no_thermo.pkl'
    outfile_path = '../data/interim/all_Monox_cpds_lookup_for_thermo.pkl'

if dataset == 'all_MetaCyc_V24_rxns':
    processed_rxns_filepath = '../data/processed/all_MetaCyc_V24_rxns_no_thermo.pkl'
    outfile_path = '../data/interim/all_MetaCyc_V24_cpds_lookup_for_thermo.pkl'

if dataset == 'all_MetaCyc_V26_rxns':
    processed_rxns_filepath = '../data/processed/all_MetaCyc_V26_rxns_no_thermo.pkl'
    outfile_path = '../data/interim/all_MetaCyc_V26_cpds_lookup_for_thermo.pkl'

if dataset == 'all_MetaCyc_V27_rxns':
    processed_rxns_filepath = '../data/processed/all_MetaCyc_V27_rxns_no_thermo.pkl'
    outfile_path = '../data/interim/all_MetaCyc_V27_cpds_lookup_for_thermo.pkl'

if processed_rxns_filepath:
    with open(processed_rxns_filepath, 'rb') as file:
        processed_rxns = pickle.load(file)

# Extract all unique SMILES strings from processed reactions file
all_cpds_list = []
for rxn in processed_rxns:
    for key in ['Substrates','LHS_cofactors','Products','RHS_cofactors']:
        for cpd_smiles in rxn[key]:
            if cpd_smiles not in all_cpds_list:
                all_cpds_list.append(cpd_smiles)

print(f"Total number of unique compounds:{len(all_cpds_list)}")

# Get Inchi keys for each unique SMILES string
all_cpds_dict = thermo.add_inchikeys_to_SMILESlist(all_cpds_list, verbose = True)

if __name__ == '__main__':

    chunk_size = 50
    max_workers = 4

    print(f"\nStarting to search eQuilibrator database to see if parsed compounds already exist")

    start_time = time.time()  # Record the start time
    processed_compounds = thermo.lookup_compounds_in_parallel(compounds = all_cpds_dict,
                                                                URI_EQ = URI_EQ,
                                                                chunk_size = chunk_size,
                                                                max_workers = max_workers)

    print(processed_compounds)

    with open(outfile_path, 'wb') as file:
        pickle.dump(processed_compounds,file)
    end_time = time.time()  # Record the end time

    # Calculate and print the time taken
    time_taken = end_time - start_time
    print(f"Time taken to run query eQuilibrator's compound database for all compounds: {time_taken} seconds")