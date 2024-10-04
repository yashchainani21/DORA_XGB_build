import copy
import pickle
import time
from thermo import thermo

# pick from either 'all_BKM_rxns', 'all_AdH_rxns', 'all_Monox_rxns'
dataset = "all_MetaCyc_V27_rxns"
unique_rxns_list_filepath = None
cpds_lookup_table_filepath = None
outfile_path = None

with open('./postgres_conn_str.txt','r') as f:
    URI_EQ = f.read()
# set thermodynamic parameters
# 0.1 millimolar to 100 millimolar
lb = 1e-4
ub = 1e-1

if dataset == "all_BKM_rxns":
    unique_rxns_list_filepath = '../data/processed/all_processed_rxns_no_thermo.pkl'
    cpds_lookup_table_filepath = '../data/processed/all_processed_cpds_lookup_for_thermo.pkl'
    outfile_path = '../data/processed/all_processed_rxns_with_thermo.pkl'

if dataset == "all_AdH_rxns":
    unique_rxns_list_filepath = '../data/processed/all_AdH_mined_rxns_no_thermo.pkl'
    cpds_lookup_table_filepath = '../data/processed/all_AdH_cpds_lookup_for_thermo.pkl'
    outfile_path = '../data/processed/all_AdH_rxns_with_thermo.pkl'

if dataset == "all_Monox_rxns":
    unique_rxns_list_filepath = '../data/processed/all_Monox_mined_rxns_no_thermo.pkl'
    cpds_lookup_table_filepath = '../data/processed/all_Monox_cpds_lookup_for_thermo.pkl'
    outfile_path = '../data/processed/all_Monox_rxns_with_thermo.pkl'

if dataset == "all_MetaCyc_V24_rxns":
    unique_rxns_list_filepath = '../data/processed/all_MetaCyc_V24_rxns_no_thermo.pkl'
    cpds_lookup_table_filepath = '../data/processed/all_MetaCyc_V24_cpds_lookup_for_thermo.pkl'
    outfile_path = '../data/processed/all_MetaCyc_V24_rxns_with_thermo.pkl'

if dataset == 'all_MetaCyc_V26_rxns':
    unique_rxns_list_filepath = '../data/processed/all_MetaCyc_V26_rxns_no_thermo.pkl'
    cpds_lookup_table_filepath = '../data/processed/all_MetaCyc_V26_cpds_lookup_for_thermo.pkl'
    outfile_path = '../data/processed/all_MetaCyc_V26_rxns_with_thermo.pkl'

if dataset == 'all_MetaCyc_V27_rxns':
    unique_rxns_list_filepath = '../data/processed/all_MetaCyc_V27_rxns_no_thermo.pkl'
    cpds_lookup_table_filepath = '../data/processed/all_MetaCyc_V27_cpds_lookup_for_thermo.pkl'
    outfile_path = '../data/processed/all_MetaCyc_V27_rxns_with_thermo.pkl'

# read in unique reaction list
if unique_rxns_list_filepath and cpds_lookup_table_filepath:
    with open(unique_rxns_list_filepath,'rb') as file:
        unique_rxns_list = pickle.load(file)

    with open(cpds_lookup_table_filepath,'rb') as file:
        cpds_lookup_table = pickle.load(file)

rxns_list_w_accession_ids = []

# initialize a counter to track the number of reaction strings for which all compounds have accession ids
# for these reactions, we can attempt to calculate their dG and MDF values
all_accession_ids_present_count = 0

# initialize a counter to track the number of reaction strings for which at least one compound has no accession
# for these reactions, we cannot attempt to calculate their dG and MDF values
not_all_accession_ids_present_count = 0

for rxn_entry in unique_rxns_list:
    # for the reaction equation in each reaction entry (which is a dict),
    # rewrite this in terms of accession IDs rather than in terms of SMILES
    rewritten_rxn_eq = thermo.rxn_eq_in_accession_ids(rxn_eq = rxn_entry['Reaction eq'],
                                                      cpds_lookup_table = cpds_lookup_table,
                                                      delimiter = ' = ')
    if rewritten_rxn_eq:
        all_accession_ids_present_count += 1
    else:
        not_all_accession_ids_present_count += 1

    # append the rewritten reaction string in terms of accession IDs as a new key-value pair in dictionary copy
    rxn_entry_copy = copy.copy(rxn_entry)
    rxn_entry_copy['Reaction eq in accession IDs'] = rewritten_rxn_eq
    rxns_list_w_accession_ids.append(rxn_entry_copy)

bad_cpds = ['CHEBI:48085']
for rxn in rxns_list_w_accession_ids:
    if rxn['Reaction eq in accession IDs']:
        if 'CHEBI:48085' in rxn['Reaction eq in accession IDs']:
            # reset value
            rxn['Reaction eq in accession IDs'] = None
            not_all_accession_ids_present_count += 1
            all_accession_ids_present_count -= 1

print(f"\nNumber of {dataset} reactions for which we can attempt to calculate reaction thermodynamics: {all_accession_ids_present_count}")

print(f"\nNumber of {dataset} reactions for which we cannot attempt to calculate reaction thermodynamics: {not_all_accession_ids_present_count}")

if __name__ == '__main__':

    start_time = time.time()

    rxns_w_thermo = thermo.calc_rxn_dG_and_MDF_in_parallel(rxns_list = rxns_list_w_accession_ids,
                                                           URI_EQ = URI_EQ,
                                                           chunk_size = 500,
                                                           max_workers = 4)

    end_time = time.time()

    # Calculate and print the time taken
    time_taken = end_time - start_time
    print(f"Time taken to calculate reaction dG and MDF values: {time_taken} seconds")

    # Save the updated list of reactions
    with open(outfile_path,'wb') as file:
        pickle.dump(rxns_w_thermo, file)