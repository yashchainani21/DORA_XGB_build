"""
This script processes compound thermodynamics information following an initial lookup table created previously by
01_find_cpds_in_thermo_db.py. It reads the lookup table stored in 'data/interim', which contains a list of dictionaries,
each representing a single compound with the fields 'SMILES', 'inchi_key', and 'accession ID'. This script then filters
out compounds with incomplete SMILES strings (as denoted by asterisks) or compounds containing elements other than
C, O, N, S, P, and H. With the compounds still left that have valid SMILES strings, we then extract out compounds for
which their accession IDs in eQuilibrator are none, i.e. these compounds have not been decomposed and added to the db.
As such, these compounds with valid SMILES strings but that lack accession ids are added to eQuilibrator
The final lookup table is then saved in data/processed/
"""

from thermo import thermo
import pickle

with open('./postgres_conn_str.txt','r') as f:
    URI_EQ = f.read()

# Read in compounds that have already been queried against the eQuilibrator database
# pick from either 'all_BKM_rxns', 'all_AdH_rxns', 'all_Monox_rxns'
dataset = "all_MetaCyc_V27_rxns"

existing_thermo_lookup_filepath = None
outfile_path_for_valid_cpds = None
outfile_path_for_invalid_cpds = None

if dataset == "all_BKM_rxns": # all BRENDA, KEGG, and MetaCyc reactions
    existing_thermo_lookup_filepath = '../data/interim/all_processed_cpds_lookup_for_thermo.pkl'
    outfile_path_for_valid_cpds = '../data/processed/all_processed_cpds_lookup_for_thermo.pkl'
    outfile_path_for_invalid_cpds = '../archives/invalid_SMILES_from_all_BKM_rxns.pkl'

if dataset == "all_AdH_rxns": # all alcohol dehydrogenase reactions
    existing_thermo_lookup_filepath = '../data/interim/all_AdH_cpds_lookup_for_thermo.pkl'
    outfile_path_for_valid_cpds = '../data/processed/all_AdH_cpds_lookup_for_thermo.pkl'
    outfile_path_for_invalid_cpds = '../archives/invalid_SMILES_from_all_AdH_cpds.pkl'

if dataset == "all_Monox_rxns": # all monooxygenase reactions
    existing_thermo_lookup_filepath = '../data/interim/all_Monox_cpds_lookup_for_thermo.pkl'
    outfile_path_for_valid_cpds = '../data/processed/all_Monox_cpds_lookup_for_thermo.pkl'
    outfile_path_for_invalid_cpds = '../archives/invalid_SMILES_from_all_Monox_cpds.pkl'

if dataset == "all_MetaCyc_V24_rxns":
    existing_thermo_lookup_filepath = '../data/interim/all_MetaCyc_V24_cpds_lookup_for_thermo.pkl'
    outfile_path_for_valid_cpds = '../data/processed/all_MetaCyc_V24_cpds_lookup_for_thermo.pkl'
    outfile_path_for_invalid_cpds = '../archives/invalid_SMILES_from_all_MetaCyc_V24_cpds.pkl'

if dataset == 'all_MetaCyc_V26_rxns':
    existing_thermo_lookup_filepath = '../data/interim/all_MetaCyc_V26_cpds_lookup_for_thermo.pkl'
    outfile_path_for_valid_cpds = '../data/processed/all_MetaCyc_V26_cpds_lookup_for_thermo.pkl'
    outfile_path_for_invalid_cpds = '../archives/invalid_SMILES_from_all_MetaCyc_V26_cpds.pkl'

if dataset == 'all_MetaCyc_V27_rxns':
    existing_thermo_lookup_filepath = '../data/interim/all_MetaCyc_V27_cpds_lookup_for_thermo.pkl'
    outfile_path_for_valid_cpds = '../data/processed/all_MetaCyc_V27_cpds_lookup_for_thermo.pkl'
    outfile_path_for_invalid_cpds = '../archives/invalid_SMILES_from_all_MetaCyc_V27_cpds.pkl'

if existing_thermo_lookup_filepath:
    with open(existing_thermo_lookup_filepath, 'rb') as file:
        thermo_cpds_lookup = pickle.load(file)

# For each compound that was queried against the eQuilibrator database
# remove compounds that have an asterisk in their SMILES
# and that contain atoms other than carbon, oxygen, nitrogen, sulfur, phosphorus, and hydrogen
# these compounds are removed EVEN IF they have accession IDs
# because the accession ID might incorrectly point to a different compound
# for instance, SMILES that contain to asterisk return an accession that actually points to water
allowed_atoms = {'C', 'O', 'N', 'S', 'P', 'H'}
filtered_compounds, removed_compounds = thermo.filter_compounds(thermo_cpds_lookup, allowed_atoms = allowed_atoms)

# Retain compounds within the filtered compounds list that DO have an accession ID
compounds_with_accession = [compound for compound in filtered_compounds if compound['accession'] is not None]

# Extract out compounds within the filtered compounds list that DO NOT have an accession ID
compounds_with_no_accession = [compound for compound in filtered_compounds if compound['accession'] is None]

print(f"\nTotal number of unique compounds across all {dataset}: {len(thermo_cpds_lookup)}")

print(f"\nNumber of filtered compounds with valid SMILES strings and an existing accession in eQuilibrator: "
      f"{len(compounds_with_accession)}")

print(f"\nNumber of filtered compounds with valid SMILES strings but no existing accession in eQuilibrator: "
      f"{len(compounds_with_no_accession)}")

print(f"\nNumber of removed compounds due to invalid SMILES strings: "
      f"{len(removed_compounds)}")

print(f"\nSaving these invalid compounds in: {outfile_path_for_invalid_cpds}")
with open(outfile_path_for_invalid_cpds,'wb') as file:
    pickle.dump(removed_compounds, file)

if __name__ == "__main__":
    # add them to eQuilibrator next using lcp.get_compounds
    newly_compounds_added = thermo.insert_compounds_in_parallel(compounds = compounds_with_no_accession,
                                                               URI_EQ = URI_EQ,
                                                               chunk_size = 10,
                                                               max_workers = 4)

    all_valid_compounds = newly_compounds_added + compounds_with_accession

    print("\nNew compounds added to eQuilibrator database")

    print(f"\nSaving all compounds with valid SMILES and accession ids in {outfile_path_for_valid_cpds}")
    with open(outfile_path_for_valid_cpds,'wb') as file:
        pickle.dump(all_valid_compounds, file)