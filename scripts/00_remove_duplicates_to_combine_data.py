"""
In this script, processed reaction entries from BRENDA, KEGG, and MetaCyc are read into memory and then combined.
To remove duplicate entries, the list of substrates, products, and cofactors as well as reaction equations are used.
Any two entries that have identical reaction equations (as formatted in 00_raw_data_preprocessing.py) as well as ...
... identical lists of substrates, products, and cofactors are counted as duplicate entries and only one entry is held.
"""
import pickle

# enter filepaths for input data
from preprocessing import preprocessing

metacyc_data_path = '../data/interim/metacyc_all_mined_rxns.pkl'
kegg_data_path = '../data/interim/kegg_all_mined_rxns.pkl'
brenda_data_path = '../data/interim/brenda_all_mined_rxns.pkl'

# enter outfile path
outfile_path = '../data/processed/all_processed_rxns_no_thermo.pkl'

with open(metacyc_data_path, 'rb') as file:
    metacyc_data = pickle.load(file)

with open(kegg_data_path, 'rb') as file:
    kegg_data = pickle.load(file)

with open(brenda_data_path, 'rb') as file:
    brenda_data = pickle.load(file)

# merge datasets
combined_data = metacyc_data + kegg_data + brenda_data

print(f'\nBefore removing duplicates, the number of unique reactions left are: {len(combined_data)}')

# remove duplicates by reaction strings as well as substrate, product, and cofactor lists
unique_data = [] # initialize a list to store only the unique dictionaries
seen_combinations = set() # store unique combinations of values

for entry in combined_data:
    identifier = (entry['Reaction eq'],
        tuple(entry['Substrates']),  # Convert list to tuple for hashability
        tuple(entry['Products']),
        tuple(entry['LHS_cofactors']),
        tuple(entry['RHS_cofactors']))

    # check if this combination has already been seen
    if identifier not in seen_combinations:
        unique_data.append(entry) # add unique dictionary to a new list
        seen_combinations.add(identifier) # remember this combination by its identifier

print(f'\nAfter removing duplicates, the number of unique reactions left are: {len(unique_data)}')

print(f'\nAttempting to remove reactions where at least one compound SMILES contains an asterisk')
num_asterisked_rxns = preprocessing.count_processed_rxns_with_asterisked_cpd_smiles(unique_data)
print(f'\n{num_asterisked_rxns} such reactions found')

cleaned_data = preprocessing.remove_processed_rxns_with_asterisked_cpd_smiles(unique_data)
print(f'\nAfter removing reactions, {len(cleaned_data)} cleaned, unique reactions are remaining')

with open(outfile_path,'wb') as file:
    pickle.dump(cleaned_data, file)
