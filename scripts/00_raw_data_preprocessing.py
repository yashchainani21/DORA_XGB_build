"""
Here, raw reaction data from BRENDA, KEGG, and MetaCyc databases used in developing JN1224MIN is processed.
For each reaction entry in each database, compound SMILES are canonicalized and stereochemistry is removed.
Reaction strings of the form "A + B = C + D" are then created for each reaction entry with these sanitized SMILES.
In developing JN1224MIN, each reaction entry in each database has been attempted to be mapped to a reaction rule.
Consequently, only mapped reactions are processed while unmapped ones are not.
Even if a reaction has been mapped to a rule but uses cofactors outside our predetermined list, it is ignored.
"""
import pickle
import pandas as pd
from preprocessing import preprocessing

db = "metacyc_V24"  # choose from "kegg", "brenda", "metacyc", or "metacyc_V24"
filepath_to_save_processed_data = f"../data/interim/{db}_all_mined_rxns.pkl"

# reaction IDs of unmapped reactions and reactions with "exotic cofactors" will be stored here
filepath_to_save_unprocessed_ids = f"../data/interim/{db}_all_mined_unprocessed_rxns.pkl"

raw_data = preprocessing.load_dataset(db)
rxn_rules_filepath = '../data/processed/JN1224MIN_rules.tsv'
cofactors_filepath = '../data/processed/expanded_cofactors_no_stereochem.tsv'
cofactors_df, cofactors_list = preprocessing.load_cofactors(cofactors_filepath)
rules_df = pd.read_csv(rxn_rules_filepath, delimiter='\t')

processed_rxns_count = 0
processed_rxn_ids_list = []
mapped_rules_of_processed_rxns = []
processed_data_storage = []

unprocessed_rxns_count = 0
unprocessed_rxn_ids_list = []
mapped_rules_of_unprocessed_rxns = []

# iterate through each reaction in the chosen database
for i, rxn_id in enumerate(raw_data.keys()):
    rxn_entry = raw_data[rxn_id]

    try: # MetaCyc and BRENDA were mined with UNIPROT IDs (No UNIPROT IDs or enzyme info is used in training)
        reactants_dict, products_dict, UNIPROT_IDs, mapped_rule = rxn_entry[0], rxn_entry[1], rxn_entry[2], rxn_entry[3]

    except IndexError: # KEGG was not mined with UNIPROT IDs
        reactants_dict, products_dict, mapped_rule = rxn_entry[0], rxn_entry[1], rxn_entry[2]

    # Only mined reactions that have been mapped to at least 1 rule in JN1224MIN are considered
    if mapped_rule != 'Unmapped':

        # ------- Mined reaction details -------
        # parse this mined reaction to extract SMILES lists of substrates, products, LHS cofactors, and RHS cofactors
        # this operation will also remove all sterochemistry and canonicalize all SMILES strings
        lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict, cofactors_list)
        rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict, cofactors_list)

        # reorder these SMILES lists in terms of either ascending or descending MW (either is fine but be consistent)
        # a ValueError will be raised here by 'reorder_smiles_by_mol_weight' if stereochemistry was not removed earlier
        # this reordering will help with removing duplicates later
        substrates_rxn = preprocessing.reorder_smiles_by_mol_weight(substrates_rxn, descending = False)
        products_rxn = preprocessing.reorder_smiles_by_mol_weight(products_rxn, descending = False)
        lhs_cofactors_rxn = preprocessing.reorder_smiles_by_mol_weight(lhs_cofactors_rxn, descending = False)
        rhs_cofactors_rxn = preprocessing.reorder_smiles_by_mol_weight(rhs_cofactors_rxn, descending = False)

        # ------- Mapped rule details -------
        # get number of substrates & products and cofactor IDs on LHS & RHS as stipulated by the mapped reaction rule
        num_substrates_rule, num_products_rule, \
        lhs_cofactor_ids_rule, rhs_cofactor_ids_rule = preprocessing.get_rule_info(rules_df = rules_df, query_rule = mapped_rule)

        # swap out lists of cofactor IDs for lists of cofactor SMILES
        lhs_cofactors_rule = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = lhs_cofactor_ids_rule,
                                                             cofactors_df = cofactors_df)
        rhs_cofactors_rule = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = rhs_cofactor_ids_rule,
                                                             cofactors_df = cofactors_df)

        # ------- Checks -------
        # number of substrates and products in mined reaction must equal that stipulated by the reaction rule
        condition_A = ( len(substrates_rxn) == num_substrates_rule )
        condition_B = ( len(products_rxn) == num_products_rule )

        # number of cofactors in mined reaction must equal the number of cofactor IDs (note IDs!) stipulated by rule
        condition_C = ( len(lhs_cofactors_rxn) == len(lhs_cofactor_ids_rule) )
        condition_D = ( len(rhs_cofactors_rxn) == len(rhs_cofactor_ids_rule) )

        # SMILES of cofactors in mined reaction must be in list of cofactor SMILES created from the rule's cofactor ids
        condition_E = all(cofactor in lhs_cofactors_rule for cofactor in lhs_cofactors_rxn)
        condition_F = all(cofactor in rhs_cofactors_rule for cofactor in rhs_cofactors_rxn)

        # for this mined reaction mapped to a reaction rule, if all conditions are satisfied, then save it for training
        if all([condition_A, condition_B, condition_C, condition_D, condition_E, condition_F]):
            print(f"\n{rxn_id} will be used for training")
            processed_rxns_count += 1

            # construct a reaction string using the substrates, products, and cofactors extracted from mined reactions
            rxn_eq = preprocessing.construct_rxn_str(substrates = substrates_rxn,
                                                    lhs_cofactors = lhs_cofactors_rxn,
                                                    products = products_rxn,
                                                    rhs_cofactors = rhs_cofactors_rxn)

            processed_entry = {"Enzymatic database": db,
                               "Reaction ID": rxn_id,
                               "Substrates": substrates_rxn,
                               "LHS_cofactors": lhs_cofactors_rxn,
                               "Products": products_rxn,
                               "RHS_cofactors": rhs_cofactors_rxn,
                               "Reaction eq": rxn_eq,
                               "Rule": mapped_rule}

            processed_data_storage.append(processed_entry)

        # if any of the above conditions are not satisfied, this reaction will not be used for training
        else:
            print(f"{rxn_id} will not be used for training")
            unprocessed_rxn_ids_list.append(rxn_id)
            mapped_rules_of_unprocessed_rxns.append(mapped_rule)
            unprocessed_rxns_count += 1

print('\n---------')
print(f'\nNumber of reactions processed: {processed_rxns_count}')
print(f'\nNumber of reactions not processed: {unprocessed_rxns_count}')

print(f'\nSaving processed data as a pickle file')
with open(filepath_to_save_processed_data, 'wb') as file:
    pickle.dump(processed_data_storage, file)

print(f'\nSaving unprocessed reaction IDs as a pickle file')
with open(filepath_to_save_unprocessed_ids, 'wb') as file:
    pickle.dump(unprocessed_rxn_ids_list, file)