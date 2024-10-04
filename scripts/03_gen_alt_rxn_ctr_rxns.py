"""
Here, negative enzymatic reactions are synthetically generated using our alternate reaction center hypothesis.
This hypothesis states that for a given reaction center that is transformed during an enzymatic reactions ...,
if other reaction centers exist that do not undergo the same enzymatic transformation ...,
... then the transformation of those alternate reaction centers represents negative reactions that must be infeasible.
This is because the enzyme involved must be specific enough that even when confronted with those reaction centers ...,
.. it does not act upon them.
"""
import pandas as pd
import pickle
import time
import copy
import random
from neg_data_utils import neg_data_utils

### ----- Defining input and output filepaths for various datasets -----

# pick from either 'all_BKM_rxns', 'all_AdH_rxns', 'all_Monox_rxns'
dataset = "all_AdH_rxns"
coreactants_filepath = '../data/raw/all_cofactors.tsv'
gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'
input_filepath = None
output_filepath = None

# define filepaths of reactions labelled by their thermodynamics (i.e. feasibility label assigned with MDF value)
# also define output filepaths once negative reactions have been synthetically generated with pickaxe
if dataset == 'all_AdH_rxns':
    input_filepath = '../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl'
    output_filepath = '../data/processed/all_AdH_rxns_with_synthetic_data.pkl'

if dataset == 'all_Monox_rxns':
    input_filepath = '../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_10.pkl'
    output_filepath = '../data/processed/all_Monox_rxns_with_synthetic_data.pkl'

if dataset == 'all_BKM_rxns':
    input_filepath = '../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl'
    output_filepath = '../data/processed/all_processed_rxns_with_synthetic_data.pkl'

# quick check to make sure we don't end up overwriting any files
assert input_filepath != output_filepath

### ----- preparing the reaction rules for negative data generation -----

# read in the generalized JN1224MIN rule set as a pandas dataframe
rules_df = pd.read_csv('../data/raw/JN1224MIN_rules.tsv', delimiter='\t')

# for each rule, count reactants and products (indicated by 'Any' and strictly excludes cofactors)
num_reactants_per_rule = [ rules_df.loc[i,'Reactants'].split(';').count('Any') for i in range(rules_df.shape[0]) ]
num_products_per_rule = [ rules_df.loc[i,'Products'].split(';').count('Any') for i in range(rules_df.shape[0]) ]

# insert in the number of reactants and number of products per rule into the rules dataframe
rules_df['num_reactants'] = num_reactants_per_rule
rules_df['num_products'] = num_products_per_rule

### ----- first, consider all monosubstrate reactions only -----

# extract all 1 substrate -> 1 product reactions (785 rules)
single_substrate_single_product_df = rules_df[ (rules_df['num_reactants']==1) & (rules_df['num_products']==1) ]
single_substrate_single_product_rules = list(single_substrate_single_product_df['Name'])

# extract all 1 substrate -> 2 product reactions (140 rules)
single_substrate_two_product_df = rules_df[ (rules_df['num_reactants']==1) & (rules_df['num_products']==2) ]
single_substrate_two_product_rules = list(single_substrate_two_product_df['Name'])

# extract all 1 substrate -> 3 product reactions (3 rules)
single_substrate_three_product_df = rules_df[ (rules_df['num_reactants']==1) & (rules_df['num_products']==3) ]
single_substrate_three_product_rules = list(single_substrate_three_product_df['Name'])

# extract all 1 substrate -> 4 product reactions (1 rule)
single_substrate_four_product_df = rules_df[(rules_df['num_reactants'] == 1) & (rules_df['num_products'] == 4)]
single_substrate_four_product_rules = list(single_substrate_four_product_df['Name'])

# ensure that each rule list is unique and that there are no overlaps across lists
assert len(set(single_substrate_single_product_rules)) == len(single_substrate_single_product_rules)
assert len(set(single_substrate_two_product_rules)) == len(single_substrate_two_product_rules)
assert len(set(single_substrate_three_product_rules)) == len(single_substrate_three_product_rules)
assert len(set(single_substrate_four_product_rules)) == len(single_substrate_four_product_rules)

intersections = set(single_substrate_single_product_rules) & set(single_substrate_two_product_rules) & set(single_substrate_three_product_rules) & set(single_substrate_four_product_rules)
assert len(intersections) == 0

### ----- next, consider multisubstrate reactions -----

# extract all multi-substrate reactions, i.e. substrate + substrate + cofactors --> stuff (295 rules)
multi_substrate_df = rules_df[rules_df['num_reactants'] > 1]
multi_substrate_rules = list(multi_substrate_df['Name'])

# now, ensure that all extracted multisubstrate rules are unique
assert len(set(multi_substrate_rules)) == len(list(multi_substrate_rules))

# ensure no overlap across all monosubstrate and multisubstrate rule lists
intersections = set(single_substrate_single_product_rules) & set(single_substrate_two_product_rules) & set(single_substrate_three_product_rules) & set(single_substrate_four_product_rules) & set(multi_substrate_rules)
assert len(intersections) == 0

# finally, check that the total number of rules adds up to 1224
num_total_rules = len(single_substrate_single_product_rules) + len(single_substrate_two_product_rules) + len(single_substrate_three_product_rules) + len(single_substrate_four_product_rules) + len(multi_substrate_rules)
assert num_total_rules == 1224

if __name__ == "__main__":

    # record time when we start synthetically generating negative reaction data
    start_time = time.time()

    # open pickle file from data/processed folder containing LABELLED REACTIONS with reaction thermodynamics
    with open(input_filepath,'rb') as file:
        thermo_labelled_rxns = pickle.load(file)

    # make a copy of this dataset - the copy will also be a list of dictionaries
    rxn_data = copy.deepcopy(thermo_labelled_rxns)

    # count the total number of thermodynamically feasible reactions
    total_thermo_feasible_rxns_count = len([x for x in rxn_data if x['feasibility_label']==1])
    rxns_parsed = 0 # counter to count the number of feasible reactions used for synthetic negative data generation

    print(f"\nThere are {len(thermo_labelled_rxns)} thermodynamically feasible "
          f"reactions out of a total of {len(rxn_data)} reactions")

    # enumerate through positive, thermodynamically feasible reactions to generate infeasible ones
    # the infeasible reactions then get added to rxn_data, which is a copy of the feasible reactions list
    for rxn in thermo_labelled_rxns:

        mapped_rule = rxn['Rule']

        # to synthetically generate negative data, consider only known positive reactions
        # or put differently, only thermodynamically feasible reactions are used to generate infeasible reactions
        if rxn['feasibility_label'] == 1:

            # for reactions mapped onto rules involving only single substrate -> single product reactions,
            # we generate alternate products then ensure the alternate products are distinct from the reported product
            if mapped_rule in single_substrate_single_product_rules:

                neg_data_list = neg_data_utils.run_pickaxe_if_single_reported_product(
                                                                  pos_rxn_entry = rxn,
                                                                  coreactants_filepath = coreactants_filepath,
                                                                  reported_substrate = rxn['Substrates'][0],
                                                                  reported_product = rxn['Products'][0],
                                                                  gen_rules_filepath = gen_rules_filepath)

                # combine list of thermodynamically labelled reaction data with synthetically generated negative data
                rxn_data += neg_data_list

                # print an update on how many reactions have been used for synthetic negative data generation
                rxns_parsed += 1
                print(f"\n{rxns_parsed} reactions out of {total_thermo_feasible_rxns_count}")

            # for reactions mapped onto rules involving a single substrate going to more than or equal to 2 products,
            # we generate alternate products that come in sets of >=2 products as well
            # each set of new products now needs to be checked against the set of reported products
            # this in turn helps to extract Pickaxe-produced reactions that are truly infeasible
            if mapped_rule in single_substrate_two_product_rules + single_substrate_three_product_rules + single_substrate_four_product_rules:

                neg_data_list = neg_data_utils.run_pickaxe_if_multiple_reported_products(
                                                                   pos_rxn_entry = rxn,
                                                                   coreactants_filepath = coreactants_filepath,
                                                                   reported_substrate = rxn['Substrates'][0],
                                                                   reported_products = rxn['Products'],
                                                                   gen_rules_filepath = gen_rules_filepath)

                # combine list of thermodynamically labelled reaction data with synthetically generated negative data
                rxn_data += neg_data_list

                # print an update on how many reactions have been used for synthetic negative data generation
                rxns_parsed += 1
                print(f"\n{rxns_parsed} reactions out of {total_thermo_feasible_rxns_count}")

            # for reactions that involve more than two substrates, we first check if the substrates are identical
            # if they are, then we only need to expand on one of them
            if mapped_rule in multi_substrate_rules:

                # let's check first if the substrates are identical
                identical_substrates = all(substrate == rxn['Substrates'][0] for substrate in rxn['Substrates'])

                # if substrates in the reported reaction are identical, we only need to expand on one of them
                if identical_substrates:
                    neg_data_list = neg_data_utils.run_pickaxe_for_dimerization_rxns(
                                                        pos_rxn_entry = rxn,
                                                        coreactants_filepath = coreactants_filepath,
                                                        reported_substrate = rxn['Substrates'][0], # one substrate
                                                        reported_products = rxn['Products'],
                                                        gen_rules_filepath = gen_rules_filepath)

                    # combine list of thermodynamically labelled reaction data with synthetically generated negative data
                    rxn_data += neg_data_list


                # for reactions involving non-identical substrates, we generate alternate products with all substrates
                # this means that if we have a reported reaction of the form A + B --> X + Y,
                # then we will run a dimerization reaction with Pickaxe of the form A + A with the corresponding rule
                # and another dimerization reaction with Pickaxe of the form B + B with the corresponding rule

                # to extrapolate this further, if we have a reported reaction of the form A + B + C --> X + Y,
                # with some rule requiring any + any + any on the substrates side,
                # such as, rule0251 with any + any + any --> any + any + water
                # then we will run three trimerization reactions with Pickaxe,
                # i.e. we will run A + A + A --> stuff, B + B + B --> stuff, and C + C + C --> stuff
                else:
                    for substrate in rxn['Substrates']:
                        neg_data_list = neg_data_utils.run_pickaxe_for_dimerization_rxns(
                                                        pos_rxn_entry = rxn,
                                                        coreactants_filepath = coreactants_filepath,
                                                        reported_substrate = substrate,
                                                        reported_products = rxn['Products'],
                                                        gen_rules_filepath = gen_rules_filepath)

                        # Once each substrate in this multi-substrate reaction has been expanded upon,
                        # combine the list of thermodynamically labelled reaction data with synthetically generated negative data
                        rxn_data += neg_data_list

                # print an update on how many reactions have been used for synthetic negative data generation
                rxns_parsed += 1
                print(f"\n{rxns_parsed} reactions out of {total_thermo_feasible_rxns_count}")

    # print total number of reactions once -ve reactions have been generated from thermodynamically feasible reactions
    print(f"\nTotal number of positive and negative reactions now: {len(rxn_data)}")

    # print time taken
    end_time = time.time()
    print(f"\nTime taken: {end_time - start_time:.2f} seconds")

    # finally, jumble up the order of reactions in rxn_data so that all negative reactions are not at the end
    random.shuffle(rxn_data)

    # save the new set of reactions
    with open(output_filepath,'wb') as file:
       pickle.dump(rxn_data,file)