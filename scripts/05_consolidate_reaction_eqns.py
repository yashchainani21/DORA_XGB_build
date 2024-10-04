# When we consolidate reaction equations from the master dataset that contains both feasible and infeasible reactions,
# most duplicates exist amongst infeasible reactions only
# but a small number of duplicates exist between feasible and infeasible reactions
# this means that there are some reactions which are feasible and have the exact same reaction string as infeasible reactions
# so as we are removing duplicates, we are going to remove duplicates only from the infeasible set
# and hold onto the reaction string in the feasible set

# Also, note that there are some rule0001 and rule0126 reactions that actually produce water instead of strictly any + any
# so think about how to parse such reaction strings
import pickle
import pandas as pd
from featurizations import featurizations
from preprocessing import preprocessing

### ----- Defining input and output filepaths for various datasets -----

# pick from either 'all_BKM_rxns', 'all_AdH_rxns', 'all_Monox_rxns'
dataset = "all_BKM_rxns_unreported_means_negative"
coreactants_filepath = '../data/raw/all_cofactors.tsv'
gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'
rules_df = pd.read_csv('../data/raw/JN1224MIN_rules.tsv', delimiter='\t')
input_filepath = None
output_filepath = None

all_cofactors_wo_stereo_filepath = '../data/processed/expanded_cofactors_no_stereochem.tsv'
cofactors_df = pd.read_csv(all_cofactors_wo_stereo_filepath, delimiter=',')
all_cofactors_wo_stereo = set(cofactors_df['SMILES'])

# define filepaths of reactions that have both thermodynamically feasible and infeasible reactions
# as well as synthetically generated infeasible reactions (via alt rxn ctr hypothesis and pickaxe)
# also define output filepaths of dataframes to save as parquet files once reactions have been collated
if dataset == 'all_AdH_rxns':
    input_filepath = '../data/processed/all_AdH_rxns_with_synthetic_data.pkl'
    output_filepath = '../data/processed/00_AdH_rxns_ready_to_fingerprint.parquet'

if dataset == 'all_Monox_rxns':
    input_filepath = '../data/processed/all_Monox_rxns_with_synthetic_data.pkl'
    output_filepath = '../data/processed/00_Monox_rxns_ready_to_fingerprint.parquet'

if dataset == 'all_BKM_rxns':
    input_filepath = '../data/processed/all_processed_rxns_with_synthetic_data.pkl'
    output_filepath = '../data/processed/00_all_rxns_ready_to_fingerprint.parquet'

if dataset == 'all_BKM_rxns_unreported_means_negative':
    input_filepath = '../data/processed/all_processed_rxns_with_synthetic_data_unreported_means_negative.pkl'
    output_filepath = '../data/processed/00_all_rxns_ready_to_fingerprint_unreported_means_negative.parquet'

with open(input_filepath,'rb') as file:
    rxns_list = pickle.load(file)

# initialize empty lists to store reaction equations and labels along with other metadata
all_rxn_IDs = []
all_rxn_eqs = []
all_rxn_rules = []
all_rxn_remarks = []
all_rxn_labels = []

# initialize trackers to count the number of reactions that don't fit the description of their mapped rule
# this helps to remove incorrectly formatted reactions and is mostly a problem with synthetically generated reactions
# for instance, a synthetic rule0001 reaction may produce water even though rule0001 reactions should not produce cofactors
# for now, we remove such incorrectly formatted reactions but if our model performs poorly, revisit this decision
num_incorrect_feasible_rxns_count = 0
num_incorrect_infeasible_rxns_count = 0

# initialize trackers to count the number of feasible and infeasible reactions that are correctly formatted
num_correct_feasible_rxns_count = 0
num_correct_infeasible_rxns_count = 0

# track the number of feasible and infeasible reactions parsed overall
num_feasible_rxns_parsed = 0
num_infeasible_rxns_parsed = 0

# parse through each reaction in the chosen dataset but store only the feasible ones first
# a reaction was labelled as feasible if and only if it was thermodynamically feasible by MDF analysis
for rxn in rxns_list:
    if rxn['remark'] == 'pos rxn by thermo':

        # check that this reaction has been labelled as being feasible
        assert rxn['feasibility_label'] == 1

        # extract the rule that this reaction has been mapped onto then get the format of the rule
        mapped_rule = rxn['Rule']
        num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df, mapped_rule)

        # now check that the format of the rule matches the substrates, products and cofactors from the reaction string
        # substrates, products, and cofactors can be extracted from a reaction string via our custom featurizations package
        try:
            # if the reaction string is correctly formatted ...
            rxn_object = featurizations.reaction(rxn_str = rxn['Reaction eq'])
            assert num_substrates == len(rxn_object.get_substrates(all_cofactors_wo_stereo))
            assert num_products == len(rxn_object.get_products(all_cofactors_wo_stereo))
            assert len(lhs_cofactors) == len(rxn_object.get_lhs_cofactors(all_cofactors_wo_stereo))
            assert len(rhs_cofactors) == len(rxn_object.get_rhs_cofactors(all_cofactors_wo_stereo))

            # ... but it has already been stored, don't store it a second time
            if rxn['Reaction eq'] not in all_rxn_eqs:
                all_rxn_IDs.append(rxn['Reaction ID'])
                all_rxn_eqs.append(rxn['Reaction eq'])
                all_rxn_labels.append(rxn['feasibility_label'])
                all_rxn_remarks.append(rxn['remark'])
                all_rxn_rules.append(rxn['Rule'])

                num_correct_feasible_rxns_count += 1

        # and if this reaction string is not correctly formatted ...
        except AssertionError:
            # we track the number of feasible reactions with this incorrect formatting
            num_incorrect_feasible_rxns_count += 1

    # update total number of feasible reactions parsed
    num_feasible_rxns_parsed += 1
    print(f"\n----Number of total reactions parsed: {num_feasible_rxns_parsed}")

print("\n----All feasible reactions parsed----")

# now, parse through reactions again but storing only infeasible ones this time
# a reaction is labelled infeasible for two reasons - either the reported reaction was thermodynamically infeasible
# or an infeasible reaction was synthetically generated from a feasible reaction
for rxn in rxns_list:
    if rxn['remark'] == 'neg rxn by thermo' or rxn['remark'] == 'alternate product from pos rxn':

        # check that this reaction has been labelled as being infeasible
        assert rxn['feasibility_label'] == 0

        # extract the rule that this reaction has been mapped onto to get the format of the rule
        try:
            mapped_rule = rxn['Rule']

        except KeyError: # different key if a negative reaction was synthetically generated
            mapped_rule = rxn['Original rule']

        num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df, mapped_rule)

        # now check that the format of the rule matches that of the reaction string
        try:
            rxn_object = featurizations.reaction(rxn_str=rxn['Reaction eq'])
            assert num_substrates == len(rxn_object.get_substrates(all_cofactors_wo_stereo))
            assert num_products == len(rxn_object.get_products(all_cofactors_wo_stereo))
            assert len(lhs_cofactors) == len(rxn_object.get_lhs_cofactors(all_cofactors_wo_stereo))
            assert len(rhs_cofactors) == len(rxn_object.get_rhs_cofactors(all_cofactors_wo_stereo))

            # if the reaction equation has already been stored in the feasible reactions loop, don't store this reaction
            # i.e. if duplicates exist between a feasible and an infeasible reaction, then only the feasible one is stored
            # this is because the number of feasible reactions dwarf that of infeasible reactions anyway
            if rxn['Reaction eq'] not in all_rxn_eqs:
                try:
                    all_rxn_IDs.append(rxn['Reaction ID'])
                except KeyError:
                    all_rxn_IDs.append(rxn['Original reaction ID'])

                try:
                    all_rxn_rules.append(rxn['Rule'])
                except KeyError:
                    all_rxn_rules.append(rxn['Original rule'])

                all_rxn_eqs.append(rxn['Reaction eq'])
                all_rxn_labels.append(rxn['feasibility_label'])
                all_rxn_remarks.append(rxn['remark'])

                num_correct_infeasible_rxns_count += 1

        except AssertionError:
            num_incorrect_infeasible_rxns_count += 1
            pass

    num_infeasible_rxns_parsed += 1
    print(f"\n----Number of total reactions parsed: {num_infeasible_rxns_parsed}")

# ensure all the stored reactions are unique
assert len(set(all_rxn_eqs)) == len(all_rxn_eqs)
print(len(set(all_rxn_eqs)))

print(f"\nNumber of feasible reactions incorrectly formatted and had to be discarded: {num_incorrect_feasible_rxns_count}")
print(f"\nNumber of infeasible reactions incorrectly formatted and had to be discarded: {num_incorrect_infeasible_rxns_count}")
print(f"\nNumber of unique feasible reactions stored: {num_correct_feasible_rxns_count}")
print(f"\nNumber of unique infeasible reactions stored: {num_correct_infeasible_rxns_count}")

print(f"\nNumber of duplicate and incorrect reactions removed: {len(rxns_list) - len(all_rxn_eqs)}")

df = pd.DataFrame({ 'Reaction eq': all_rxn_eqs,
                    'Label': all_rxn_labels,
                    'Reaction ID': all_rxn_IDs,
                    'Rule': all_rxn_rules,
                    'Remarks': all_rxn_remarks})

df.to_parquet(output_filepath)