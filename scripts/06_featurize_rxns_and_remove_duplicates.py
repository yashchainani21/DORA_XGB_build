"""
Here, we consolidate fingerprint reactions that have previously been processed and labelled by their feasibility.
Reactions that have previously been labelled as feasible correspond to known thermodynamically feasible reactions.
Meanwhile, reactions that have been labelled as infeasible could either be thermodynamically infeasible ...,
... or they could be synthetically generated infeasible reaction because an enzyme was specific enough to not ...
... transform alternate reaction centers in a given enzymatic reactions.
Reactions are fingerprinted here by first creating molecular fingerprints of participating species ...
... and subsequently, assembling these molecular fingerprints together to create reaction fingerprints.
The cofactor configuration method chosen governs how molecular fingerprints are arranged to create reaction fingerprints.
Finally, to prevent any leakage, duplicate fingerprints are removed.
"""
import time
import pandas as pd
import numpy as np
from featurizations import featurizations
from ML_utils import ML_utils

gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'
gen_rules_list = list(pd.read_csv(gen_rules_filepath,delimiter='\t')['Name'])

# ----------------- Start defining fingerprinting parameters -----------------
fp_type = 'ecfp4'
cofactor_configuration = 'by_descending_MW'
max_species = 4
dataset = 'all_BKM_rxns_unreported_means_negative'
query_rules = 'all' # or something like : ['rule0006','rule0007']
batch_num = 1
batch_size = 10000
# ----------------- End defining fingerprinting parameters -----------------

input_filepath = None
output_filepath = None

# read in cofactors list that will be needed to recognize cofactors in a reaction string
# these will help in arranging reaction fingerprints according to the chosen cofactor configuration
all_cofactors_wo_stereo_filepath = '../data/processed/expanded_cofactors_no_stereochem.tsv'
cofactors_df = pd.read_csv(all_cofactors_wo_stereo_filepath, delimiter=',')
all_cofactors_wo_stereo = set(cofactors_df['SMILES'])

# for alcohol dehydrogenase reactions, specify input filepath of reactions and metadata as well as output filepaths for fingerprints
if dataset == 'all_AdH_rxns':

    # when fingerprinting AdH reactions, only rules 0002 and 0003 are allowed
    # this is more of an additional check since the input AdH dataset would only have these rules anyway
    query_rules = ['rule0002','rule0003']

    input_filepath = '../data/processed/00_AdH_rxns_ready_to_fingerprint.parquet'
    output_filepath = f'../data/fingerprinted_data/all_AdH_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

# for monooxygenase reactions, specify input filepath of reactions and metadata as well as output filepaths for fingerprints
if dataset == 'all_Monox_rxns':

    # when fingerprinting Monox reactions, only rules 0004 and 0005 are allowed
    # this is more of an additional check since the input Monox dataset would only have these rules anyway
    query_rules = ['rule0004','rule0005']

    input_filepath = '../data/processed/00_Monox_rxns_ready_to_fingerprint.parquet'
    output_filepath = f'../data/fingerprinted_data/all_Monox_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

# for all BKM reactions, the input rules should've been defined by the user above
# but if we want to fingerprint all rules (i.e. query_rules == 'all'), then we have to do so in batches
if dataset == 'all_BKM_rxns':
    input_filepath = '../data/processed/00_all_rxns_ready_to_fingerprint.parquet' # read in all BKM reactions

    if query_rules == 'all':
        output_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_{batch_num}.parquet'
    else:
        query_rules_str = '_'.join(query_rules)
        output_filepath = f'../data/fingerprinted_data/{query_rules_str}_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

if dataset == 'all_BKM_rxns_unreported_means_negative':
    input_filepath = '../data/processed/00_all_rxns_ready_to_fingerprint_unreported_means_negative_JN.parquet'
    output_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_{batch_num}_unreported_means_negative_JN.parquet'

rxns_df = pd.read_parquet(input_filepath)

# when fingerprinting all BKM reactions, we feed in the entire list of generalized rules so that all reactions get fingerprinted
if query_rules == 'all':
    # try-except here because the unreported is negative dataset doesn't have a 'rule' entry
    try:
        rxns_df = rxns_df[rxns_df['Rule'].isin(gen_rules_list)] # include all generalized rules, hence we use gen_rules_list
    except KeyError:
        pass

    # when featurizing all reactions, we need to do this in batches
    start_index = (batch_num - 1) * 10000 # define start index to start fingerprinting
    end_index = batch_num * 10000 # define end index to stop fingerprinting
    rxns_df = rxns_df.iloc[start_index:end_index, :]

# otherwise, we extract only the reactions mapped to the query rules that we want to fingerprint
else:
    rxns_df = rxns_df[rxns_df['Rule'].isin(query_rules)]

# initialize a fingerprint of zeros depending on the cofactor configuration and fingerprint type chosen
initial_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_configuration,
                                                   fp_type = fp_type,
                                                   max_species = max_species)

# reassign this initial fingerprint to a new variable
# this will form the base upon which all reaction fingerprints will be vertically stacked
all_fps = initial_fp

# initialize trackers to count the number of reactions that could and could not be successfully featurized
num_reactions_featurized = 0
num_reactions_not_featurized = 0

# record the start time that we begin fingerprinting
start_time = time.time()

for i in range(0,rxns_df.shape[0]):
    rxn_str = rxns_df.iloc[i,:]['Reaction eq']
    reaction_object = featurizations.reaction(rxn_str)

    try:
        rxn_fp = reaction_object.rxn_2_fp_w_positioning(fp_type = fp_type,
                                                        cofactor_positioning = cofactor_configuration,
                                                        max_species = max_species,
                                                        all_cofactors_wo_stereo = all_cofactors_wo_stereo)

        # convert any nan values in the fingerprint to zeros before vertical stacking (useful for mordred descriptors)
        rxn_fp = np.nan_to_num(rxn_fp, copy=True)
        all_fps = np.vstack((all_fps, rxn_fp))
        num_reactions_featurized += 1

    except Exception as e:
        num_reactions_not_featurized += 1
        print(e)

    print(f"Number of reactions featurized: {num_reactions_featurized}")
    print(f"Number of reactions not featurized: {num_reactions_not_featurized}")
    print('-------')

# remove the first row, which was the dummy row formed by the previously initialized fingerprint
all_fps = all_fps[1:]

# convert the 2D numpy array of fingerprints into a pandas dataframe
column_names = [ f'feature_{i}' for i in range(0,fp_length) ]
all_fps_df = pd.DataFrame(all_fps, columns = column_names)

# concatenate the dataframe of fingerprints with the previously loaded in dataframe of reaction strings and metadata
all_fps_df.reset_index(drop=True, inplace=True)
rxns_df.reset_index(drop=True, inplace=True)
final_df = pd.concat([all_fps_df, rxns_df],axis=1) # we reset and drop indexes to ensure the dataframes are aligned before their merger

# if not fingerprinting absolutely all BKM reactions,
if query_rules != 'all':

    # then, remove any duplicates on the basis of reaction fingerprints only (keep first occurrence and remove all else)
    final_df_unique = final_df.drop_duplicates( keep = 'first', subset = [ f'feature_{i}' for i in range(0, fp_length) ] )

    print(f"\nNumber of duplicates counted using reaction fingerprints only: {final_df.shape[0] - final_df_unique.shape[0]}")
    print(f'\nNumber of featurized reactions remaining: {final_df_unique.shape[0]}')

    # finally, shuffle the rows of the dataframe then save this unique dataframe
    final_df_unique = final_df_unique.sample(frac = 1, random_state = 42)
    final_df_unique.to_parquet(output_filepath)

    # print the final shape of the dataframe and the time taken for this operation
    end_time = time.time()
    print(f'\nTotal time taken for fingerprinting and removing duplicates: {end_time - start_time:.2f} seconds')
    print(f'\nFinal shape of the fingerprints and reaction metadata dataframe: {final_df_unique.shape}')

# but if fingerprinting all BKM reactions in batches, then we will take care of duplicates and shuffling rows later
else:
    final_df.to_parquet(output_filepath)

    # print the final shape of the dataframe and the time taken for this operation
    end_time = time.time()
    print(f'\nTotal time taken for fingerprinting and removing duplicates: {end_time - start_time:.2f} seconds')
    print(f'\nFinal shape of the fingerprints and reaction metadata dataframe: {final_df.shape}')


