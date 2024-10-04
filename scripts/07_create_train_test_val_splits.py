"""
The goal of this script is to create training, testing, and validation sets. Ideally, the distribution
of feasible/ infeasible reactions within each of these sets should mirror the distribution of
feasible/ infeasible reactions in our entire dataset. In order to ensure this, we perform stratified splits.
Additionally, we also want each reaction rule to be represented in the train, test, and validation sets.
This in turn ensures that all types of reactions are seen by the model while preventing the model from
overfitting to the most common reaction type and underperforming on rarer types.

Here, we perform train/ test/ val splits of 80/ 10/ 10 on input reactions on the level of reaction rules for all
generalized reaction rules that have at least 10 examples mapped to them. For generalized reaction rules
that have less than 10 examples mapped to them, these examples are all lumped together and an 80/ 10/ 10
train/ test/ val split is then performed on these lumped rare reactions. This approach maximizes the number
of reaction rules that make it to each set.
"""
import pandas as pd
from sklearn.model_selection import train_test_split
from ML_utils import ML_utils
import dask.dataframe as dd

# ----------------- Start defining fingerprinting parameters -----------------
fp_type = 'ecfp4'
cofactor_configuration = 'by_descending_MW'
max_species = 4
dataset = 'all_BKM_rxns_unreported_means_negative'
query_rules = 'all' # either 'all' or a list of rules, e.g. ['rule0006','rule0007']
random_state = 42
# ----------------- End defining fingerprinting parameters -----------------

# calculate the length of the reaction fingerprint given the fingerprint type and cofactor configuration
_, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_configuration,
                                          fp_type = fp_type,
                                          max_species = max_species)

# initialize None for various filepaths before we define them further
input_filepath = None
train_fps_filepath = None
test_fps_filepath = None
val_fps_filepath = None

batch1_input_filepath = None
batch2_input_filepath = None
batch3_input_filepath = None
batch4_input_filepath = None
batch5_input_filepath = None
batch6_input_filepath = None
batch7_input_filepath = None
batch8_input_filepath = None
batch9_input_filepath = None
batch10_input_filepath = None
batch11_input_filepath = None

# specify the input filepaths of the reactions and their metadata as well as the output filepaths for fingerprints
if dataset == 'all_AdH_rxns':

    # when fingerprinting AdH reactions, only rules 0002 and 0003 are allowed
    # this is more of an additional check since the input AdH dataset would only have these rules anyway
    query_rules = ['rule0002','rule0003']

    input_filepath = f'../data/fingerprinted_data/all_AdH_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_AdH_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_AdH_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_AdH_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

if dataset == 'all_Monox_rxns':

    # when fingerprinting Monox reactions, only rules 0004 and 0005 are allowed
    # this is more of an additional check since the input Monox dataset would only have these rules anyway
    query_rules = ['rule0004','rule0005']

    input_filepath = f'../data/fingerprinted_data/all_Monox_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

    train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_Monox_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_Monox_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_Monox_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

if dataset == 'all_BKM_rxns': # for all BKM reactions, the input rules should've been defined by the user above

    if query_rules == 'all':

        # ---- define input filepaths when reading in all batches of absolutely all BKM reactions ----
        batch1_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_1.parquet'
        batch2_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_2.parquet'
        batch3_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_3.parquet'
        batch4_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_4.parquet'
        batch5_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_5.parquet'
        batch6_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_6.parquet'
        batch7_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_7.parquet'
        batch8_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_8.parquet'
        batch9_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_9.parquet'
        batch10_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_10.parquet'
        batch11_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_11.parquet'
        batch12_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_12.parquet'
        batch13_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_13.parquet'
        batch14_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_14.parquet'

        # ---- define output filepaths for the consolidated BKM reactions ----
        train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_BKM_rxns_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
        test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_BKM_rxns_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
        val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_BKM_rxns_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    else:
        query_rules_str = '_'.join(query_rules)

        input_filepath = f'../data/fingerprinted_data/{query_rules_str}_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

        train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/{query_rules_str}_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
        test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/{query_rules_str}_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
        val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/{query_rules_str}_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

if dataset == 'all_BKM_rxns_unreported_means_negative':

    if query_rules == 'all':
        # ---- define input filepaths when reading in all batches of absolutely all BKM reactions ----
        batch1_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_1_unreported_means_negative_JN.parquet'
        batch2_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_2_unreported_means_negative_JN.parquet'
        batch3_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_3_unreported_means_negative_JN.parquet'
        batch4_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_4_unreported_means_negative_JN.parquet'
        batch5_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_5_unreported_means_negative_JN.parquet'
        batch6_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_6_unreported_means_negative_JN.parquet'
        batch7_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_7_unreported_means_negative_JN.parquet'
        batch8_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_8_unreported_means_negative_JN.parquet'
        batch9_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_9_unreported_means_negative_JN.parquet'
        batch10_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_10_unreported_means_negative_JN.parquet'
        batch11_input_filepath = f'../data/fingerprinted_data/all_BKM_rxns_{fp_type}_fingerprints_max_species_{max_species}_{cofactor_configuration}_batch_11_unreported_means_negative_JN.parquet'

        # ---- define output filepaths for the consolidated BKM reactions ----
        train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_BKM_rxns_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.parquet'
        test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_BKM_rxns_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.parquet'
        val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_BKM_rxns_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.parquet'

# if specific query rules have been requested rather than absolutely all 1224 generalized reaction rules,
if query_rules != 'all':

    # then we read in the corresponding input fingerprints
    rxns_with_fps_df = pd.read_parquet(input_filepath)

    # and we extract out only reactions mapped to the query rules specified above
    rxns_with_fps_df = rxns_with_fps_df[rxns_with_fps_df['Rule'].isin(query_rules)]

# but if fingerprinting absolutely all BKM reactions
else:
    # first, we need to combine the batched fingerprints into a single dataframe
    # here, given the size of the data, we use dask rather than pandas
    rxns_with_fps_df_batch_1 = dd.read_parquet(batch1_input_filepath)
    rxns_with_fps_df_batch_2 = dd.read_parquet(batch2_input_filepath)
    rxns_with_fps_df_batch_3 = dd.read_parquet(batch3_input_filepath)
    rxns_with_fps_df_batch_4 = dd.read_parquet(batch4_input_filepath)
    rxns_with_fps_df_batch_5 = dd.read_parquet(batch5_input_filepath)
    rxns_with_fps_df_batch_6 = dd.read_parquet(batch6_input_filepath)
    rxns_with_fps_df_batch_7 = dd.read_parquet(batch7_input_filepath)
    rxns_with_fps_df_batch_8 = dd.read_parquet(batch8_input_filepath)
    rxns_with_fps_df_batch_9 = dd.read_parquet(batch9_input_filepath)
    rxns_with_fps_df_batch_10 = dd.read_parquet(batch10_input_filepath)
    rxns_with_fps_df_batch_11 = dd.read_parquet(batch11_input_filepath)

    # vertically stack the dask dataframes
    combined_fps_df = dd.concat(
        [rxns_with_fps_df_batch_1,
         rxns_with_fps_df_batch_2,
         rxns_with_fps_df_batch_3,
         rxns_with_fps_df_batch_4,
         rxns_with_fps_df_batch_5,
         rxns_with_fps_df_batch_6,
         rxns_with_fps_df_batch_7,
         rxns_with_fps_df_batch_8,
         rxns_with_fps_df_batch_9,
         rxns_with_fps_df_batch_10,
         rxns_with_fps_df_batch_11],
         axis = 0,  # Stack vertically (row-wise concatenation)
         ignore_index = True)  # Reset the index in the new DataFrame


    # for our combined dask dataframe, .drop_duplicates will not be executed immediately and we need to add .compute()
    rxns_with_fps_df_unique = combined_fps_df.drop_duplicates(keep='first',
                                                              subset= ( [f'feature_{i}' for i in range(0, fp_length)] + ['Label'] )).compute()

    total_rows = combined_fps_df.shape[0].compute()
    unique_rows = rxns_with_fps_df_unique.shape[0]
    print(f"\nNumber of duplicates counted using reaction fingerprints only: {total_rows - unique_rows}")
    print(f'\nNumber of featurized reactions remaining: {unique_rows}')

    rxns_with_fps_df_unique.dropna(subset='Label', inplace=True)

    print(f'\nNumber of unique rows left after dropping NaN labels: {rxns_with_fps_df_unique.shape[0]}')

    # finally, shuffle the rows of the dataframe then save this unique dataframe
    rxns_with_fps_df = rxns_with_fps_df_unique.sample(frac=1, random_state=42)

# initialize empty DataFrames for train, validation, and test sets
train_df = pd.DataFrame(columns = rxns_with_fps_df.columns)
val_df = pd.DataFrame(columns = rxns_with_fps_df.columns)
test_df = pd.DataFrame(columns = rxns_with_fps_df.columns)

MIN_SAMPLES_PER_RULE = 10

# initialize another empty dataframe for collecting reactions with less than the minimum number of samples per rule
rare_data_df = pd.DataFrame(columns = rxns_with_fps_df.columns)

# for each unique reaction rule in this reactions dataframe,
try:
    for rule in rxns_with_fps_df['Rule'].unique():

        # extract the subset of the reactions dataframe containing reactions mapped to this specific rule only
        subset_df = rxns_with_fps_df[rxns_with_fps_df['Rule'] == rule]

        # if the number of reactions mapped to this rule is less than the user defined minimum,
        if subset_df.shape[0] < MIN_SAMPLES_PER_RULE:
            # then we lump this subset under a dataframe holding all such rare reactions
            rare_data_df = pd.concat([rare_data_df, subset_df])

        # but if the number of reactions mapped to this rule is more than the user defined minimum
        else:
            # then we perform a stratified split into train, test, and val sets within this particular reaction rule
            try:
                # perform stratified split for test (10% of total data)
                s_train_and_val, s_test = train_test_split(subset_df,
                                                           test_size = 0.1,
                                                           stratify = subset_df['Label'],
                                                           random_state = random_state)

                # with the remaining data, perform a stratified split to get validation (10%) and training (80%) sets
                s_train, s_val = train_test_split(s_train_and_val,
                                                  test_size = 1/9,
                                                  stratify = s_train_and_val['Label'],
                                                  random_state = random_state)

            # if we can't perform a stratified split then we perform a regular train/test split
            except ValueError:
                # perform regular split for test (10% of total data)
                s_train_and_val, s_test = train_test_split(subset_df,
                                                           test_size = 0.1,
                                                           random_state = random_state)

                # with the remaining data, perform a stratified split to get validation (10%) and training (80%) sets
                s_train, s_val = train_test_split(s_train_and_val,
                                                  test_size = 1/9,
                                                  random_state=random_state)

            train_df = pd.concat([train_df, s_train])
            val_df = pd.concat([val_df, s_val])
            test_df = pd.concat([test_df, s_test])

except KeyError:
    # perform stratified split for test (10% of total data)
    s_train_and_val, s_test = train_test_split(rxns_with_fps_df,
                                               test_size = 0.1,
                                               stratify = rxns_with_fps_df['Label'],
                                               random_state = random_state)

    # with the remaining data, perform a stratified split to get validation (10%) and training (80%) sets
    s_train, s_val = train_test_split(s_train_and_val,
                                      test_size = 1 / 9,
                                      stratify = s_train_and_val['Label'],
                                      random_state = random_state)

    print(f"\nTrain size: {s_train.shape}")
    print(f"\nValidation size: {s_val.shape}")
    print(f"\nTest size: {s_test.shape}")

    pd.DataFrame(s_train).to_parquet(train_fps_filepath)
    pd.DataFrame(s_val).to_parquet(val_fps_filepath)
    pd.DataFrame(s_test).to_parquet(test_fps_filepath)

    exit()

# now, we handle the aggregated rare reaction rules
if not rare_data_df.empty:
    try:
        # perform stratified split for test (10% of total data)
        r_train_and_val, r_test = train_test_split(rare_data_df,
                                                   test_size = 0.1,
                                                   stratify = rare_data_df['Label'],
                                                   random_state = random_state)

        # with the remaining data, perform a stratified split to get validation (10%) and training (80%) sets
        r_train, r_val = train_test_split(r_train_and_val,
                                          test_size = 1/9,
                                          stratify = r_train_and_val['Label'],
                                          random_state = random_state)

    # if we can't perform a stratified split then we perform a regular train/test split
    except ValueError:
        r_train_and_val, r_test = train_test_split(rare_data_df,
                                                   test_size = 0.1,
                                                   random_state=random_state)
        r_train, r_val = train_test_split(r_train_and_val,
                                          test_size = 1/9,
                                          random_state=random_state)

    train_df = pd.concat([train_df, r_train])
    val_df = pd.concat([val_df, r_val])
    test_df = pd.concat([test_df, r_test])

print(f"\nTrain size: {train_df.shape}")
print(f"\nValidation size: {val_df.shape}")
print(f"\nTest size: {test_df.shape}")

train_df.to_parquet(train_fps_filepath)
val_df.to_parquet(val_fps_filepath)
test_df.to_parquet(test_fps_filepath)
