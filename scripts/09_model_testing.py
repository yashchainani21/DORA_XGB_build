import json
from ML_utils import ML_utils
import pandas as pd
import pickle

# ----------------- Start defining fingerprinting parameters -----------------
fp_type = 'ecfp4'
cofactor_configuration = 'by_descending_MW'
max_species = 4
dataset = 'all_BKM_rxns_unreported_means_negative'
query_rules = 'all'  # either 'all' or a list like ['rule0006','rule0007']
model_type = 'XGBoost'
random_state = 42
# ----------------- End defining fingerprinting parameters -----------------

test_fps_filepath = None
input_model_filepath = None
model_test_results_filepath = None

# specify the input filepaths of the reactions and their metadata as well as the output filepaths for fingerprints
if dataset == 'all_AdH_rxns':
    test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_AdH_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

    # specify filepath to load in saved AdH model
    input_model_filepath = f'../models/indiv_models/all_AdH_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

    # assemble a .json filepath to save the performance of the optimized AdH model on the test set
    model_test_results_filepath = f'../models/performance_results/all_AdH_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_test_performance_results.json'

if dataset == 'all_Monox_rxns':
    test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_Monox_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

    # specify filepath to load in saved Monox model
    input_model_filepath = f'../models/indiv_models/all_Monox_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

    # assemble a .json filepath to save the performance of the optimized Monox model on the test set
    model_test_results_filepath = f'../models/performance_results/all_Monox_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_test_performance_results.json'

if dataset == 'all_BKM_rxns':

    if query_rules == 'all':
        test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_BKM_rxns_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

        # specify filepath to load in the complete, consolidated model trained on all BKM reactions
        input_model_filepath = f'../models/consolidated_models/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

        # assemble a .json filepath to save the performance of the optimized, consolidated model on the validation set
        model_test_results_filepath = f'../models/performance_results/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_test_performance_results.json'

    else:

        query_rules_str = '_'.join(query_rules)

        test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/{query_rules_str}_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

        # specify filepath to load in model specific to certain query rules
        input_model_filepath = f'../models/params/{query_rules_str}_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

        # assemble a .json filepath to save the performance of this specific model on the test set
        model_test_results_filepath = f'../models/performance_results/{query_rules_str}_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_test_performance_results.json'

if dataset == 'all_BKM_rxns_unreported_means_negative':

    if query_rules == 'all':
        test_fps_filepath = f'../data/fingerprinted_data/testing_fingerprints/all_BKM_rxns_{fp_type}_test_fingerprints_max_species_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.parquet'

        # specify filepath to load in the complete, consolidated model trained on all BKM reactions
        input_model_filepath = f'../models/consolidated_models/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.pkl'

        # assemble a .json filepath to save the performance of the optimized, consolidated model on the validation set
        model_test_results_filepath = f'../models/performance_results/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_test_performance_results_unreported_means_negative_JN.json'

_, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_configuration,
                                        fp_type = fp_type,
                                        max_species = max_species)

# load in fingerprints and data
test_df = pd.read_parquet(test_fps_filepath)
test_fps = test_df.iloc[:, 0:fp_length].to_numpy()
test_labels = test_df['Label'].to_numpy()

# load in model
with open(input_model_filepath, 'rb') as model_file:
    loaded_model = pickle.load(model_file)

### test model

# start by predicting binary labels on test data
y_predicted_binary_labels = loaded_model.predict(test_fps)

# then predict probabilities on test data
y_predicted_probabilities = loaded_model.predict_proba(test_fps)[:, 1]

# calculate bootstrapped metrics on the test set to evaluate model
mean_auprc, (lower_auprc, upper_auprc) = ML_utils.bootstrap_auprc(y_true = test_labels,
                                                                  y_pred_proba = y_predicted_probabilities
                                                                  , n_iterations = 1000)

mean_precision, (lower_precision, upper_precision) = ML_utils.bootstrap_precision(y_true = test_labels,
                                                                                  y_pred_binary = y_predicted_binary_labels,
                                                                                  n_iterations = 1000)

mean_accuracy, (lower_accuracy, upper_accuracy) = ML_utils.bootstrap_accuracy(y_true = test_labels,
                                                                              y_pred_binary = y_predicted_binary_labels,
                                                                              n_iterations = 1000)

mean_f1, (lower_f1, upper_f1) = ML_utils.bootstrap_f1(y_true = test_labels,
                                                      y_pred_binary = y_predicted_binary_labels,
                                                      n_iterations = 1000)

mean_recall, (lower_recall, upper_recall) = ML_utils.bootstrap_recall(y_true = test_labels,
                                                                      y_pred_binary = y_predicted_binary_labels,
                                                                      n_iterations = 1000)

# store the optimized model's performance on this validation set
results_dict = {'mean AUPRC': mean_auprc,
                'AUPRC lower CI': lower_auprc,
                'AUPRC upper CI': upper_auprc,

                'mean Precision': mean_precision,
                'Precision lower CI': lower_precision,
                'Precision upper CI': upper_precision,

                'mean F1': mean_f1,
                'F1 lower CI': lower_f1,
                'F1 upper CI': upper_f1,

                'mean Recall': mean_recall,
                'Recall lower CI': lower_recall,
                'Recall upper CI': upper_recall,

                'mean Accuracy': mean_accuracy,
                'Accuracy lower CI': lower_accuracy,
                'Accuracy upper CI': upper_accuracy}

with open(model_test_results_filepath, 'w') as file:
    json.dump(results_dict, file, indent=4)
