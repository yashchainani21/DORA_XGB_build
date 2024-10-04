import pickle
import time
import json
from ML_utils import ML_utils
import pandas as pd
import dask.dataframe as dd
from imblearn.over_sampling import SMOTE

# ----------------- Start defining fingerprinting parameters -----------------
fp_type = 'ecfp4'
cofactor_configuration = 'by_descending_MW'
max_species = 4
dataset = 'all_BKM_rxns_unreported_means_negative'
query_rules = 'all' # either 'all' or specific rules like ['rule0006','rule0007']
random_state = 42
model_type = 'XGBoost'
model_used_for_tuning = 'XGBoost'
n_jobs = 5
# ----------------- End defining fingerprinting parameters -----------------

train_fps_filepath = None
val_fps_filepath = None
best_params_filepath = None
output_model_filepath = None
model_val_results_filepath = None
feasibility_model = None

# specify the input filepaths of the reactions and their metadata as well as the output filepaths for fingerprints
if dataset == 'all_AdH_rxns':
    train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_AdH_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_AdH_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

    # assemble a .json filepath to save the AdH model's optimized hyperparameters
    best_params_filepath = f'../models/params/all_AdH_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.json'

    # assemble a .pkl filepath to save the final AdH model
    output_model_filepath = f'../models/indiv_models/all_AdH_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

    # assemble a .json filepath to save the performance of the optimized AdH model on the validation set
    model_val_results_filepath = f'../models/performance_results/all_AdH_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_val_performance_results.json'

if dataset == 'all_Monox_rxns':
    train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_Monox_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
    val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_Monox_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

    # assemble a .json filepath to save the Monox model's optimized hyperparameters
    best_params_filepath = f'../models/params/all_Monox_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.json'

    # assemble a .pkl filepath to save the final Monox model
    output_model_filepath = f'../models/indiv_models/all_Monox_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

    # assemble a .json filepath to save the performance of the optimized Monox model on the validation set
    model_val_results_filepath = f'../models/performance_results/all_Monox_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_val_performance_results.json'

if dataset == 'all_BKM_rxns_unreported_means_negative':
    if query_rules == 'all':
        train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_BKM_rxns_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.parquet'
        val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_BKM_rxns_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.parquet'

        # assemble the hyperparameter filepath
        best_params_filepath = f'../models/params/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.json'

        # assemble a .pkl filepath to save the final complete, consolidated model
        output_model_filepath = f'../models/consolidated_models/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_unreported_means_negative_JN.pkl'

        # assemble a .json filepath to save the performance of the optimized, consolidated model on the validation set
        model_val_results_filepath = f'../models/performance_results/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_val_performance_results_unreported_means_negative_JN.json'

if dataset == 'all_BKM_rxns':

    if query_rules == 'all':
        train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/all_BKM_rxns_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
        val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/all_BKM_rxns_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

        # assemble the hyperparameter filepath
        best_params_filepath = f'../models/params/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.json'

        # assemble a .pkl filepath to save the final complete, consolidated model
        output_model_filepath = f'../models/consolidated_models/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

        # assemble a .json filepath to save the performance of the optimized, consolidated model on the validation set
        model_val_results_filepath = f'../models/performance_results/all_BKM_rxns_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_val_performance_results.json'

    else:

        query_rules_str = '_'.join(query_rules)

        train_fps_filepath = f'../data/fingerprinted_data/training_fingerprints/{query_rules_str}_{fp_type}_train_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'
        val_fps_filepath = f'../data/fingerprinted_data/validation_fingerprints/{query_rules_str}_{fp_type}_val_fingerprints_max_species_{max_species}_{cofactor_configuration}.parquet'

        # assemble the hyperparameter filepath
        best_params_filepath = f'../models/params/{query_rules_str}_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.json'

        # assemble a .pkl filepath to save this model specific to certain query rules
        output_model_filepath = f'../models/params/{query_rules_str}_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}.pkl'

        # assemble a .json filepath to save the performance of this specific model on the validation set
        model_val_results_filepath = f'../models/performance_results/{query_rules_str}_{fp_type}_{model_type}_{max_species}_{cofactor_configuration}_val_performance_results.json'

# calculate the length of the fingerprints we need to extract out depending on cofactor config and fp type
_, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_configuration,
                                        fp_type = fp_type,
                                        max_species = max_species)

# if reading in all data spanning all query rules, use dask
if dataset == 'all_BKM_rxns' and query_rules == 'all':

    training_df = dd.read_parquet(train_fps_filepath)

    training_fps = training_df.iloc[:, 0:fp_length].compute()
    training_fps = training_fps.to_numpy()

    training_labels = training_df['Label'].compute()
    training_labels = training_labels.to_numpy()

    val_df = dd.read_parquet(val_fps_filepath)

    val_fps = val_df.iloc[:, 0:fp_length].compute()
    val_fps = val_fps.to_numpy()

    val_labels = val_df['Label'].compute()
    val_labels = val_labels.to_numpy()

elif dataset == 'all_BKM_rxns_unreported_means_negative' and query_rules == 'all':

    training_df = dd.read_parquet(train_fps_filepath)

    training_fps = training_df.iloc[:, 0:fp_length].compute()
    training_fps = training_fps.to_numpy()

    training_labels = training_df['Label'].compute()
    training_labels = training_labels.to_numpy()

    val_df = dd.read_parquet(val_fps_filepath)

    val_fps = val_df.iloc[:, 0:fp_length].compute()
    val_fps = val_fps.to_numpy()

    val_labels = val_df['Label'].compute()
    val_labels = val_labels.to_numpy()

# otherwise, pandas is sufficient for specific rules
else:

    training_df = pd.read_parquet(train_fps_filepath)
    training_fps = training_df.iloc[:,0:fp_length].to_numpy()
    training_labels = training_df['Label'].to_numpy()

    val_df = pd.read_parquet(val_fps_filepath)
    val_fps = val_df.iloc[:,0:fp_length].to_numpy()
    val_labels = val_df['Label'].to_numpy()

# perform SMITE oversampling of the minority class (positive reaction examples)
print('\nPerforming SMOTE oversampling of the minority class')

try:
    smote = SMOTE(sampling_strategy = 0.5, k_neighbors = 10, random_state = random_state)
    training_fps, training_labels = smote.fit_resample(X = training_fps, y = training_labels)
except:
    pass

print('\nFinished SMOTE')

# run hyperparameter optimization
print(f'\nStarting Bayesian hyperparameter optimization. Model selected for tuning: {model_used_for_tuning}')
start_time = time.time()
opt_params = ML_utils.run_bayesian_hyperparameter_search(model_type = model_used_for_tuning,
                                                         X_train = training_fps,
                                                         y_train = training_labels,
                                                         X_val = val_fps,
                                                         y_val = val_labels,
                                                         random_state = random_state,
                                                         init_points = 5,
                                                         n_iter = 20,
                                                         n_jobs = n_jobs)
end_time = time.time()
print(f'\nFinished Bayesian hyperparameter optimization. Time taken: {end_time - start_time:.2f} seconds')

# save the optimized hyperparameters to a json file
with open(best_params_filepath,'w') as json_file:
    json.dump(opt_params, json_file, indent = 4)

# train ML model
if model_type == 'XGBoost':
    feasibility_model = ML_utils.train_XGBoost_model(random_state = random_state,
                                                        n_jobs = n_jobs,
                                                        opt_hyerparams = opt_params,
                                                        model_used_for_tuning = model_used_for_tuning,
                                                        X_train = training_fps,
                                                        y_train = training_labels)

    # Save the trained XGBoost model
    with open(output_model_filepath, 'wb') as model_file:
        pickle.dump(feasibility_model, model_file)

if model_type == 'Logistic':
    feasibility_model = ML_utils.train_logistic_model(random_state = random_state,
                                                      n_jobs = n_jobs,
                                                      opt_hyperparams = opt_params,
                                                      X_train = training_fps,
                                                      y_train = training_labels)

    # save the trained logistic model
    with open(output_model_filepath, 'wb') as model_file:
        pickle.dump(feasibility_model, model_file)

if model_type == 'SVM':
    feasibility_model = ML_utils.train_SVM_model(random_state = random_state,
                                                 opt_hyperparams = opt_params,
                                                 X_train = training_fps,
                                                 y_train = training_labels)

    # save the trained SVM model
    with open(output_model_filepath, 'wb') as model_file:
        pickle.dump(feasibility_model, model_file)

if model_type == 'Random_forest':
    feasibility_model = ML_utils.train_Random_forest_model(random_state = random_state,
                                                           n_jobs=n_jobs,
                                                           opt_hyperparams = opt_params,
                                                           X_train = training_fps,
                                                           y_train = training_labels)

    # save the trained random forest model
    with open(output_model_filepath, 'wb') as model_file:
        pickle.dump(feasibility_model, model_file)

### Evaluate ML model on validation set

# start by predicting binary labels on validation data
y_predicted_binary_labels = feasibility_model.predict(val_fps)

# then predict probabilities on validation data
y_predicted_probabilities = feasibility_model.predict_proba(val_fps)[:, 1]

# calculate bootstrapped metrics on the validation set to evaluate model
mean_auprc, (lower_auprc, upper_auprc) = ML_utils.bootstrap_auprc(y_true = val_labels,
                                                                  y_pred_proba = y_predicted_probabilities
                                                                  , n_iterations = 1000)

mean_precision, (lower_precision, upper_precision) = ML_utils.bootstrap_precision(y_true = val_labels,
                                                                                  y_pred_binary = y_predicted_binary_labels,
                                                                                  n_iterations = 1000)

mean_accuracy, (lower_accuracy, upper_accuracy) = ML_utils.bootstrap_accuracy(y_true = val_labels,
                                                                              y_pred_binary = y_predicted_binary_labels,
                                                                              n_iterations = 1000)

mean_f1, (lower_f1, upper_f1) = ML_utils.bootstrap_f1(y_true = val_labels,
                                                      y_pred_binary = y_predicted_binary_labels,
                                                      n_iterations = 1000)

mean_recall, (lower_recall, upper_recall) = ML_utils.bootstrap_recall(y_true = val_labels,
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

with open(model_val_results_filepath, 'w') as file:
    json.dump(results_dict, file, indent=4)



