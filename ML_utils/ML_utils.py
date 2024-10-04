import numpy as np
from typing import Tuple, Optional
from bayes_opt import BayesianOptimization
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, average_precision_score

try:
    import lightgbm as lgb
except:
    pass

def initialize_rxn_fp(cofactor_positioning_method: str, fp_type: str,
                      max_species: int) -> Optional[Tuple[np.ndarray, int]]:
    """
    Construct the initial dummy reaction fingerprint based on the cofactor positioning method chosen
    This dummy fingerprint will form the base on which all other reaction fingerprints will be stacked upon
    :param cofactor_positioning_method:
    :param ecfp_bits:
    :param max_species:
    :return:
    """

    # The number of bits of a feature vector depending on the type of fingerprinting method chosen
    nBits = None

    if fp_type == "morgan" or fp_type == "ecfp4":
        nBits = 2048

    if fp_type == "atom_pair":
        nBits = 2048

    if fp_type == "MACCS":
        nBits = 167

    if fp_type == "mordred":
        nBits = 1613

    if fp_type == "MAP4":
        nBits = 2048

    if fp_type == "MHFP4":
        nBits = 2048

    # for 'add_concat', all fingerprints on LHS are added and concatenated with the sum of all fingerprints on RHS
    if cofactor_positioning_method == "add_concat":
        initialized_fp = np.zeros(nBits * 2)
        fp_length = nBits * 2
        return initialized_fp, fp_length

    # for 'add_subtract', all fingerprints on LHS are added and then subtracted from the sum of all fingerprints on LHS
    if cofactor_positioning_method == "add_subtract":
        initialized_fp = np.zeros(nBits)
        fp_length = nBits
        return initialized_fp, fp_length

    # if arranging species at random, all fingerprinting positions are randomized but the number of species and padding is constant
    if cofactor_positioning_method == "half_random":
        initialized_fp = np.zeros(nBits * max_species * 2)
        fp_length = nBits * max_species * 2
        return initialized_fp, fp_length

    if cofactor_positioning_method == "full_random":
        initialized_fp = np.zeros(nBits * max_species * 2)
        fp_length = nBits * max_species * 2
        return initialized_fp, fp_length

    # if arranging species by molecular weight, the number of species and padding is constant
    if cofactor_positioning_method == "by_ascending_MW":
        initialized_fp = np.zeros(nBits * max_species * 2)
        fp_length = nBits * max_species * 2
        return initialized_fp, fp_length

    if cofactor_positioning_method == "by_descending_MW":
        initialized_fp = np.zeros(nBits * max_species * 2)
        fp_length = nBits * max_species * 2
        return initialized_fp, fp_length

    else:
        return None

def XGBC_objective(X_train: np.ndarray,
                   y_train: np.ndarray,
                   X_val: np.ndarray,
                   y_val: np.ndarray,
                   n_jobs: int = 1):
    """
    Objective function for XGBoost hyperparameter optimization via a Bayesian algorithm. This function will be
    passed to an instantiated BayesianOptimization object from the BayesOpt package
    """
    def objective(learning_rate: float, max_leaves: int, max_depth: int, reg_alpha: float, reg_lambda: float,
                   n_estimators: int, min_child_weight: float, colsample_bytree: float, colsample_bylevel: float,
                   colsample_bynode: float, subsample: float, scale_pos_weight: float) -> float:

        params = {'learning_rate': learning_rate,
                  'max_leaves': int(max_leaves),
                  'max_depth': int(max_depth),
                  'reg_alpha': reg_alpha,
                  'reg_lambda': reg_lambda,
                  'n_estimators': int(n_estimators),
                  'min_child_weight': min_child_weight,
                  'colsample_bytree': colsample_bytree,
                  'colsample_bylevel': colsample_bylevel,
                  'colsample_bynode': colsample_bynode,
                  'subsample': subsample,
                  'scale_pos_weight': scale_pos_weight,
                  'objective': 'binary:logistic',
                  'eval_metric': 'logloss',
                  'n_jobs': n_jobs}

        # train XGBoost classifier on training data
        model = XGBClassifier(**params)
        model.fit(X_train, y_train)

        # then evaluate on validation data by predicting probabilities using validation fingerprints
        y_val_predicted_probabilities = model.predict_proba(X_val)[:, 1]

        # finally, calculate the AUPRC score between the validation labels and the validation predicted probabilities
        auprc = average_precision_score(y_val, y_val_predicted_probabilities)
        return auprc

    return objective

def LGBC_objective(X_train: np.ndarray,
                   y_train: np.ndarray,
                   X_val: np.ndarray,
                   y_val: np.ndarray,
                   n_jobs: int = 1):

    def objective(num_leaves: int, learning_rate: float, max_depth: int, reg_alpha: float,
                   reg_lambda: float, n_estimators: int, min_child_weight, colsample_bytree: float,
                   subsample: float, scale_pos_weight: float):
        params = {
            'num_leaves': int(num_leaves),
            'learning_rate': learning_rate,
            'max_depth': int(max_depth),
            'reg_alpha': reg_alpha,
            'reg_lambda': reg_lambda,
            'n_estimators': int(n_estimators),
            'min_child_weight': int(min_child_weight),
            'colsample_bytree': colsample_bytree,
            'subsample': subsample,
            'scale_pos_weight': scale_pos_weight,
            'objective': 'binary',
        }
        model = lgb.LGBMClassifier(n_jobs=n_jobs, **params)
        model.fit(X_train, y_train)
        y_val_predicted_probabilities = model.predict_proba(X_val)[:, 1]
        auprc = average_precision_score(y_val, y_val_predicted_probabilities)
        return auprc
    return objective

def logistic_regression_objective(X_train: np.ndarray,
                                  y_train: np.ndarray,
                                  X_val: np.ndarray,
                                  y_val: np.ndarray,
                                  n_jobs: int = 1):
    def objective(C):
        # Logistic Regression expects penalty as a string and C as a float
        model = LogisticRegression(C = float(C), solver = 'liblinear', n_jobs = n_jobs)
        model.fit(X_train, y_train)
        y_val_predicted = model.predict_proba(X_val)[:, 1]
        auprc = average_precision_score(y_val, y_val_predicted)
        return auprc

    return objective

def svm_objective(X_train: np.ndarray,
                  y_train: np.ndarray,
                  X_val: np.ndarray,
                  y_val: np.ndarray,
                  n_jobs: int = 1):
    def objective(C, gamma):
        # SVC with RBF kernel
        model = SVC(C=float(C), gamma=float(gamma), probability=True)
        model.fit(X_train, y_train)
        y_val_predicted = model.predict_proba(X_val)[:, 1]
        auprc = average_precision_score(y_val, y_val_predicted)
        return auprc
    return objective

def random_forest_objective(X_train: np.ndarray,
                            y_train: np.ndarray,
                            X_val: np.ndarray,
                            y_val: np.ndarray,
                            n_jobs: int = 1):
    def objective(n_estimators, max_depth, min_samples_split):
        model = RandomForestClassifier(
            n_estimators=int(n_estimators),
            max_depth=int(max_depth),
            min_samples_split=int(min_samples_split),
            n_jobs=n_jobs)
        model.fit(X_train, y_train)
        y_val_predicted = model.predict_proba(X_val)[:, 1]
        auprc = average_precision_score(y_val, y_val_predicted)
        return auprc
    return objective

def run_bayesian_hyperparameter_search(model_type: str,
                                       X_train: np.ndarray,
                                       y_train: np.ndarray,
                                       X_val: np.ndarray,
                                       y_val: np.ndarray,
                                       random_state: int,
                                       init_points: int,
                                       n_iter: int,
                                       n_jobs: int = 1) -> dict:
    """
    Conduct hyperparameter optimization using Bayesian optimization for various ML architectures
    """
    if model_type == 'XGBoost':
        objective = XGBC_objective(X_train, y_train, X_val, y_val, n_jobs)

        # Define the bounds for each hyperparameter
        pbounds = {
            'learning_rate': (0.1, 0.5),
            'max_leaves': (20, 300),
            'max_depth': (1, 15),
            'reg_alpha': (0, 1.0),
            'reg_lambda': (0, 1.0),
            'n_estimators': (20, 300),
            'min_child_weight': (2, 10),
            'colsample_bytree': (0.5, 1.0),
            'colsample_bylevel': (0.5, 1.0),
            'colsample_bynode': (0.5, 1.0),
            'subsample': (0.4, 1.0),
            'scale_pos_weight': (1, 5)
        }

        optimizer = BayesianOptimization(f = objective,
                                         pbounds = pbounds,
                                         random_state = random_state)

        optimizer.maximize(
            init_points = init_points,  # number of randomly chosen points to sample the target function before fitting the GP
            n_iter = n_iter)  # total number of times the process is to be repeated

        best_params = optimizer.max['params']
        best_score = optimizer.max['target']

        print(f"Best AUPRC: {best_score:.4f} achieved with {best_params}")

        return best_params

    if model_type == 'LGBM':
        objective = LGBC_objective(X_train, y_train, X_val, y_val, n_jobs)

        pbounds = { 'learning_rate': (0.01, 1.0),#(0.1, 0.5),
                    'num_leaves': (20, 300),
                    'max_depth': (1, 15),
                    'reg_alpha': (0, 1.0),
                    'reg_lambda': (0, 1.0),
                    'n_estimators': (20, 300),
                    'min_child_weight': (1,10),#(2, 10),
                    'colsample_bytree': (0, 1.0),#(0.5, 1.0),
                    'subsample': (0, 1.0),#(0.4, 1.0),
                    'scale_pos_weight': (1, 5)}

        optimizer = BayesianOptimization(f = objective, pbounds = pbounds, random_state = random_state)
        optimizer.maximize(init_points = init_points, n_iter = n_iter)
        best_params = optimizer.max['params']
        best_score = optimizer.max['target']
        print(f"Best AUPRC: {best_score:.4f} achieved with {best_params}")
        return best_params

    if model_type == 'Logistic':
        objective = logistic_regression_objective(X_train, y_train, X_val, y_val, n_jobs)
        pbounds = {'C': (0.01, 10)} # define bounds for regularization strength
        optimizer = BayesianOptimization(f=objective, pbounds=pbounds, random_state=random_state)
        optimizer.maximize(init_points=init_points, n_iter=n_iter)
        best_params = optimizer.max['params']
        best_score = optimizer.max['target']
        print(f"Best AUPRC: {best_score:.4f} achieved with {best_params}")
        return best_params

    if model_type == 'SVM':
        objective = svm_objective(X_train, y_train, X_val, y_val, n_jobs)
        pbounds = {'C': (0.1, 100),'gamma': (0.01, 1)}
        optimizer = BayesianOptimization(f=objective, pbounds=pbounds, random_state=random_state)
        optimizer.maximize(init_points=init_points, n_iter=n_iter)
        best_params = optimizer.max['params']
        best_score = optimizer.max['target']
        print(f"Best AUPRC: {best_score:.4f} achieved with {best_params}")
        return best_params

    if model_type == 'Random_forest':
        objective = random_forest_objective(X_train, y_train, X_val, y_val, n_jobs)
        pbounds = {'n_estimators': (100, 1000), 'max_depth': (5, 50), 'min_samples_split': (2, 20)}
        optimizer = BayesianOptimization(f=objective, pbounds=pbounds, random_state=random_state)
        optimizer.maximize(init_points=init_points, n_iter=n_iter)
        best_params = optimizer.max['params']
        best_score = optimizer.max['target']
        print(f"Best AUPRC: {best_score:.4f} achieved with {best_params}")
        return best_params

def train_XGBoost_model(random_state: int,
                        n_jobs: int,
                        opt_hyerparams: dict,
                        model_used_for_tuning: str,
                        X_train: np.ndarray,
                        y_train: np.ndarray) -> any:

    if model_used_for_tuning == 'XGBoost':

        feasibility_xgboost = XGBClassifier(objective = 'binary:logistic',
                                            random_state = random_state,
                                            n_jobs = int(n_jobs),
                                            max_leaves = int(opt_hyerparams['max_leaves']),
                                            learning_rate = opt_hyerparams['learning_rate'],
                                            max_depth = int(opt_hyerparams['max_depth']),
                                            reg_alpha = opt_hyerparams['reg_alpha'],
                                            reg_lambda = opt_hyerparams['reg_lambda'],
                                            n_estimators = int(opt_hyerparams['n_estimators']),
                                            min_child_weight = opt_hyerparams['min_child_weight'],
                                            colsample_bytree = opt_hyerparams['colsample_bytree'],
                                            colsample_bylevel = opt_hyerparams['colsample_bylevel'],
                                            colsample_bynode = opt_hyerparams['colsample_bynode'],
                                            subsample = opt_hyerparams['subsample'],
                                            scale_pos_weight = opt_hyerparams['scale_pos_weight'])

        feasibility_xgboost.fit(X_train, y_train)

        return feasibility_xgboost

    if model_used_for_tuning == 'lightGBM':
        feasibility_xgboost = XGBClassifier(objective = 'binary:logistic',
                                            random_state = random_state,
                                            n_jobs = int(n_jobs),
                                            learning_rate = opt_hyerparams['learning_rate'],
                                            max_leaves = int(opt_hyerparams['num_leaves']),
                                            max_depth = int(opt_hyerparams['max_depth']),
                                            reg_alpha = opt_hyerparams['reg_alpha'],
                                            reg_lambda = opt_hyerparams['reg_lambda'],
                                            n_estimators = int(opt_hyerparams['n_estimators']),
                                            min_child_weight = opt_hyerparams['min_child_weight'],
                                            colsample_bytree = opt_hyerparams['colsample_bytree'],
                                            colsample_bylevel = 1.0,
                                            colsample_bynode = 1.0,
                                            scale_pos_weight = opt_hyerparams['scale_pos_weight'],
                                            subsample = opt_hyerparams['subsample'])

        feasibility_xgboost.fit(X_train, y_train)
        return feasibility_xgboost

def train_logistic_model(random_state: int,
                         n_jobs: int,
                         opt_hyperparams: dict,
                         X_train: np.ndarray,
                         y_train: np.ndarray) -> any:
    """
    Trains a Logistic Regression model using optimal hyperparameters obtained from Bayesian optimization.

    Parameters:
        random_state (int): Seed used by the random number generator.
        n_jobs (int): Number of CPU cores used when parallelizing over classes.
        opt_hyperparams (dict): Dictionary containing the optimal hyperparameters.
        X_train (np.ndarray): Training data features.
        y_train (np.ndarray): Training data labels.

    Returns:
        Trained Logistic Regression model.
    """

    # Instantiate the Logistic Regression model with optimal hyperparameters
    logistic_model = LogisticRegression(
        C = opt_hyperparams['C'],  # Regularization strength
        penalty = 'l2',  # or opt_hyperparams['penalty'] if penalty type was also optimized
        solver = 'liblinear',  # Suitable solver for small datasets and for 'l1' and 'l2'
        random_state = random_state,
        n_jobs = n_jobs)

    # Fit the model on the training data
    logistic_model.fit(X_train, y_train)

    return logistic_model

def train_SVM_model(random_state: int,
                    opt_hyperparams: dict,
                    X_train: np.ndarray,
                    y_train: np.ndarray) -> any:
    """
    Trains a Support Vector Machine (SVM) model using optimal hyperparameters obtained from Bayesian optimization.

    Parameters:
        random_state (int): Seed used by the random number generator.
        n_jobs (int): Number of CPU cores used for the training process.
        opt_hyperparams (dict): Dictionary containing the optimal hyperparameters.
        model_used_for_tuning (str): Specifies the model type to ensure correct hyperparameter assignment.
        X_train (np.ndarray): Training data features.
        y_train (np.ndarray): Training data labels.

    Returns:
        Trained SVM model.
    """
    # Instantiate the SVM model with the RBF kernel and optimal hyperparameters
    svm_model = SVC(
        C = opt_hyperparams['C'],  # Penalty parameter
        gamma = opt_hyperparams['gamma'],  # Kernel coefficient for RBF
        kernel = 'rbf',  # You can make this dynamic if needed
        probability = True,  # To enable probability estimates
        random_state = random_state)

    # Fit the model on the training data
    svm_model.fit(X_train, y_train)

    return svm_model

def train_Random_forest_model(random_state: int,
                              n_jobs: int,
                              opt_hyperparams: dict,
                              X_train: np.ndarray,
                              y_train: np.ndarray) -> any:
    """
    Trains a Random Forest model using optimal hyperparameters obtained from Bayesian optimization.

    Parameters:
        random_state (int): Seed used by the random number generator.
        n_jobs (int): Number of CPU cores used for the training process.
        opt_hyperparams (dict): Dictionary containing the optimal hyperparameters.
        model_used_for_tuning (str): Specifies the model type to ensure correct hyperparameter assignment.
        X_train (np.ndarray): Training data features.
        y_train (np.ndarray): Training data labels.

    Returns:
        Trained Random Forest model.
     """
    # Instantiate the Random Forest classifier with optimal hyperparameters
    random_forest_model = RandomForestClassifier(
                                    n_estimators=int(opt_hyperparams['n_estimators']),
                                    max_depth=int(opt_hyperparams['max_depth']),
                                    min_samples_split=int(opt_hyperparams['min_samples_split']),
                                    random_state=random_state,
                                    n_jobs=n_jobs)

    # Fit the model on the training data
    random_forest_model.fit(X_train, y_train)

    return random_forest_model

def bootstrap_auprc(y_true, y_pred_proba, n_iterations=1000):
    """
    Calculate the mean Area Under the Precision-Recall Curve (AUPRC) and its 95% confidence intervals
    based on bootstrap resampling.

    Parameters:
    - y_true (array-like): True binary labels.
    - y_pred_proba (array-like): Predicted probabilities for the positive class.
    - n_iterations (int, optional): Number of bootstrap iterations. Default is 1000.

    Returns:
    - tuple: mean AUPRC and its 95% confidence interval (lower, upper).
    """
    scores = []
    n_samples = len(y_true)

    for _ in range(n_iterations):
        # Randomly sample indices with replacement
        indices = np.random.choice(n_samples, n_samples, replace=True)
        score = average_precision_score(y_true[indices], y_pred_proba[indices])
        scores.append(score)

    # Calculate mean and confidence intervals
    mean_score = np.mean(scores)
    lower = np.percentile(scores, 2.5)  # 2.5 percentile
    upper = np.percentile(scores, 97.5)  # 97.5 percentile

    print(f"\nMean AUPRC: {mean_score:.4f}, 95% CI: ({lower:.4f}, {upper:.4f})")

    return mean_score, (lower, upper)

def bootstrap_precision(y_true, y_pred_binary, n_iterations=1000):
    """
    Calculate the mean precision score and its 95% confidence intervals based on
    bootstrap resampling. If there are no positive samples, returns 0.

    Parameters:
    - y_true (array-like): True binary labels.
    - y_pred_binary (array-like): Predicted binary labels.
    - n_iterations (int, optional): Number of bootstrap iterations. Default is 1000.

    Returns:
    - tuple: mean precision and its 95% confidence interval (lower, upper).
    """
    scores = []
    n_samples = len(y_true)

    for _ in range(n_iterations):
        # Randomly sample indices with replacement
        indices = np.random.choice(n_samples, n_samples, replace=True)
        try:
            score = precision_score(y_true[indices], y_pred_binary[indices])
            scores.append(score)
        except ValueError:
            scores.append(0.0)  # or continue to skip

    # Calculate mean and confidence intervals
    mean_score = np.mean(scores)
    lower = np.percentile(scores, 2.5)  # 2.5 percentile
    upper = np.percentile(scores, 97.5)  # 97.5 percentile

    print(f"\nMean precision: {mean_score:.4f}, 95% CI: ({lower:.4f}, {upper:.4f})")

    return mean_score, (lower, upper)

def bootstrap_accuracy(y_true, y_pred_binary, n_iterations=1000):
    """
    Calculate the mean accuracy score and its 95% confidence intervals based on
    bootstrap resampling.

    Parameters:
    - y_true (array-like): True binary labels.
    - y_pred_binary (array-like): Predicted binary labels.
    - n_iterations (int, optional): Number of bootstrap iterations. Default is 1000.

    Returns:
    - tuple: mean accuracy and its 95% confidence interval (lower, upper).
    """
    scores = []
    n_samples = len(y_true)

    for _ in range(n_iterations):
        # Randomly sample indices with replacement
        indices = np.random.choice(n_samples, n_samples, replace=True)
        score = accuracy_score(y_true[indices], y_pred_binary[indices])
        scores.append(score)
        scores.append(0.0)  # or continue to skip

    # Calculate mean and confidence intervals
    mean_score = np.mean(scores)
    lower = np.percentile(scores, 2.5)  # 2.5 percentile
    upper = np.percentile(scores, 97.5)  # 97.5 percentile

    print(f"\nMean accuracy: {mean_score:.4f}, 95% CI: ({lower:.4f}, {upper:.4f})")

    return mean_score, (lower, upper)

def bootstrap_f1(y_true, y_pred_binary, n_iterations=1000):
    """
    Calculate the mean F1 score and its 95% confidence intervals based on bootstrap resampling.

    Parameters:
    - y_true (array-like): True binary labels.
    - y_pred_binary (array-like): Predicted binary labels.
    - n_iterations (int, optional): Number of bootstrap iterations. Default is 1000.

    Returns:
    - tuple: mean F1 score and its 95% confidence interval (lower, upper).
    """
    scores = []
    n_samples = len(y_true)

    for _ in range(n_iterations):
        # Randomly sample indices with replacement
        indices = np.random.choice(n_samples, n_samples, replace=True)
        score = f1_score(y_true[indices], y_pred_binary[indices], zero_division=0)
        scores.append(score)

    # Calculate mean and confidence intervals
    mean_score = np.mean(scores)
    lower = np.percentile(scores, 2.5)  # 2.5 percentile
    upper = np.percentile(scores, 97.5)  # 97.5 percentile

    print(f"\nMean F1: {mean_score:.4f}, 95% CI: ({lower:.4f}, {upper:.4f})")

    return mean_score, (lower, upper)

def bootstrap_recall(y_true, y_pred_binary, n_iterations=1000):
    """
    Calculate the mean recall score and its 95% confidence intervals based on bootstrap resampling.

    Parameters:
    - y_true (array-like): True binary labels.
    - y_pred_binary (array-like): Predicted binary labels.
    - n_iterations (int, optional): Number of bootstrap iterations. Default is 1000.

    Returns:
    - tuple: mean recall and its 95% confidence interval (lower, upper).
    """
    scores = []
    n_samples = len(y_true)

    for _ in range(n_iterations):
        # Randomly sample indices with replacement
        indices = np.random.choice(n_samples, n_samples, replace=True)
        score = recall_score(y_true[indices], y_pred_binary[indices], zero_division=0)
        scores.append(score)

    # Calculate mean and confidence intervals
    mean_score = np.mean(scores)
    lower = np.percentile(scores, 2.5)  # 2.5 percentile
    upper = np.percentile(scores, 97.5)  # 97.5 percentile

    print(f"\nMean Recall: {mean_score:.4f}, 95% CI: ({lower:.4f}, {upper:.4f})")

    return mean_score, (lower, upper)