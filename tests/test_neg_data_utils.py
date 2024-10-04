"""
All tests here are for functions in our custom neg_data_utils package
For tests that actually check the synthetically generated negative reactions, refer to:
--> 04/data/all_non_fingerprinted_data_tests.py
"""

from minedatabase.pickaxe import Pickaxe
from neg_data_utils import neg_data_utils
import pandas as pd
import pickle
import os

# Change directory and run this from in the scripts' folder - since we will generate negative data from this folder also
os.chdir('../../scripts')

# -------- Test to select specific reaction rules for a Pickaxe expansion and create their .tsv files  --------
def test_select_rule_1_for_expansion():

    mapped_rule = 'rule0001'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0001.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_2_for_expansion():

    mapped_rule = 'rule0002'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0002.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_3_for_expansion():

    mapped_rule = 'rule0003'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0003.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_4_for_expansion():

    mapped_rule = 'rule0004'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0004.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_5_for_expansion():

    mapped_rule = 'rule0005'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0005.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_83_for_expansion():

    mapped_rule = 'rule0083'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0083.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_798_for_expansion():

    mapped_rule = 'rule0798'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0798.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_999_for_expansion():

    mapped_rule = 'rule0999'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule0999.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_1024_for_expansion():

    mapped_rule = 'rule1024'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule1024.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

def test_select_rule_1224_for_expansion():

    mapped_rule = 'rule1224'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule = mapped_rule,
                                             gen_rules_filepath = gen_rules_filepath)

    assert rule_filepath == "../data/interim/rules_for_pickaxe/rule1224.tsv"

    df = pd.read_csv(rule_filepath, delimiter='\t').drop(labels=['Unnamed: 0'], axis=1)
    nrows, ncols = df.shape

    # assert the number of rows is 1 because this functions creates a dataframe for a single rule only
    assert nrows == 1

    # number of columns for this single rule dataframe should still be the same as the number of columns in JN1224MIN
    assert ncols == pd.read_csv(gen_rules_filepath, delimiter='\t').shape[1]

# -------- Test to create tsv files of compounds that will undergo a single-step Pickaxe expansion  --------
def test_create_tsv_for_input_cpd_1():
    input_cpd_filepath = neg_data_utils.create_tsv_for_input_cpd(input_cpd_SMILES='CC')

    assert input_cpd_filepath == '../data/interim/cpds_for_pickaxe/CC.tsv'

    df = pd.read_csv(input_cpd_filepath,delimiter='\t')
    nrows, ncols = df.shape
    assert nrows == 1
    assert ncols == 2

def test_create_tsv_for_input_cpd_2():
    input_cpd_filepath = neg_data_utils.create_tsv_for_input_cpd(input_cpd_SMILES='CCC')

    assert input_cpd_filepath == '../data/interim/cpds_for_pickaxe/CCC.tsv'

    df = pd.read_csv(input_cpd_filepath,delimiter='\t')
    nrows, ncols = df.shape
    assert nrows == 1
    assert ncols == 2

def test_create_tsv_for_input_cpd_3():
    input_cpd_filepath = neg_data_utils.create_tsv_for_input_cpd(input_cpd_SMILES='C1=CC=CC=C1')

    assert input_cpd_filepath == '../data/interim/cpds_for_pickaxe/C1=CC=CC=C1.tsv'

    df = pd.read_csv(input_cpd_filepath,delimiter='\t')
    nrows, ncols = df.shape
    assert nrows == 1
    assert ncols == 2

def test_create_tsv_for_input_cpd_4():
    input_cpd_filepath = neg_data_utils.create_tsv_for_input_cpd(input_cpd_SMILES='C12=CC=CC=C1CCCC2')

    assert input_cpd_filepath == '../data/interim/cpds_for_pickaxe/C12=CC=CC=C1CCCC2.tsv'

    df = pd.read_csv(input_cpd_filepath,delimiter='\t')
    nrows, ncols = df.shape
    assert nrows == 1
    assert ncols == 2

def test_create_tsv_for_input_cpd_5():
    input_cpd_filepath = neg_data_utils.create_tsv_for_input_cpd(input_cpd_SMILES='Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1[nH]c(=O)[nH]c(=O)c1N2')

    assert input_cpd_filepath == '../data/interim/cpds_for_pickaxe/Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1[nH]c(=O)[nH]c(=O)c1N2.tsv'

    df = pd.read_csv(input_cpd_filepath, delimiter='\t')
    nrows, ncols = df.shape
    assert nrows == 1
    assert ncols == 2

# -------- Test to create a pandas dataframe of compounds generated after a Pickaxe expansion  --------
def test_create_cpds_df_aft_expansion_rule0002_1():

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'
    mapped_rule = 'rule0002'
    reported_substrate = 'O=C(CO)C(O)C(O)C(O)CO' # contains 5 hydroxyl groups so expect 5 products with rule0002

    # create .tsv file for the rule that this reported reaction has been mapped onto (should be exactly one rule)
    input_rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule=mapped_rule,
                                                                   gen_rules_filepath=gen_rules_filepath)

    # create .tsv file for the substrate that will be used in Pickaxe expansions to enumerate negative products
    input_cpd_filepath = neg_data_utils.create_tsv_for_input_cpd(input_cpd_SMILES=reported_substrate)

    # initialize Pickaxe object then run single-step metabolic expansion
    pk = Pickaxe(coreactant_list=coreactants_filepath, rule_list=input_rule_filepath)
    pk.load_compound_set(compound_file=input_cpd_filepath)  # load starting compound
    pk.transform_all(processes=1, generations=1)  # expand on this substrate for one step only

    compounds_df = neg_data_utils.create_cpds_df_aft_expansion(pk)

    # first, a simple check to ensure that there is no compound beyond a single generation
    assert all(compounds_df["Generation"] <= 1)

    # then, I manually checked that only 5 new products should be generated for this reaction
    products_df = compounds_df[(compounds_df['Type'] == 'Predicted') & (compounds_df['Generation'] == 1)]

    # as such, check that exactly 5 products are generated
    assert products_df.shape[0] == 5

def test_create_cpds_df_aft_expansion_rule0002_2():

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'
    mapped_rule = 'rule0002'
    reported_substrate = 'OCC(O)C(O)C(O)C(O)CO'

    # create .tsv file for the rule that this reported reaction has been mapped onto (should be exactly one rule)
    input_rule_filepath = neg_data_utils.select_rule_for_expansion(mapped_rule=mapped_rule,
                                                                   gen_rules_filepath=gen_rules_filepath)

    # create .tsv file for the substrate that will be used in Pickaxe expansions to enumerate negative products
    input_cpd_filepath = neg_data_utils.create_tsv_for_input_cpd(input_cpd_SMILES=reported_substrate)

    # initialize Pickaxe object then run single-step metabolic expansion
    pk = Pickaxe(coreactant_list=coreactants_filepath, rule_list=input_rule_filepath)
    pk.load_compound_set(compound_file=input_cpd_filepath)  # load starting compound
    pk.transform_all(processes=1, generations=1)  # expand on this substrate for one step only

    compounds_df = neg_data_utils.create_cpds_df_aft_expansion(pk)

    # first, a simple check to ensure that there is no compound beyond a single generation
    assert all(compounds_df["Generation"] <= 1)

    # then, I manually checked that only 6 new products should be generated for this reaction
    products_df = compounds_df[(compounds_df['Type'] == 'Predicted') & (compounds_df['Generation'] == 1)]

    # as such, check that exactly 3 products are generated
    assert products_df.shape[0] == 3

# -------- Test to check 'run_pickaxe_if_single_reported_product' function --------

### Test synthetically generating infeasible reactions from feasible alcohol dehydrogenase reactions

def test_gen_neg_products_if_single_reported_product_1_AdH_rxns():

    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[0]

    neg_rxns_list = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                          coreactants_filepath = coreactants_filepath,
                                                          reported_substrate = pos_rxn_entry['Substrates'][0],
                                                          reported_product = pos_rxn_entry['Products'][0],
                                                          gen_rules_filepath = gen_rules_filepath)

    # this reaction should not produce any alternate products because it involves reduction of 'O=C(CO)C(O)C(O)C(O)CO'
    # which has only one reactive site anyway
    assert neg_rxns_list == []

def test_gen_neg_products_if_single_reported_product_2_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[2]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath = coreactants_filepath,
                                                                             reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                             reported_product = pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate just one new negative reaction
    assert len(syn_gen_neg_data) == 1

    # check that original reaction metadata has not been changed
    for neg_data_entry in syn_gen_neg_data:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert syn_gen_neg_data == [{'Original enzymatic database': 'metacyc',
                                  'Original reaction ID': '1.1.1.51-RXN',
                                  'Original substrates': ['CC12CCC3C(CCC4=CC(=O)CCC43C)C1CCC2=O'],
                                  'Original products': ['CC12CCC(=O)C=C1CCC1C2CCC2(C)C(O)CCC12'],
                                  'Alternate products': ['CC12CCC3C(CCC4=CC(O)CCC43C)C1CCC2=O'],
                                  'LHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
                                  'RHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
                                  'Reaction eq': 'CC12CCC3C(CCC4=CC(=O)CCC43C)C1CCC2=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = CC12CCC3C(CCC4=CC(O)CCC43C)C1CCC2=O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
                                  'Original rule': 'rule0003',
                                  'feasibility_label': 0,
                                  'remark': 'alternate product from pos rxn'}]

def test_gen_neg_products_if_single_reported_product_3_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[3]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate no new reactions
    assert len(syn_gen_neg_data) == 0

    # check that original reaction metadata has not been changed
    assert syn_gen_neg_data == []

def test_gen_neg_products_if_single_reported_product_4_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[4]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate no new reactions
    assert len(syn_gen_neg_data) == 0

    # check that original reaction metadata has not been changed
    assert syn_gen_neg_data == []

def test_gen_neg_products_if_single_reported_product_5_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[15]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate one negative reactions
    assert len(syn_gen_neg_data) == 1

    # check that original reaction metadata has not been changed
    for neg_data_entry in syn_gen_neg_data:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert syn_gen_neg_data == [{'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-7784_reverse',
  'Original substrates': ['Oc1cc(O)c2c(c1)OC(c1cc(O)c(O)c(O)c1)C(O)C2O'],
  'Original products': ['O=C1c2c(O)cc(O)cc2OC(c2cc(O)c(O)c(O)c2)C1O'],
  'Alternate products': ['O=C1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1'],
  'LHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'RHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'Reaction eq': 'Oc1cc(O)c2c(c1)OC(c1cc(O)c(O)c(O)c1)C(O)C2O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1 = O=C1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1 + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1',
  'Original rule': 'rule0002',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'}]

def test_gen_neg_products_if_single_reported_product_6_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[28]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate three negative reactions
    assert len(syn_gen_neg_data) == 3

    # check that original reaction metadata has not been changed
    for neg_data_entry in syn_gen_neg_data:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
              'Original reaction ID': 'RXN-16095',
              'Original substrates': ['CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
              'Original products': ['CCCCCC=CCC=CCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
              'Alternate products': ['CCCCCC=CCC=CCCCCCCCC(=O)CC(O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
              'LHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
              'RHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
              'Reaction eq': 'CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CCCCCC=CCC=CCCCCCCCC(=O)CC(O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
              'Original rule': 'rule0003',
              'feasibility_label': 0,
              'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

    assert {'Original enzymatic database': 'metacyc',
          'Original reaction ID': 'RXN-16095',
          'Original substrates': ['CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
          'Original products': ['CCCCCC=CCC=CCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
          'Alternate products': ['CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
          'LHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
          'RHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
          'Reaction eq': 'CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
          'Original rule': 'rule0003',
          'feasibility_label': 0,
          'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-16095',
  'Original substrates': ['CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
  'Original products': ['CCCCCC=CCC=CCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
  'Alternate products': ['CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
  'LHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CCCCCC=CCC=CCCCCCCCC(=O)CC(=O)SCCNC(O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0003',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

def test_gen_neg_products_if_single_reported_product_7_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[47]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate three negative reactions
    assert len(syn_gen_neg_data) == 2

    # check that original reaction metadata has not been changed
    for neg_data_entry in syn_gen_neg_data:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
              'Original reaction ID': 'RXN-9796_reverse',
              'Original substrates': ['CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(O)CCC4(C)C3CC(O)C12C'],
              'Original products': ['CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(=O)CCC4(C)C3CC(O)C12C'],
              'Alternate products': ['CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(O)CCC4(C)C3CC(=O)C12C'],
              'LHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
              'RHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
              'Reaction eq': 'CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(O)CCC4(C)C3CC(O)C12C + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(O)CCC4(C)C3CC(=O)C12C + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1',
              'Original rule': 'rule0002',
              'feasibility_label': 0,
              'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

    assert {'Original enzymatic database': 'metacyc',
              'Original reaction ID': 'RXN-9796_reverse',
              'Original substrates': ['CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(O)CCC4(C)C3CC(O)C12C'],
              'Original products': ['CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(=O)CCC4(C)C3CC(O)C12C'],
              'Alternate products': ['CC(C)CCCC(C)C1CCC2C3C(=O)CC4CC(O)CCC4(C)C3CC(O)C12C'],
              'LHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
              'RHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
              'Reaction eq': 'CC(C)CCCC(C)C1CCC2C3C(O)CC4CC(O)CCC4(C)C3CC(O)C12C + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = CC(C)CCCC(C)C1CCC2C3C(=O)CC4CC(O)CCC4(C)C3CC(O)C12C + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1',
              'Original rule': 'rule0002',
              'feasibility_label': 0,
              'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

def test_gen_neg_products_if_single_reported_product_8_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[68]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate no new reactions
    assert len(syn_gen_neg_data) == 0

    # check that original reaction metadata has not been changed
    assert syn_gen_neg_data == []

def test_gen_neg_products_if_single_reported_product_9_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[350]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate one new negative reaction
    assert len(syn_gen_neg_data) == 1

    # check that original reaction metadata has not been changed
    for neg_data_entry in syn_gen_neg_data:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
          'Original reaction ID': 'RXN-17409',
          'Original substrates': ['CC(CCC(=O)O)C1CCC2C3CCC4CC(=O)CCC4(C)C3CCC12C'],
          'Original products': ['CC(CCC(=O)O)C1CCC2C3CCC4CC(O)CCC4(C)C3CCC12C'],
          'Alternate products': ['CC(CCC(O)O)C1CCC2C3CCC4CC(=O)CCC4(C)C3CCC12C'],
          'LHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
          'RHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
          'Reaction eq': 'CC(CCC(=O)O)C1CCC2C3CCC4CC(=O)CCC4(C)C3CCC12C + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = CC(CCC(O)O)C1CCC2C3CCC4CC(=O)CCC4(C)C3CCC12C + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
          'Original rule': 'rule0003',
          'feasibility_label': 0,
          'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

def test_gen_neg_products_if_single_reported_product_10_AdH_rxns():
    with open('../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[500]

    syn_gen_neg_data = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                                             coreactants_filepath=coreactants_filepath,
                                                                             reported_substrate=pos_rxn_entry['Substrates'][0],
                                                                             reported_product=pos_rxn_entry['Products'][0],
                                                                             gen_rules_filepath=gen_rules_filepath)

    # this reaction should generate four new negative reactions
    assert len(syn_gen_neg_data) == 4

    # check that original reaction metadata has not been changed
    for neg_data_entry in syn_gen_neg_data:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
              'Original reaction ID': 'RXN-15924_reverse',
              'Original substrates': ['OCC1OC(O)C(O)C(O)C1O'],
              'Original products': ['O=C1OC(CO)C(O)C(O)C1O'],
              'Alternate products': ['O=C1C(CO)OC(O)C(O)C1O'],
              'LHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
              'RHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
              'Reaction eq': 'OCC1OC(O)C(O)C(O)C1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1 = O=C1C(CO)OC(O)C(O)C1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1',
              'Original rule': 'rule0002',
              'feasibility_label': 0,
              'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

    assert {'Original enzymatic database': 'metacyc',
          'Original reaction ID': 'RXN-15924_reverse',
          'Original substrates': ['OCC1OC(O)C(O)C(O)C1O'],
          'Original products': ['O=C1OC(CO)C(O)C(O)C1O'],
          'Alternate products': ['O=C1C(O)OC(CO)C(O)C1O'],
          'LHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
          'RHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
          'Reaction eq': 'OCC1OC(O)C(O)C(O)C1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1 = O=C1C(O)OC(CO)C(O)C1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1',
          'Original rule': 'rule0002',
          'feasibility_label': 0,
          'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

    assert {'Original enzymatic database': 'metacyc',
              'Original reaction ID': 'RXN-15924_reverse',
              'Original substrates': ['OCC1OC(O)C(O)C(O)C1O'],
              'Original products': ['O=C1OC(CO)C(O)C(O)C1O'],
              'Alternate products': ['O=CC1OC(O)C(O)C(O)C1O'],
              'LHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
              'RHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
              'Reaction eq': 'OCC1OC(O)C(O)C(O)C1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1 = O=CC1OC(O)C(O)C(O)C1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1',
              'Original rule': 'rule0002',
              'feasibility_label': 0,
              'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

    assert {'Original enzymatic database': 'metacyc',
          'Original reaction ID': 'RXN-15924_reverse',
          'Original substrates': ['OCC1OC(O)C(O)C(O)C1O'],
          'Original products': ['O=C1OC(CO)C(O)C(O)C1O'],
          'Alternate products': ['O=C1C(O)C(O)OC(CO)C1O'],
          'LHS_cofactors': ['NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
          'RHS_cofactors': ['NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
          'Reaction eq': 'OCC1OC(O)C(O)C(O)C1O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1 = O=C1C(O)C(O)OC(CO)C1O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1',
          'Original rule': 'rule0002',
          'feasibility_label': 0,
          'remark': 'alternate product from pos rxn'} in syn_gen_neg_data

def test_gen_neg_products_if_single_reported_product_1_Monox_rxns():

    with open('../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_10.pkl','rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[3]

    neg_rxns_list = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                          coreactants_filepath = coreactants_filepath,
                                                          reported_substrate = pos_rxn_entry['Substrates'][0],
                                                          reported_product = pos_rxn_entry['Products'][0],
                                                          gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 15 new negative reactions
    assert len(neg_rxns_list) == 15

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(CO)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(CO)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3C(CCC3(C)O)C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3C(CCC3(C)O)C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert  {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(CO)CC3C(C)CCC3C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(CO)CC3C(C)CCC3C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1C2=CCC(C)(O)C3CCC(C)C3CC2(C)CC1O'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1C2=CCC(C)(O)C3CCC(C)C3CC2(C)CC1O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3(O)C(C)CCC3C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3(O)C(C)CCC3C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC1CCC2C1CC1(C)CCC(O)(C(C)C)C1=CCC2(C)O'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC1CCC2C1CC1(C)CCC(O)(C(C)C)C1=CCC2(C)O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC(O)=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC(O)=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3C(C)CCC3(O)C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3C(C)CCC3(O)C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert  {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3C(C)C(O)CC3C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3C(C)C(O)CC3C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)C1=CCC(C)(O)C1CCC(C)C1C2O'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)C1=CCC(C)(O)C1CCC(C)C1C2O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)C(O)C=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)C(O)C=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(O)(CO)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3C(C)CCC3C(O)(CO)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CC(O)C2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CC(O)C2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC(C)C1CCC2(C)CC3C(CO)CCC3C(C)(O)CC=C12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC(C)C1CCC2(C)CC3C(CO)CCC3C(C)(O)CC=C12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-15431_reverse',
  'Original substrates': ['CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12'],
  'Original products': ['CC(C)C1CCC2(C)CC3C(C)CC(O)C3C(C)(O)CC=C12'],
  'Alternate products': ['CC1CCC2C1CC1(C)CCC(C(C)(C)O)C1=CCC2(C)O'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'CC(C)C1CCC2(C)CC3C(C)CCC3C(C)(O)CC=C12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC1CCC2C1CC1(C)CCC(C(C)(C)O)C1=CCC2(C)O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_reported_product_2_Monox_rxns():

    with open('../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_10.pkl','rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[5]

    neg_rxns_list = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                          coreactants_filepath = coreactants_filepath,
                                                          reported_substrate = pos_rxn_entry['Substrates'][0],
                                                          reported_product = pos_rxn_entry['Products'][0],
                                                          gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 15 new negative reactions
    assert len(neg_rxns_list) == 9

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
          'Original reaction ID': 'RXN-8707_reverse',
          'Original substrates': ['NCCc1c[nH]c2ccccc12'],
          'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
          'Alternate products': ['NCCc1c(O)[nH]c2ccccc12'],
          'LHS_cofactors': ['O=O',
           'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
          'RHS_cofactors': ['O',
           'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
          'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCCc1c(O)[nH]c2ccccc12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
          'Original rule': 'rule0004',
          'feasibility_label': 0,
          'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NCC(O)c1c[nH]c2ccccc12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCC(O)c1c[nH]c2ccccc12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NCCc1c[nH]c2cc(O)ccc12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCCc1c[nH]c2cc(O)ccc12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NCCC1(O)C=Nc2ccccc21'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCCC1(O)C=Nc2ccccc21 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NCCC1=CN=C2C=CC=CC12O'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCCC1=CN=C2C=CC=CC12O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NC(O)Cc1c[nH]c2ccccc12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NC(O)Cc1c[nH]c2ccccc12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NCCC1=C2C=CC=CC2(O)N=C1'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCCC1=C2C=CC=CC2(O)N=C1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NCCc1c[nH]c2c(O)cccc12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCCc1c[nH]c2c(O)cccc12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
  'Original reaction ID': 'RXN-8707_reverse',
  'Original substrates': ['NCCc1c[nH]c2ccccc12'],
  'Original products': ['NCCc1c[nH]c2ccc(O)cc12'],
  'Alternate products': ['NCCc1c[nH]c2cccc(O)c12'],
  'LHS_cofactors': ['O=O',
   'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
  'RHS_cofactors': ['O',
   'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
  'Reaction eq': 'NCCc1c[nH]c2ccccc12 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = NCCc1c[nH]c2cccc(O)c12 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
  'Original rule': 'rule0004',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_reported_product_3_Monox_rxns():

    with open('../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_10.pkl','rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[9]

    neg_rxns_list = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                          coreactants_filepath = coreactants_filepath,
                                                          reported_substrate = pos_rxn_entry['Substrates'][0],
                                                          reported_product = pos_rxn_entry['Products'][0],
                                                          gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 3 new negative reactions
    assert len(neg_rxns_list) == 3

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-16949_reverse',
 'Original substrates': ['Cc1ccc(C(=O)O)cc1O'],
 'Original products': ['Cc1cc(O)c(C(=O)O)cc1O'],
 'Alternate products': ['Cc1ccc(C(=O)O)c(O)c1O'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'Cc1ccc(C(=O)O)cc1O + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = Cc1ccc(C(=O)O)c(O)c1O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-16949_reverse',
 'Original substrates': ['Cc1ccc(C(=O)O)cc1O'],
 'Original products': ['Cc1cc(O)c(C(=O)O)cc1O'],
 'Alternate products': ['Cc1c(O)cc(C(=O)O)cc1O'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'Cc1ccc(C(=O)O)cc1O + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = Cc1c(O)cc(C(=O)O)cc1O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-16949_reverse',
 'Original substrates': ['Cc1ccc(C(=O)O)cc1O'],
 'Original products': ['Cc1cc(O)c(C(=O)O)cc1O'],
 'Alternate products': ['O=C(O)c1ccc(CO)c(O)c1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'Cc1ccc(C(=O)O)cc1O + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = O=C(O)c1ccc(CO)c(O)c1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_reported_product_4_Monox_rxns():

    with open('../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_10.pkl','rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[101]

    neg_rxns_list = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                          coreactants_filepath = coreactants_filepath,
                                                          reported_substrate = pos_rxn_entry['Substrates'][0],
                                                          reported_product = pos_rxn_entry['Products'][0],
                                                          gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 3 new negative reactions
    assert len(neg_rxns_list) == 4

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C=C(O)c1ccccc1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C=C(O)c1ccccc1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C(O)=Cc1ccccc1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C(O)=Cc1ccccc1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C=Cc1cccc(O)c1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C=Cc1cccc(O)c1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C=Cc1ccc(O)cc1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C=Cc1ccc(O)cc1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_reported_product_5_Monox_rxns():

    with open('../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_10.pkl','rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[101]

    neg_rxns_list = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                          coreactants_filepath = coreactants_filepath,
                                                          reported_substrate = pos_rxn_entry['Substrates'][0],
                                                          reported_product = pos_rxn_entry['Products'][0],
                                                          gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 3 new negative reactions
    assert len(neg_rxns_list) == 4

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C=C(O)c1ccccc1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C=C(O)c1ccccc1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C(O)=Cc1ccccc1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C(O)=Cc1ccccc1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C=Cc1cccc(O)c1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C=Cc1cccc(O)c1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'TRANS-CINNAMATE-2-MONOOXYGENASE-RXN_reverse',
 'Original substrates': ['O=C(O)C=Cc1ccccc1'],
 'Original products': ['O=C(O)C=Cc1ccccc1O'],
 'Alternate products': ['O=C(O)C=Cc1ccc(O)cc1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'O=C(O)C=Cc1ccccc1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = O=C(O)C=Cc1ccc(O)cc1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_reported_product_6_Monox_rxns():
    with open('../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_10.pkl','rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    pos_rxn_entry = processed_rxns[105]

    neg_rxns_list = neg_data_utils.run_pickaxe_if_single_reported_product(pos_rxn_entry = pos_rxn_entry,
                                                          coreactants_filepath = coreactants_filepath,
                                                          reported_substrate = pos_rxn_entry['Substrates'][0],
                                                          reported_product = pos_rxn_entry['Products'][0],
                                                          gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 23 new negative reactions
    assert len(neg_rxns_list) == 21

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    # some spot checks (I didn't check for all new 21 negative reactions generated)
    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-15108_reverse',
 'Original substrates': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'Original products': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)(O)C(C(C)O)OC(=O)C2C)O1'],
 'Alternate products': ['CC1CC(N(C)CO)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC1CC(N(C)CO)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-15108_reverse',
 'Original substrates': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'Original products': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)(O)C(C(C)O)OC(=O)C2C)O1'],
 'Alternate products': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)(O)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)(O)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-15108_reverse',
 'Original substrates': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'Original products': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)(O)C(C(C)O)OC(=O)C2C)O1'],
 'Alternate products': ['CC1CC(N(C)C)C(O)C(OC2C(C)C(=O)OC(C(C)O)C(C)C=CC(=O)C(C)C(O)C2C)O1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC1CC(N(C)C)C(O)C(OC2C(C)C(=O)OC(C(C)O)C(C)C=CC(=O)C(C)C(O)C2C)O1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-15108_reverse',
 'Original substrates': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'Original products': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)(O)C(C(C)O)OC(=O)C2C)O1'],
 'Alternate products': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C(O)=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C(O)=CC(C)C(C(C)O)OC(=O)C2C)O1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-15108_reverse',
 'Original substrates': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1'],
 'Original products': ['CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)(O)C(C(C)O)OC(=O)C2C)O1'],
 'Alternate products': ['CC1CC(N(C)C)C(O)C(OC2C(C)C(=O)OC(C(C)O)C(C)C=CC(=O)C(C)CC2(C)O)O1'],
 'LHS_cofactors': ['O=O',
  'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'],
 'RHS_cofactors': ['O',
  'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'],
 'Reaction eq': 'CC1CC(N(C)C)C(O)C(OC2C(C)CC(C)C(=O)C=CC(C)C(C(C)O)OC(=O)C2C)O1 + O=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = CC1CC(N(C)C)C(O)C(OC2C(C)C(=O)OC(C(C)O)C(C)C=CC(=O)C(C)CC2(C)O)O1 + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1',
 'Original rule': 'rule0004',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_substrate_and_two_products_rule7_rxn1():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    rxns_checked = 0
    syn_gen_neg_data = None

    # this should not generate any negative reactions
    for rxn in processed_rxns:
        if rxn['Reaction ID'] == 'RXN-12253':
            syn_gen_neg_data = neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = rxn,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = rxn['Substrates'][0],
                                                                                        reported_products = rxn['Products'],
                                                                                        gen_rules_filepath=gen_rules_filepath)
            rxns_checked += 1

    assert rxns_checked == 1
    assert syn_gen_neg_data == []

def test_gen_neg_products_if_single_substrate_and_two_products_rule7_rxn2():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = None
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == 'LACTOSE6P-HYDROXY-RXN_reverse':
            pos_rxn_entry = rxn
            neg_rxns_list= neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = pos_rxn_entry,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                                        reported_products = pos_rxn_entry['Products'],
                                                                                        gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 1 new negative reactions
    assert len(neg_rxns_list) == 1

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
              'Original reaction ID': 'LACTOSE6P-HYDROXY-RXN_reverse',
              'Original substrates': ['O=P(O)(O)OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O'],
              'Original products': ['OCC1OC(O)C(O)C(O)C1O',
               'O=P(O)(O)OCC1OC(O)C(O)C(O)C1O'],
              'Alternate products': ['OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O', 'O=P(O)(O)O'],
              'LHS_cofactors': ['O'],
              'RHS_cofactors': [],
              'Reaction eq': 'O=P(O)(O)OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O + O = OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O + O=P(O)(O)O',
              'Original rule': 'rule0007',
              'feasibility_label': 0,
              'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_substrate_and_two_products_rule7_rxn3():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = None
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == 'RXN-16612_reverse':
            pos_rxn_entry = rxn
            neg_rxns_list= neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = pos_rxn_entry,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                                        reported_products = pos_rxn_entry['Products'],
                                                                                        gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 1 new negative reactions
    assert len(neg_rxns_list) == 1

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'metacyc',
          'Original reaction ID': 'RXN-16612_reverse',
          'Original substrates': ['OCC1OC(OCC2OC(OC3C(CO)OC(O)C(O)C3O)C(O)C(O)C2O)C(O)C(O)C1O'],
          'Original products': ['OCC1OC(O)C(O)C(O)C1O',
           'OCC1OC(OCC2OC(O)C(O)C(O)C2O)C(O)C(O)C1O'],
          'Alternate products': ['OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O',
           'OCC1OC(O)C(O)C(O)C1O'],
          'LHS_cofactors': ['O'],
          'RHS_cofactors': [],
          'Reaction eq': 'OCC1OC(OCC2OC(OC3C(CO)OC(O)C(O)C3O)C(O)C(O)C2O)C(O)C(O)C1O + O = OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O',
          'Original rule': 'rule0007',
          'feasibility_label': 0,
          'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_substrate_and_two_products_rule7_rxn4():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = None
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == '3.1.1.85_2':
            pos_rxn_entry = rxn
            neg_rxns_list= neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = pos_rxn_entry,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                                        reported_products = pos_rxn_entry['Products'],
                                                                                        gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 2 new negative reactions
    assert len(neg_rxns_list) == 2

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'brenda',
             'Original reaction ID': '3.1.1.85_2',
             'Original substrates': ['CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OCC(NC)C(C)=O'],
             'Original products': ['CCCCO',
              'CNC(COP(=O)(O)OCC(C)(C)C(O)C(=O)NCCC(=O)NCCSC(=O)CCCCCC(=O)O)C(C)=O'],
             'Alternate products': ['CNC(COP(=O)(O)O)C(C)=O',
              'CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)CO'],
             'LHS_cofactors': ['O'],
             'RHS_cofactors': [],
             'Reaction eq': 'CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OCC(NC)C(C)=O + O = CNC(COP(=O)(O)O)C(C)=O + CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)CO',
             'Original rule': 'rule0007',
             'feasibility_label': 0,
             'remark': 'alternate product from pos rxn'} in neg_rxns_list

    assert {'Original enzymatic database': 'brenda',
 'Original reaction ID': '3.1.1.85_2',
 'Original substrates': ['CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OCC(NC)C(C)=O'],
 'Original products': ['CCCCO',
  'CNC(COP(=O)(O)OCC(C)(C)C(O)C(=O)NCCC(=O)NCCSC(=O)CCCCCC(=O)O)C(C)=O'],
 'Alternate products': ['CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)O',
  'CNC(CO)C(C)=O'],
 'LHS_cofactors': ['O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OCC(NC)C(C)=O + O = CCCCOC(=O)CCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)O + CNC(CO)C(C)=O',
 'Original rule': 'rule0007',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_substrate_and_two_products_rule7_rxn5():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = None
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == '3.4.17.1_43':
            pos_rxn_entry = rxn
            neg_rxns_list= neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = pos_rxn_entry,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                                        reported_products = pos_rxn_entry['Products'],
                                                                                        gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate no new negative reactions
    assert len(neg_rxns_list) == 0

def test_gen_neg_products_if_single_substrate_and_two_products_rule7_rxn6():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = None
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == '3.2.1.59_1':
            pos_rxn_entry = rxn
            neg_rxns_list= neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = pos_rxn_entry,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                                        reported_products = pos_rxn_entry['Products'],
                                                                                        gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 1 new negative reaction1
    assert len(neg_rxns_list) == 1

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    assert {'Original enzymatic database': 'brenda',
  'Original reaction ID': '3.2.1.59_1',
  'Original substrates': ['OCC1OC(OC2C(O)C(CO)OC(OC3C(O)C(CO)OC(OC4C(O)C(CO)OC(OC5C(O)C(O)OC(CO)C5O)C4O)C3O)C2O)C(O)C(O)C1O'],
  'Original products': ['OCC1OC(O)C(O)C(O)C1O',
   'OCC1OC(OC2C(O)C(CO)OC(OC3C(O)C(CO)OC(OC4C(O)C(O)OC(CO)C4O)C3O)C2O)C(O)C(O)C1O'],
  'Alternate products': ['OCC1OC(OC2C(O)C(CO)OC(OC3C(O)C(O)OC(CO)C3O)C2O)C(O)C(O)C1O',
   'OCC1OC(OC2C(O)C(O)OC(CO)C2O)C(O)C(O)C1O'],
  'LHS_cofactors': ['O'],
  'RHS_cofactors': [],
  'Reaction eq': 'OCC1OC(OC2C(O)C(CO)OC(OC3C(O)C(CO)OC(OC4C(O)C(CO)OC(OC5C(O)C(O)OC(CO)C5O)C4O)C3O)C2O)C(O)C(O)C1O + O = OCC1OC(OC2C(O)C(CO)OC(OC3C(O)C(O)OC(CO)C3O)C2O)C(O)C(O)C1O + OCC1OC(OC2C(O)C(O)OC(CO)C2O)C(O)C(O)C1O',
  'Original rule': 'rule0007',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'} in neg_rxns_list

def test_gen_neg_products_if_single_substrate_and_two_products_rule105_rxn7():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = None
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == 'RXN-13207_reverse':
            pos_rxn_entry = rxn
            neg_rxns_list= neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = pos_rxn_entry,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                                        reported_products = pos_rxn_entry['Products'],
                                                                                        gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 10 new negative reactions
    assert len(neg_rxns_list) == 10

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    # Finally, check the presence of alternate products and confirm that they are unique
    d1 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1C(O)OC(CO)C(O)C1O',
  'NC1CC(N)C(OP(=O)(O)OP(=O)(O)O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1C(O)OC(CO)C(O)C1O + NC1CC(N)C(OP(=O)(O)OP(=O)(O)O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d2 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O',
  'NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + NC1C(OP(=O)(O)OP(=O)(O)O)OC(CO)C(O)C1O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d3 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OP(=O)(O)OP(=O)(O)O)C1O',
  'O=P(O)(O)OCC1OC(O)C(O)C1O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OP(=O)(O)OP(=O)(O)O)C1O + O=P(O)(O)OCC1OC(O)C(O)C1O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d4 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OP(=O)(O)OP(=O)(O)O)C1O',
  'O=P(O)(O)OCC1OC(O)C(O)C1O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OP(=O)(O)OP(=O)(O)O)C1O + O=P(O)(O)OCC1OC(O)C(O)C1O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d5 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OP(=O)(O)OP(=O)(O)O)C1O',
  'O=P(O)(O)OCC1OC(O)C(O)C1O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OP(=O)(O)OP(=O)(O)O)C1O + O=P(O)(O)OCC1OC(O)C(O)C1O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d6 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OC2OC(CO)C(O)C(OP(=O)(O)OP(=O)(O)O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O',
  'O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OC2OC(CO)C(O)C(OP(=O)(O)OP(=O)(O)O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d7 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OC2OC(COP(=O)(O)OP(=O)(O)O)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O',
  'O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OC2OC(COP(=O)(O)OP(=O)(O)O)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d8 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OC2OC(CO)C(OP(=O)(O)OP(=O)(O)O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O',
  'O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OC2OC(CO)C(OP(=O)(O)OP(=O)(O)O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d9 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OP(=O)(O)OP(=O)(O)O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1OC1OC(CO)C(O)C(O)C1N',
  'O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OP(=O)(O)OP(=O)(O)O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1OC1OC(CO)C(O)C(O)C1N + O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d10 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-13207_reverse',
 'Original substrates': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O'],
 'Original products': ['NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(O)C1O',
  'O=P(O)(O)OCC1OC(OP(=O)(O)OP(=O)(O)O)C(O)C1O'],
 'Alternate products': ['NC1CC(N)C(OP(=O)(O)OP(=O)(O)O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1OC1OC(CO)C(O)C(O)C1N',
  'O'],
 'LHS_cofactors': ['O=P(O)(O)OP(=O)(O)O'],
 'RHS_cofactors': [],
 'Reaction eq': 'NC1CC(N)C(OC2OC(CO)C(O)C(O)C2N)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1O + O=P(O)(O)OP(=O)(O)O = NC1CC(N)C(OP(=O)(O)OP(=O)(O)O)C(OC2OC(COP(=O)(O)O)C(O)C2O)C1OC1OC(CO)C(O)C(O)C1N + O',
 'Original rule': 'rule0105',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}

    assert d1 in neg_rxns_list
    assert d2 in neg_rxns_list
    assert d3 in neg_rxns_list
    assert d4 in neg_rxns_list
    assert d5 in neg_rxns_list
    assert d6 in neg_rxns_list
    assert d7 in neg_rxns_list
    assert d8 in neg_rxns_list
    assert d9 in neg_rxns_list
    assert d10 in neg_rxns_list

    # Create a list of all alternate products from each alternate reaction dictionary
    all_alt_product_lists = [ dd['Alternate products'] for dd in [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10] ]

    all_alt_product_sets = []

    for alt_product_list in all_alt_product_lists:
        alt_product_set = set(alt_product_list)
        all_alt_product_sets.append(alt_product_set)

    assert len(all_alt_product_sets) == len(all_alt_product_lists)

def test_gen_neg_products_if_single_substrate_and_three_products_rule265_rxn8():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = None
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == 'RXN-7885_reverse':
            pos_rxn_entry = rxn
            neg_rxns_list= neg_data_utils.run_pickaxe_if_multiple_reported_products(pos_rxn_entry = pos_rxn_entry,
                                                                                        coreactants_filepath = coreactants_filepath,
                                                                                        reported_substrate = pos_rxn_entry['Substrates'][0],
                                                                                        reported_products = pos_rxn_entry['Products'],
                                                                                        gen_rules_filepath = gen_rules_filepath)

    # this reaction should generate 19 new negative reactions
    assert len(neg_rxns_list) == 19

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:

        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    d1 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=O',
  'CC(=O)C=O',
  'CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=O + CC(=O)C=O + CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d2 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=O',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=O + CC(=O)C=CC1=C(C)CC(O)CC1(C)C + CC(=O)C=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d3 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(=O)C=CC=C(C)C=O',
  'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=O',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(=O)C=CC=C(C)C=O + CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=O + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d4 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
  'CC(C=O)=CC=O',
  'CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O + CC(C=O)=CC=O + CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d5 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=O)=CC=CC=O',
  'CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=O)=CC=CC=O + CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d6 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=O)=CC=CC=C(C)C=O',
  'CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
  'CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=O)=CC=CC=C(C)C=O + CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O + CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d7 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=O)=CC=CC=C(C)C=CC=O',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=O)=CC=CC=C(C)C=CC=O + CC(=O)C=CC1=C(C)CC(O)CC1(C)C + CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d8 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=O)=CC=CC=C(C)C=CC=C(C)C=O',
  'CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=O)=CC=CC=C(C)C=CC=C(C)C=O + CC(C=CC1=C(C)CC(O)CC1(C)C)=CC=O + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d9 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=O)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'O=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(=O)C=CC1=C(C)CC(O)CC1(C)C + CC(C=O)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d10 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=O)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=O)=CC=O',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=O)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + CC(C=O)=CC=O + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d11 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
  'CC(C=CC=O)=CC=O',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O + CC(C=CC=O)=CC=O + CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d12 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=O)=CC=CC=O',
  'CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=O)=CC=CC=O + CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d13 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=O)=CC=CC(C)=CC=CC=C(C)C=CC=O',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=O)=CC=CC(C)=CC=CC=C(C)C=CC=O + CC(=O)C=CC1=C(C)CC(O)CC1(C)C + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d14 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=O)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=O',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=O)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + CC(=O)C=O + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d15 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
  'CC(C=O)=CC=CC(C)=CC=O',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O + CC(C=O)=CC=CC(C)=CC=O + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d16 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=O)=CC=CC(C)=CC=CC=O',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + CC(C=O)=CC=CC(C)=CC=CC=O + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d17 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=O)=CC=CC(C)=CC=CC=C(C)C=CC=C(C)C=O',
  'CC1=C(C=O)C(C)(C)CC(O)C1',
  'CC1=C(C=O)C(C)(C)CC(O)C1'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=O)=CC=CC(C)=CC=CC=C(C)C=CC=C(C)C=O + CC1=C(C=O)C(C)(C)CC(O)C1 + CC1=C(C=O)C(C)(C)CC(O)C1',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d18 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O',
  'CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'O=CC=O'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=O + CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=CC=O',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d19 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-7885_reverse',
 'Original substrates': ['CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'Original products': ['CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(C=CC=O)=CC=CC=C(C)C=CC=O'],
 'Alternate products': ['O=CC=CC=O',
  'CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
  'CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C'],
 'LHS_cofactors': ['O=O', 'O=O'],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C)=CC=CC=C(C)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + O=O + O=O = O=CC=CC=O + CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C + CC(=O)C=CC=C(C)C=CC1=C(C)CC(O)CC1(C)C',
 'Original rule': 'rule0265',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}

    assert d1 in neg_rxns_list
    assert d2 in neg_rxns_list
    assert d3 in neg_rxns_list
    assert d4 in neg_rxns_list
    assert d5 in neg_rxns_list
    assert d6 in neg_rxns_list
    assert d7 in neg_rxns_list
    assert d8 in neg_rxns_list
    assert d9 in neg_rxns_list
    assert d10 in neg_rxns_list
    assert d11 in neg_rxns_list
    assert d12 in neg_rxns_list
    assert d13 in neg_rxns_list
    assert d14 in neg_rxns_list
    assert d15 in neg_rxns_list
    assert d16 in neg_rxns_list
    assert d17 in neg_rxns_list
    assert d18 in neg_rxns_list
    assert d19 in neg_rxns_list

    # Create a list of all alternate products from each alternate reaction dictionary
    all_alt_product_lists = [ dd['Alternate products'] for dd in [d1,d2,d3,d4,d5,d6,d7,d8,d9,
                                                                  d10,d11,d12,d13,d14,d15,d16,d17,d18,d19] ]

    all_alt_product_sets = []

    for alt_product_list in all_alt_product_lists:
        alt_product_set = set(alt_product_list)
        all_alt_product_sets.append(alt_product_set)

    assert len(all_alt_product_sets) == len(all_alt_product_lists)

def test_gen_neg_products_if_multisubstrate_rxns_rule0001_rxn9():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = []
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == '3.4.21.62_48':

            pos_rxn_entry = rxn

            for substrate in rxn['Substrates']:
                neg_data_list = neg_data_utils.run_pickaxe_for_dimerization_rxns(
                    pos_rxn_entry = pos_rxn_entry,
                    coreactants_filepath = coreactants_filepath,
                    reported_substrate = substrate,
                    reported_products = rxn['Products'],
                    gen_rules_filepath = gen_rules_filepath)

                neg_rxns_list += neg_data_list

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:
        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(
            set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    d1 = {'Original enzymatic database': 'brenda',
  'Original reaction ID': '3.4.21.62_48',
  'Original substrates': ['CCCCO', 'CCOC(=O)C(C)NC(=O)OCc1ccccc1'],
  'Original products': ['CCO', 'CCCCOC(=O)C(C)NC(=O)OCc1ccccc1'],
  'Alternate products': ['CCCCOCCCC', 'O'],
  'LHS_cofactors': [],
  'RHS_cofactors': [],
  'Reaction eq': 'CCCCO + CCCCO = CCCCOCCCC + O',
  'Original rule': 'rule0001',
  'feasibility_label': 0,
  'remark': 'alternate product from pos rxn'}

    assert d1 in neg_rxns_list

    # confirm that only 1 negative product is generated
    assert len(neg_rxns_list) == 1

def test_gen_neg_products_if_multisubstrate_rxns_rule0001_rxn10():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = []
    pos_rxn_entry = None
    rxn_count = 0

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == 'RXN-11301':

            pos_rxn_entry = rxn
            rxn_count += 1

            for substrate in rxn['Substrates']:
                neg_data_list = neg_data_utils.run_pickaxe_for_dimerization_rxns(
                    pos_rxn_entry = pos_rxn_entry,
                    coreactants_filepath = coreactants_filepath,
                    reported_substrate = substrate,
                    reported_products = rxn['Products'],
                    gen_rules_filepath = gen_rules_filepath)

                neg_rxns_list += neg_data_list

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:
        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(
            set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    # check that for the reaction equation in each alternate reaction dictionary,
    # there are exactly two reactants (which must be identical under Pickaxe) and two products
    for neg_data_entry in neg_rxns_list:
        rxn_eq = neg_data_entry['Reaction eq']
        assert len(rxn_eq.split(' = ')[0].split(' + ')) == 2  # there should be 2 substrates
        assert len(set(rxn_eq.split(' = ')[0].split(' + '))) == 1  # but they are identical

        assert len(rxn_eq.split(' = ')[1].split(' + ')) == 2  # there should be 2 products per this rule

        assert neg_data_entry['Original rule'] == 'rule0001'  # all rules should be rule0126

    assert rxn_count == 1

    d1 = {'Original enzymatic database': 'metacyc',
            'Original reaction ID': 'RXN-11301',
            'Original substrates': [
                'CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O',
                'CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
            'Original products': [
                'CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O',
                'CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
            'Alternate products': [
                'CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C',
                'O=P(O)(O)OP(=O)(O)O'],
            'LHS_cofactors': [],
            'RHS_cofactors': [],
            'Reaction eq': 'CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O + CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O = CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C + O=P(O)(O)OP(=O)(O)O',
            'Original rule': 'rule0001',
            'feasibility_label': 0,
            'remark': 'alternate product from pos rxn'}
    d2 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-11301',
 'Original substrates': ['CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O',
  'CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
 'Original products': ['CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O',
  'CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
 'Alternate products': ['CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C',
  'O=P(O)(O)OP(=O)(O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O + CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O = CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C + O=P(O)(O)OP(=O)(O)O',
 'Original rule': 'rule0001',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d3 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-11301',
 'Original substrates': ['CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O',
  'CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
 'Original products': ['CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O',
  'CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
 'Alternate products': ['CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C',
  'O=P(O)(O)OP(=O)(O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O + CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O = CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C + O=P(O)(O)OP(=O)(O)O',
 'Original rule': 'rule0001',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d4 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': 'RXN-11301',
 'Original substrates': ['CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O',
  'CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
 'Original products': ['CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O',
  'CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
 'Alternate products': ['CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)OP(=O)(O)O',
  'O=P(O)(O)OP(=O)(O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O + CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O = CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)OP(=O)(O)O + O=P(O)(O)OP(=O)(O)O',
 'Original rule': 'rule0001',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d5 = {'Original enzymatic database': 'metacyc',
     'Original reaction ID': 'RXN-11301',
     'Original substrates': ['CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)(O)OP(=O)(O)O',
      'CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
     'Original products': ['CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O',
      'CC(=O)NC1C(OC2C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O'],
     'Alternate products': ['CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)OCC4OC(OC5C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C5OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C(O)C4OC4OC(CO)C(OC5OC(CO)C(O)C(O)C5NC(C)=O)C(OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C4NC(C)=O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O',
      'O'],
     'LHS_cofactors': [],
     'RHS_cofactors': [],
     'Reaction eq': 'CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O + CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O = CC(=O)NC1C(OC2C(CO)OC(OC3C(CO)OC(OC4C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C4OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)OCC4OC(OC5C(CO)OC(OP(=O)(O)OP(=O)(O)OCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)C(NC(C)=O)C5OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C(O)C4OC4OC(CO)C(OC5OC(CO)C(O)C(O)C5NC(C)=O)C(OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C4NC(C)=O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)C(NC(C)=O)C3O)C(NC(C)=O)C2OC(C)C(=O)NC(C)C(=O)NC(CCC(=O)NC(CCCC(N)C(=O)O)C(=O)NC(C)C(=O)NC(C)C(=O)O)C(=O)O)OC(CO)C(O)C1O + O',
     'Original rule': 'rule0001',
     'feasibility_label': 0,
     'remark': 'alternate product from pos rxn'}

    assert d1 in neg_rxns_list
    assert d2 in neg_rxns_list
    assert d3 in neg_rxns_list
    assert d4 in neg_rxns_list
    assert d5 in neg_rxns_list

    # confirm that 284 new negative reactions are generated
    assert len(neg_rxns_list) == 284

    # Create a list of all alternate products from each alternate reaction dictionary
    all_alt_product_lists = [ dd['Alternate products'] for dd in [d1, d2, d3, d4,d5] ]

    # Check that each alternate reaction dictionary has a unique set of products
    all_alt_product_sets = []

    for alt_product_list in all_alt_product_lists:
        alt_product_set = set(alt_product_list)
        all_alt_product_sets.append(alt_product_set)

    assert len(all_alt_product_sets) == len(all_alt_product_lists)

def test_gen_neg_products_if_multisubstrate_rxns_rule0126_rxn11():

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl', 'rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    neg_rxns_list = []
    pos_rxn_entry = None

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == '2.1.3.1-RXN':

            pos_rxn_entry = rxn

            for substrate in rxn['Substrates']:
                neg_data_list = neg_data_utils.run_pickaxe_for_dimerization_rxns(
                    pos_rxn_entry = pos_rxn_entry,
                    coreactants_filepath = coreactants_filepath,
                    reported_substrate = substrate,
                    reported_products = rxn['Products'],
                    gen_rules_filepath = gen_rules_filepath)

                neg_rxns_list += neg_data_list

    # check that original reaction metadata has not been changed
    for neg_data_entry in neg_rxns_list:
        assert neg_data_entry['Original enzymatic database'] == pos_rxn_entry['Enzymatic database']
        assert neg_data_entry['Original reaction ID'] == pos_rxn_entry['Reaction ID']
        assert neg_data_entry['Original substrates'] == pos_rxn_entry['Substrates']
        assert neg_data_entry['Original products'] == pos_rxn_entry['Products']
        assert neg_data_entry['Original rule'] == pos_rxn_entry['Rule']
        assert neg_data_entry['remark'] == 'alternate product from pos rxn'

        # but the alternate products cannot be the same as the original products
        assert neg_data_entry['Alternate products'] != pos_rxn_entry['Products']
        assert neg_data_entry['Alternate products'] != neg_data_entry['Original products']

        # to really double-check this, there should be no intersections between the alternate products and originals
        intersecting_list = list(
            set(neg_data_entry['Alternate products']) & set(neg_data_entry['Original products']))
        assert intersecting_list == [] or len(intersecting_list) < len(set(neg_data_entry['Alternate products']))

        # and the reaction equations of the negative and positive reactions certainly cannot be the same
        assert neg_data_entry['Reaction eq'] != pos_rxn_entry['Reaction eq']

        # finally, check that all synthetically generated reactions were labelled as infeasible
        assert neg_data_entry['feasibility_label'] == 0

    # check that for the reaction equation in each alternate reaction dictionary,
    # there are exactly two reactants (which must be identical under Pickaxe) and two products
    for neg_data_entry in neg_rxns_list:
        rxn_eq = neg_data_entry['Reaction eq']
        assert len(rxn_eq.split(' = ')[0].split(' + ')) == 2 # there should be 2 substrates
        assert len(set(rxn_eq.split(' = ')[0].split(' + '))) == 1 # but they are identical

        assert len(rxn_eq.split(' = ')[1].split(' + ')) == 2 # there should be 2 products per this rule

        assert neg_data_entry['Original rule'] == 'rule0126' # all rules should be rule0126

    # the first substrate, O=C(O)CC(=O)C(=O)O, in this reaction generated 6 negative reactions
    # 2nd substrate, CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O ...
    # ... in this reaction generated 319 negative reactions
    assert len(neg_rxns_list) == 6 + 319

    d1 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': '2.1.3.1-RXN',
 'Original substrates': ['O=C(O)CC(=O)C(=O)O',
  'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Original products': ['CC(=O)C(=O)O',
  'CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Alternate products': ['C',
  'CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCSC(=O)CCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O = C + CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCSC(=O)CCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O',
 'Original rule': 'rule0126',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d2 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': '2.1.3.1-RXN',
 'Original substrates': ['O=C(O)CC(=O)C(=O)O',
  'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Original products': ['CC(=O)C(=O)O',
  'CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Alternate products': ['O=CO', 'O=C(O)C(=O)CC(C(=O)O)C(=O)C(=O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'O=C(O)CC(=O)C(=O)O + O=C(O)CC(=O)C(=O)O = O=CO + O=C(O)C(=O)CC(C(=O)O)C(=O)C(=O)O',
 'Original rule': 'rule0126',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d3 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': '2.1.3.1-RXN',
 'Original substrates': ['O=C(O)CC(=O)C(=O)O',
  'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Original products': ['CC(=O)C(=O)O',
  'CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Alternate products': ['O=C(O)C(=O)C(C(=O)O)C(=O)O', 'CC(=O)C(=O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'O=C(O)CC(=O)C(=O)O + O=C(O)CC(=O)C(=O)O = O=C(O)C(=O)C(C(=O)O)C(=O)O + CC(=O)C(=O)O',
 'Original rule': 'rule0126',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d4 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': '2.1.3.1-RXN',
 'Original substrates': ['O=C(O)CC(=O)C(=O)O',
  'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Original products': ['CC(=O)C(=O)O',
  'CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Alternate products': ['O=C(O)C(=O)C(C(=O)O)C(=O)C(=O)O', 'CC(=O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'O=C(O)CC(=O)C(=O)O + O=C(O)CC(=O)C(=O)O = O=C(O)C(=O)C(C(=O)O)C(=O)C(=O)O + CC(=O)O',
 'Original rule': 'rule0126',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d5 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': '2.1.3.1-RXN',
 'Original substrates': ['O=C(O)CC(=O)C(=O)O',
  'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Original products': ['CC(=O)C(=O)O',
  'CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Alternate products': ['O=CC(=O)O', 'O=C(O)CC(C(=O)O)C(=O)C(=O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'O=C(O)CC(=O)C(=O)O + O=C(O)CC(=O)C(=O)O = O=CC(=O)O + O=C(O)CC(C(=O)O)C(=O)C(=O)O',
 'Original rule': 'rule0126',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d6 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': '2.1.3.1-RXN',
 'Original substrates': ['O=C(O)CC(=O)C(=O)O',
  'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Original products': ['CC(=O)C(=O)O',
  'CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Alternate products': ['O=CCC(=O)O', 'O=C(O)C(=O)C(C(=O)O)C(=O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'O=C(O)CC(=O)C(=O)O + O=C(O)CC(=O)C(=O)O = O=CCC(=O)O + O=C(O)C(=O)C(C(=O)O)C(=O)O',
 'Original rule': 'rule0126',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}
    d7 = {'Original enzymatic database': 'metacyc',
 'Original reaction ID': '2.1.3.1-RXN',
 'Original substrates': ['O=C(O)CC(=O)C(=O)O',
  'CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Original products': ['CC(=O)C(=O)O',
  'CC(C(=O)O)C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O'],
 'Alternate products': ['O=CO', 'O=C(O)CC(=O)C(C(=O)O)C(=O)C(=O)O'],
 'LHS_cofactors': [],
 'RHS_cofactors': [],
 'Reaction eq': 'O=C(O)CC(=O)C(=O)O + O=C(O)CC(=O)C(=O)O = O=CO + O=C(O)CC(=O)C(C(=O)O)C(=O)C(=O)O',
 'Original rule': 'rule0126',
 'feasibility_label': 0,
 'remark': 'alternate product from pos rxn'}

    assert d1 in neg_rxns_list
    assert d2 in neg_rxns_list
    assert d3 in neg_rxns_list
    assert d4 in neg_rxns_list
    assert d5 in neg_rxns_list
    assert d6 in neg_rxns_list
    assert d7 in neg_rxns_list


if __name__ == "__main__":

    with open('../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_10.pkl','rb') as file:
        processed_rxns = pickle.load(file)

    coreactants_filepath = '../data/raw/all_cofactors.tsv'
    gen_rules_filepath = '../data/raw/JN1224MIN_rules.tsv'

    for rxn in processed_rxns:
        if rxn['Reaction ID'] == '2.1.3.1-RXN':

            syn_gen_neg_data = neg_data_utils.run_pickaxe_for_dimerization_rxns(
                                                                            pos_rxn_entry = rxn,
                                                                            coreactants_filepath = coreactants_filepath,
                                                                            reported_substrate = rxn['Substrates'][0],
                                                                            reported_products = rxn['Products'],
                                                                            gen_rules_filepath = gen_rules_filepath)

            print(syn_gen_neg_data)