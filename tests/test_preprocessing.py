import pytest
import pandas as pd

try:
    from preprocessing import preprocessing
except ImportError:
    import sys
    sys.path.append('../../preprocessing')
    import preprocessing

# load in the data needed for tests
cofactors_filepath = '../../data/processed/expanded_cofactors_no_stereochem.tsv'
cofactors_df = pd.read_csv(cofactors_filepath, delimiter=',')
cofactors_list = list(pd.read_csv(cofactors_filepath, delimiter=',')['SMILES'])
rules_df = pd.read_csv('../../data/raw/JN1224MIN_rules.tsv', delimiter='\t')

# --------- Tests for 'load_cofactors' function in preprocessing package ---------
def test_load_cofactors():
    cofactors_df, cofactors_list = preprocessing.load_cofactors(cofactors_filepath = cofactors_filepath)
    # there are 43 unique cofactors
    # ATP serves as both a pyrophosphate donor (in which case it becomes AMP, which is then a pyrophosphate acceptor)
    # and a phosphate donor (in which case it becaomes ADP, which is then a phosphate acceptor)
    # PPI serves as both as pyrophosphate cofactor itself and a prenyl acceptor
    # looks like oxygen appears twice and is a bug in the original cofactors list also but it is not fatal
    assert len(cofactors_list) == 48
    for cofactor_smiles in cofactors_list:
        assert '@' not in cofactor_smiles # since the '@' symbol indicates chirality

    assert type(cofactors_df) == pd.DataFrame
    assert type(cofactors_list) == list

# --------- Tests for 'load_dataset' function in preprocessing package ---------
def test_load_dataset_1(db = 'brenda'):

    # ---- Test loading in BRENDA ----
    BRENDA_db = preprocessing.load_dataset(db = db)

    for key in BRENDA_db.keys():
        key = key.split(' ', 1)[0] # remove everything after space (e.g. "2.7.1.110 (transferred to EC 2.7.11.3)_0")
        key_prefix = key.split('_')[0]
        numbers_str = key_prefix.split('.') # since BRENDA's identifiers are EC numbers of the form i.j.k.l_xxx
        assert len(numbers_str) == 4 # check that all four EC numbers are present
        for num_str in numbers_str:
            assert num_str.isdigit() # confirm that i,j,k and l are valid digits
            num = int(num_str)
            assert num >= 0 # confirm that i,j,k, and l are integers

def test_load_dataset_2(db = 'kegg'):

    # ---- Test loading in KEGG ----
    KEGG_db = preprocessing.load_dataset(db = db)

    for key in KEGG_db.keys():
        assert key[0] == 'R' # KEGG keys are of the form 'Rxxxxx', e.g. 'R00031' or 'R00031_reverse'
        assert len(key.split('_')[0]) == 6 # Check that anything before the '_reverse' is of the form 'Rxxxxx'

def test_load_dataset_3(db = 'metacyc'):

    # ---- Test loading in MetaCyc ----
    MetaCyc_db = preprocessing.load_dataset(db = db)

    for key in MetaCyc_db.keys(): # MetaCyc keys often have a 'RXN' in them
        assert 'RXN' in key or key.count('.') == 0 # and keys without a 'RXN' should NOT have any EC numbers

def test_load_dataset_4(db = 'metacyc_V24'):
    # ---- Test loading in MetaCyc V24 ----
    MetaCyc_V24_db = preprocessing.load_dataset(db = db)

    for key in MetaCyc_V24_db.keys():
        assert 'RXN' in key or key.count('.') == 0

# --------- Tests for 'get_rule_info' in preprocessing package ---------
def test_get_rule_info_1(rules_df = rules_df, query_rule = 'rule0001'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 2
    assert num_products == 2
    assert lhs_cofactors == []
    assert rhs_cofactors == []

def test_get_rules_info_2(rules_df = rules_df, query_rule = 'rule0002'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['NAD_CoF']
    assert rhs_cofactors == ['NADH_CoF']

def test_get_rules_info_3(rules_df = rules_df, query_rule = 'rule0003'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['NADH_CoF']
    assert rhs_cofactors == ['NAD_CoF']

def test_get_rules_info_4(rules_df = rules_df, query_rule = 'rule0004'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['NADH_CoF','O2']
    assert rhs_cofactors == ['NAD_CoF','WATER']

def test_get_rules_info_5(rules_df = rules_df, query_rule = 'rule0005'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['NAD_CoF','WATER']
    assert rhs_cofactors == ['NADH_CoF','O2']

def test_get_rules_info_6(rules_df = rules_df, query_rule = 'rule0006'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 2
    assert num_products == 1
    assert lhs_cofactors == []
    assert rhs_cofactors == ['WATER']

def test_get_rules_info_7(rules_df = rules_df, query_rule = 'rule0007'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 2
    assert lhs_cofactors == ['WATER']
    assert rhs_cofactors == []

def test_get_rules_info_8(rules_df = rules_df, query_rule = 'rule0008'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['GLUCOSYL_ACCEPTOR_CoF']
    assert rhs_cofactors == ['GLUCOSYL_DONOR_CoF']

def test_get_rules_info_9(rules_df = rules_df, query_rule = 'rule0009'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['GLUCOSYL_DONOR_CoF']
    assert rhs_cofactors == ['GLUCOSYL_ACCEPTOR_CoF']

def test_get_rules_info_10(rules_df = rules_df, query_rule = 'rule0010'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['METHYL_ACCEPTOR_CoF']
    assert rhs_cofactors == ['METHYL_DONOR_CoF']

def test_get_rules_info_11(rules_df = rules_df, query_rule = 'rule0243'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df=rules_df,
                                                                                             query_rule=query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['WATER','WATER']
    assert rhs_cofactors == ['NH3']

def test_get_rules_info_12(rules_df = rules_df, query_rule = 'rule0269'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df=rules_df,
                                                                                             query_rule=query_rule)
    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['CoA', 'NAD_CoF']
    assert rhs_cofactors == ['NADH_CoF','CO2']

def test_get_rules_info_13(rules_df = rules_df, query_rule = 'rule0695'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)
    assert num_substrates == 2
    assert num_products == 1
    assert lhs_cofactors == ['NADH_CoF']
    assert rhs_cofactors == ['CoA', 'NAD_CoF', 'WATER']

def test_get_rules_info_14(rules_df = rules_df, query_rule = 'rule0932'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df = rules_df,
                                                                                             query_rule = query_rule)

    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['FADH2_CoF','CO2']
    assert rhs_cofactors == ['FAD_CoF']

def test_get_rules_info_15(rules_df = rules_df, query_rule = 'rule1224'):
    num_substrates, num_products, lhs_cofactors, rhs_cofactors = preprocessing.get_rule_info(rules_df=rules_df,
                                                                                             query_rule=query_rule)

    assert num_substrates == 1
    assert num_products == 1
    assert lhs_cofactors == ['WATER']
    assert rhs_cofactors == ['CO2','NH3']

# --------- Tests for 'lookup_cofactor_with_code' in preprocessing package ---------

# Check correct SMILES is returned when the cofactor ID in a reaction rule is 'PYROPHOSPHATE_DONOR_CoF
def test_lookup_cofactor_by_ID_1():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'PYROPHOSPHATE_DONOR_CoF',cofactors_df = cofactors_df)
    assert cofactor_smi == 'Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O'

# Check correct SMILES is returned when the cofactor ID in a reaction rule is 'FAD_CoF'
def test_lookup_cofactor_by_ID_2():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'FAD_CoF', cofactors_df = cofactors_df)
    assert cofactor_smi == \
           'Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)c2cc1C'

# Check correct SMILES is returned when the cofactor ID in a reaction rule is 'FADH2_CoF'
def test_lookup_cofactor_by_ID_3():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'FADH2_CoF', cofactors_df = cofactors_df)
    assert cofactor_smi == \
           'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1[nH]c(=O)[nH]c(=O)c1N2'

# Check correct SMILES is returned when the cofactor ID in a reaction rule is 'NAD_CoF' (should return NADP+ SMILES)
# NADP+ is the official cofactor used in our Pickaxe V2.0 expansions
def test_lookup_cofactor_by_ID_4():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'NAD_CoF', cofactors_df = cofactors_df)
    assert cofactor_smi == 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'

# Check correct SMILES is returned when the cofactor ID in a reaction rule is 'NADH_CoF' (should return NADPH SMILES)
# NADPH is the official cofactor used in our Pickaxe V2.0 expansions
def test_lookup_cofactor_by_ID_5():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'NADH_CoF', cofactors_df = cofactors_df)
    assert cofactor_smi == 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1'

# Check correct SMILES are returned when the cofactor ID is 'NAD' (should return NAD+ SMILES)
# this will actually never happen because NAD+ is NOT an official cofactor in Pickaxe V2.0
# and pickaxe always defaults to NADP+
# but since reactions recorded in databases actually have both NAD+/NADH as well as NADP+/NADPH,
# we need to find a way to allow that either cofactor is fine when present in a mined reaction
def test_lookup_cofactor_by_ID_6():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'NAD', cofactors_df = cofactors_df)
    assert cofactor_smi == 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1'

# Check correct SMILES are returned when the cofactor ID is 'NADH' (should return NADH SMILES)
# again this will never happen because NADH is NOT an official cofactor in Pickaxe V2.0 and Pickaxe defaults to NADPH
# so we again need to allow that either NADH or NADPH is a valid cofactor when present in a mined reaction
def test_lookup_cofactor_by_ID_7():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'NADH', cofactors_df = cofactors_df)
    assert cofactor_smi == 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1'

# Check correct SMILES are returned when the cofactor ID is 'WATER'
def test_lookup_cofactor_by_ID_8():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'WATER', cofactors_df = cofactors_df)
    assert cofactor_smi == 'O'

# Check correct SMILES are returned when the cofactor ID is 'NH3'
def test_lookup_cofactor_by_ID_9():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'NH3', cofactors_df = cofactors_df)
    assert cofactor_smi == 'N'

# Check correct SMILES are returned when the cofactor ID is 'CoA'
def test_lookup_cofactor_by_ID_10():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'CoA', cofactors_df = cofactors_df)
    assert cofactor_smi == 'CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS'

# Check correct SMILES are returned when the cofactor ID is 'CO2'
def test_lookup_cofactor_by_ID_11():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'CO2', cofactors_df = cofactors_df)
    assert cofactor_smi == 'O=C=O'

# Check correct SMILES are returned when the cofactor ID is 'O2'
def test_lookup_cofactor_by_ID_12():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID='O2', cofactors_df=cofactors_df)
    assert cofactor_smi == 'O=O'

# Check correct SMILES are returned when the cofactor ID is 'OXYGEN' (yes oxygen appears twice in list from JN1224MIN)
def test_lookup_cofactor_by_ID_13():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID='OXYGEN', cofactors_df=cofactors_df)
    assert cofactor_smi == 'O=O'

# Check correct SMILES are returned when the cofactor ID is 'PYROPHOSPHATE_DONOR_CoF' (this is ATP)
def test_lookup_cofactor_by_ID_14():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'PYROPHOSPHATE_DONOR_CoF', cofactors_df = cofactors_df)
    assert cofactor_smi == 'Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O'

# Check correct SMILES are returned when the cofactor ID is 'PHOSPHATE_DONOR_CoF' (this is also ATP)
def test_lookup_cofactor_by_ID_15():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'PHOSPHATE_DONOR_CoF', cofactors_df = cofactors_df)
    assert cofactor_smi == 'Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O'

# Check correct SMILES are returned when the cofactor ID is 'PRENYL_ACCEPTOR_CoF' (this is PPI)
def test_lookup_cofactor_by_ID_16():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'PRENYL_ACCEPTOR_CoF', cofactors_df = cofactors_df)
    assert cofactor_smi == 'O=P(O)(O)OP(=O)(O)O'

# Check correct SMILES are returned when the cofactor ID is 'PPI' (this of course is PPI too)
def test_lookup_cofactor_by_ID_17():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID = 'PPI', cofactors_df = cofactors_df)
    assert cofactor_smi == 'O=P(O)(O)OP(=O)(O)O'

def test_lookup_cofactor_by_ID_18():
    cofactor_smi = preprocessing.lookup_cofactor_by_ID(ID='AMINO_CoF', cofactors_df=cofactors_df)
    assert cofactor_smi == 'NC(CCC(=O)O)C(=O)O'

# --------- Tests for 'cofactor_ids_to_smiles' in preprocessing package ---------
def test_cofactor_ids_to_smiles_1():
    cofactor_ids_list = ['NAD_CoF']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_smiles_list) == 2 # the cofactor SMILES list should contain NAD+ and NADP+ SMILES
    assert 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1' in cofactor_smiles_list
    assert 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1' in cofactor_smiles_list

def test_cofactor_ids_to_smiles_2():
    cofactor_ids_list = ['NADH_CoF']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_smiles_list) == 2 # the cofactor SMILES list should contain NADH and NADPH SMILES
    assert 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1' in cofactor_smiles_list
    assert 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1' in cofactor_smiles_list

def test_cofactor_ids_to_smiles_3():
    cofactor_ids_list = ['WATER']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_smiles_list) == 1
    assert 'O' in cofactor_smiles_list

def test_cofactor_ids_to_smiles_4():
    cofactor_ids_list = []
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_smiles_list) == 0
    assert cofactor_smiles_list == []

def test_cofactor_ids_to_smiles_5():
    cofactor_ids_list = ['NH3']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_smiles_list) == 1
    assert cofactor_smiles_list == ['N']

def test_cofactor_ids_to_smiles_6():
    cofactor_ids_list = ['WATER','WATER']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_ids_list) == 2
    assert cofactor_smiles_list == ['O','O']

def test_cofactor_ids_to_smiles_7():
    cofactor_ids_list = ['CoA', 'NAD_CoF']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_smiles_list) == 3

    assert 'CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS' \
           in cofactor_smiles_list # CoA should be present (note there is a difference between CoA and Acetyl-CoA)

    assert 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1'\
           in cofactor_smiles_list # NADP+ should be present

    assert 'NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1' in\
           cofactor_smiles_list # NAD+ should be present

def test_cofactor_ids_to_smiles_8():
    cofactor_ids_list = ['NADH_CoF','CO2']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)

    assert len(cofactor_smiles_list) == 3

    assert 'O=C=O' in cofactor_smiles_list # CO2 should be present

    assert 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1' in \
           cofactor_smiles_list # NADH should be present

    assert 'NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1' in \
           cofactor_smiles_list # NADPH should be present

def test_cofactor_ids_to_smiles_9():
    cofactor_ids_list = ['FADH2_CoF', 'CO2']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list = cofactor_ids_list,
                                                                cofactors_df = cofactors_df)
    assert len(cofactor_smiles_list) == 2
    assert 'O=C=O' in cofactor_smiles_list
    assert 'Cc1cc2c(cc1C)N(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC1OC(n3cnc4c(N)ncnc43)C(O)C1O)c1[nH]c(=O)[nH]c(=O)c1N2' in cofactor_smiles_list

def test_cofactor_ids_to_smiles_10():
    cofactor_ids_list = ['FAD_CoF']
    cofactor_smiles_list = preprocessing.cofactor_ids_to_smiles(cofactor_ids_list=cofactor_ids_list,
                                                                cofactors_df=cofactors_df)
    assert len(cofactor_smiles_list) == 1
    assert 'Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(CC(O)C(O)C(O)COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)c2cc1C' in cofactor_smiles_list

# --------- Tests for 'fix_cofactors' in preprocessing package ---------

# first, test all the SMILES strings in the cofactors list have no stereochemistry
def test_cofactors_list_in_processed_data():
    # there are 43 unique cofactors
    # ATP serves as both a pyrophosphate donor (in which case it becomes AMP, which is then a pyrophosphate acceptor)
    # and a phosphate donor (in which case it becaomes ADP, which is then a phosphate acceptor)
    # PPI serves as both as pyrophosphate cofactor itself and a prenyl acceptor
    # looks like oxygen appears twice and is a bug in the original cofactors list also but it is not fatal
    assert len(cofactors_list) == 48
    for cofactor_smiles in cofactors_list:
        assert '@' not in cofactor_smiles # since the '@' symbol indicates chirality

# test fixing incorrect NAD+ SMILES string that is missing an atom to correct the NAD+ SMILES
nad_plus_incorrect_smiles = "*OC1C(O)C(COP(=O)(O)OP(=O)(O)OCC2OC([n+]3cccc(C(N)=O)c3)C(O)C2O)OC1n1cnc2c(N)ncnc21"
def test_fix_cofactors_1():
    nad_plus_correct_smiles = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = nad_plus_incorrect_smiles,
                                       cofactors_list = cofactors_list) == nad_plus_correct_smiles

# test fixing incorrect NADH SMILES string that is missing an atom to the correct NADH SMILES
nadh_incorrect_smiles = "*OC1C(O)C(COP(=O)(O)OP(=O)(O)OCC2OC(N3C=CCC(C(N)=O)=C3)C(O)C2O)OC1n1cnc2c(N)ncnc21"
def test_fix_cofactors_2():
    nadh_correct_smiles = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = nadh_incorrect_smiles,
                                       cofactors_list = cofactors_list) == nadh_correct_smiles

# test replacing glucosyl acceptor not in our cofactors list with UDP (our official glucosyl acceptor)
GDP_smiles = "Nc1nc(=O)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)O)C(O)C3O)c2[nH]1"
def test_fix_cofactors_3():
    UDP_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = GDP_smiles,
                                       cofactors_list = cofactors_list) == UDP_smiles

# test replacing glucosyl acceptor not in our cofactors list with UDP (our official glucosyl acceptor)
CDP_smiles = "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)n1"
def test_fix_cofactors_4():
    UDP_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = CDP_smiles,
                                       cofactors_list = cofactors_list) == UDP_smiles

# test replacing glucosyl acceptor not in our cofactors list with UDP (our official glucosyl acceptor)
NDP_smiles = "CC1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O"
def test_fix_cofactors_5():
    UDP_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = NDP_smiles,
                                       cofactors_list = cofactors_list) == UDP_smiles

# test replacing glucosyl acceptor not in our cofactors list with UDP (our official glucosyl acceptor)
dTDP_smiles = "Cc1cn(C2CC(O)C(COP(=O)(O)OP(=O)(O)O)O2)c(=O)[nH]c1=O"
def test_fix_Cofactors_6():
    UDP_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = dTDP_smiles,
                                       cofactors_list = cofactors_list) == UDP_smiles

# test replacing glucosyl donor not in our cofactors list with CPD_1275 (our official glucosyl donor)
GDP_mannose_smiles = "Nc1nc(=O)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)OC4OC(CO)C(O)C(O)C4O)C(O)C3O)c2[nH]1"
def test_fix_cofactors_7():
    CPD_12575_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = GDP_mannose_smiles,
                                       cofactors_list = cofactors_list) == CPD_12575_smiles

# test replacing glucosyl donor not in our cofactors list with CPD_1275 (our official glucosyl donor)
CDP_alpha_D_glucose_smiles = "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)n1"
def test_fix_cofactors_8():
    CPD_12575_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = CDP_alpha_D_glucose_smiles,
                                       cofactors_list = cofactors_list) == CPD_12575_smiles

# test replacing glucosyl donor not in our cofactors list with CPD_1275 (our official glucosyl donor)
NDP_glucose_smiles = "CC1OC(COP(=O)(O)OP(=O)(O)OC2OC(CO)C(O)C(O)C2O)C(O)C1O"
def test_fix_cofactors_9():
    CPD_12575_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = NDP_glucose_smiles,
                                       cofactors_list = cofactors_list) == CPD_12575_smiles

# test replacing glucosyl donor not in our cofactors list with CPD_1275 (our official glucosyl donor)
dTDP_alpha_D_glucose_smiles = "Cc1cn(C2CC(O)C(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)O2)c(=O)[nH]c1=O"
def test_fix_cofactors_10():
    CPD_12575_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = dTDP_alpha_D_glucose_smiles,
                                       cofactors_list = cofactors_list) == CPD_12575_smiles

# test fixing a phosphate acceptor that is missing an atom and also not in our cofactor list with ADP
# where ADP is our official phosphate acceptor
phosphate_acceptor_incorrect = "*C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O"
def test_fix_factors_11():
    # phosphate acceptor that is officially in our cofactors list
    ADP_smiles = "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = phosphate_acceptor_incorrect,
                                       cofactors_list = cofactors_list) == ADP_smiles

# test fixing phosphate donor that is missing an atom and also not in our cofactor list with ATP
# where ATP is our official phosphate donor
phosphate_donor_incorrect_1 = "*C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O" # missing 1 atom
def test_fix_cofactors_12():
    # phosphate donor that is officially in our cofactors list
    ATP_smiles = "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles=phosphate_donor_incorrect_1,
                                       cofactors_list = cofactors_list) == ATP_smiles

# test fixing phosphate donor that is missing two atoms and also not in our cofactor list with ATP
# where ATP is our official phosphate donor
phosphate_donor_incorrect_2 = "*C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1*"
def test_fix_cofactors_13():
    # phosphate donor that is officially in our cofactors list
    ATP_smiles = "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O"
    assert preprocessing.fix_cofactors(mined_cofactor_smiles=phosphate_donor_incorrect_2,
                                       cofactors_list = cofactors_list) == ATP_smiles

# test retaining the SMILES string of all cofactors that already exist in our cofactors list
def test_fix_cofactors_14():
    for cofactor_smiles in cofactors_list:
        assert preprocessing.fix_cofactors(mined_cofactor_smiles = cofactor_smiles,
                                           cofactors_list = cofactors_list) == cofactor_smiles

# test returning None when a cofactor that is not in our list or not accounted for is provided
random_non_cofactor_1 = 'CC'
def test_fix_cofactors_15():
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = random_non_cofactor_1,
                                       cofactors_list = cofactors_list) is None

# test returning None when a cofactor that is not in our list or not accounted for is provided
random_non_cofactor_2 = 'CCC'
def test_fix_cofactors_16():
    assert preprocessing.fix_cofactors(mined_cofactor_smiles = random_non_cofactor_2,
                                       cofactors_list = cofactors_list) is None

# --------- Tests for 'process_reactants_dict' in preprocessing package ---------
def test_process_reactants_dict_1():
    reactants_dict = {"METHYLAMINE:0": "CN",
                      "NADPH:0": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1",
                      "PYRUVATE:0": "CC(=O)C(=O)O"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1"]
    assert substrates_rxn == ["CN","CC(=O)C(=O)O"]

def test_process_reactants_dict_2():
    reactants_dict = {"ALPHA-N-ACETYLNEURAMINYL-ETCETERA:0": "*O[C@H]1O[C@H](CO)[C@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@]3(C(=O)O)C[C@H](O)[C@@H](NC(C)=O)[C@H]([C@H](O)[C@H](O)CO)O3)[C@H]2O)[C@H]1NC(C)=O",
                       "CMP-N-ACETYL-NEURAMINATE:0": "CC(=O)N[C@@H]1[C@@H](O)C[C@@](OP(=O)(O)OC[C@H]2O[C@@H](n3ccc(N)nc3=O)[C@H](O)[C@@H]2O)(C(=O)O)O[C@H]1[C@H](O)[C@H](O)CO"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == []
    assert substrates_rxn == ["*OC1OC(CO)C(O)C(OC2OC(CO)C(O)C(OC3(C(=O)O)CC(O)C(NC(C)=O)C(C(O)C(O)CO)O3)C2O)C1NC(C)=O",
                              "CC(=O)NC1C(O)CC(OP(=O)(O)OCC2OC(n3ccc(N)nc3=O)C(O)C2O)(C(=O)O)OC1C(O)C(O)CO"]

def test_process_reactants_dict_3():
    reactants_dict = {"BENZYL-ALCOHOL:0": "OCc1ccccc1",
                      "CARBON-DIOXIDE:0": "O=C=O",
                      "LEU:0": "CC(C)C[C@H](N)C(=O)O"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["O=C=O"]
    assert len(substrates_rxn) == 2
    assert "OCc1ccccc1" in substrates_rxn
    assert "CC(C)CC(N)C(=O)O" in substrates_rxn

def test_process_reactants_dict_4():
    reactants_dict = {"CPD-16854:0": "CC(=O)NC1C(=O)C(O)C(CO)OC1OP(=O)(O)OP(=O)(O)OCC1OC(n2ccc(=O)[nH]c2=O)C(O)C1O",
                      "NADH:0": "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]
    assert substrates_rxn == ["CC(=O)NC1C(=O)C(O)C(CO)OC1OP(=O)(O)OP(=O)(O)OCC1OC(n2ccc(=O)[nH]c2=O)C(O)C1O"]

def test_process_reactants_dict_5():
    reactants_dict = {"CPD-10791:0": "CC(=O)C(O)C(O)C(=O)COP(=O)(O)O", "GAP:0": "O=CC(O)COP(=O)(O)O"}
    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == []
    assert substrates_rxn == ["CC(=O)C(O)C(O)C(=O)COP(=O)(O)O","O=CC(O)COP(=O)(O)O"]

def test_process_reactants_dict_6():
    reactants_dict = {"CPD-13090:0": "CCC=C\\CC1OC1C/C=C\\CCCCCCCC(=O)O", "WATER:0": "O"}
    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["O"]
    assert substrates_rxn == ["CCC=CCC1OC1CC=CCCCCCCCC(=O)O"]

def test_process_reactants_dict_7():
    reactants_dict = {
        "CPDQT-295:0": "CS(=O)CCCCCCC/C(=N\\OS(=O)(=O)O)S[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
        "NADP:0": "NC(=O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1",
        "WATER:0": "O"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)
    assert lhs_cofactors_rxn == \
           ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1",
            "O"]
    assert substrates_rxn == \
           ["CS(=O)CCCCCCCC(=NOS(=O)(=O)O)SC1OC(CO)C(O)C(O)C1O"]

def test_process_reactants_dict_8():
    reactants_dict = {"CPD-11975:0": "CC(=O)N[C@H]1[C@@H](O[C@@H]2[C@H](O)[C@H](OP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]2O)O[C@H](CO)[C@@H](O)[C@@H]1O",
                      "UDP:0": "O=c1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"] # UDP
    assert substrates_rxn == ["CC(=O)NC1C(OC2C(O)C(O)C(O)C(OP(=O)(O)O)C2O)OC(CO)C(O)C1O"]

def test_process_reactants_dict_9():
    reactants_dict = {"CPD-12575:0": "O=c1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H]2O)c(=O)[nH]1",
                      "CPD-8665:0": "CC(/C=C/C=C(\\C)C(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)C(=O)O[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1"]  # CPD12575 (glucosyl donor)
    assert substrates_rxn == ["CC(C=CC=C(C)C(=O)OC1OC(CO)C(O)C(O)C1O)=CC=CC=C(C)C=CC=C(C)C(=O)OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O"]

def test_process_reactants_dict_10():
    reactants_dict = {"LONG-CHAIN-KETONE:0": "*CC(=O)C*",
                      "Reduced-Factor-F420:0": "*C(=O)[C@H](C)OP(=O)(O)OC[C@@H](O)[C@@H](O)[C@@H](O)CN1c2cc(O)ccc2Cc2c1[nH]c(=O)[nH]c2=O"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict=reactants_dict,
                                                                             cofactors_list=cofactors_list)

    assert lhs_cofactors_rxn == ["*C(=O)C(C)OP(=O)(O)OCC(O)C(O)C(O)CN1c2cc(O)ccc2Cc2c1[nH]c(=O)[nH]c2=O"]
    assert substrates_rxn == ["*CC(=O)C*"]

def test_process_reactants_dict_11():
    reactants_dict = {"(-)-perillyl alcohol:0": "C=C(C)C1CC=C(CO)CC1",
                      "NAD+:0": "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]
    assert substrates_rxn == ["C=C(C)C1CC=C(CO)CC1"]

def test_process_reactants_dict_12():
    reactants_dict = {"NADPH:0":
                          "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1",
                      "acetyl-CoA:0":
                          "CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1",
                              "CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O"]

    assert substrates_rxn == []

def test_process_reactants_dict_13():
    reactants_dict = {"2-oxo-3-hexenedioate:0": "O=C(O)CC=CC(=O)C(=O)O",
                      "ammonia:0": "N"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == ["N"]
    assert substrates_rxn == ["O=C(O)CC=CC(=O)C(=O)O"]

def test_process_reactants_dict_14():
    reactants_dict = {"cpd00228:0": "O[C@@H](CS)[C@@H](O)CS",
                      "cpd03477:0": "C/C(=C\\CC12OC1(C)C(=O)c1ccccc1C2=O)CCC[C@H](C)CCC[C@H](C)CCCC(C)C"}

    lhs_cofactors_rxn, substrates_rxn = preprocessing.process_reactants_dict(reactants_dict = reactants_dict,
                                                                             cofactors_list = cofactors_list)

    assert lhs_cofactors_rxn == []
    assert substrates_rxn == ["OC(CS)C(O)CS",
                              "CC(=CCC12OC1(C)C(=O)c1ccccc1C2=O)CCCC(C)CCCC(C)CCCC(C)C"]

# --------- Tests for 'process_products_dict' in preprocessing package ---------
def test_process_products_dict_1():
    products_dict = {
        "CPD-13667:0": "CC1(C)CC[C@]2(C)CC[C@]3(C)C(=CC(=O)C4[C@@]5(C)CC[C@H](O)C(C)(C)C5CC[C@]43C)C2C1",
        "NADPH:0": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1",
        "NADPH:1": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1",
        "NADPH:2": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1",
        "OXYGEN-MOLECULE:0": "O=O",
        "OXYGEN-MOLECULE:1": "O=O",
        "OXYGEN-MOLECULE:2": "O=O"}

    rhs_cofactors_list, products_list = preprocessing.process_products_dict(products_dict = products_dict,
                                                                            cofactors_list = cofactors_list)

    assert rhs_cofactors_list == [
        "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1",
        "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1",
        "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1",
        "O=O",
        "O=O",
        "O=O"]

    assert products_list == ["CC1(C)CCC2(C)CCC3(C)C(=CC(=O)C4C5(C)CCC(O)C(C)(C)C5CCC43C)C2C1"]

def test_process_products_dict_2():
    products_dict = {"CPD-199:0": "CC(C)C[C@H](NC(=O)OCc1ccccc1)C(=O)O",
                     "WATER:0": "O"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["O"]
    assert products_rxn == ["CC(C)CC(NC(=O)OCc1ccccc1)C(=O)O"]

def test_process_products_dict_3():
    products_dict = {"NAD:0": "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1",
                     "UDP-N-ACETYL-D-GLUCOSAMINE:0": "CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]
    assert products_rxn == ["CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O"]

def test_process_products_dict_4():
    products_dict = {"FRUCTOSE-16-DIPHOSPHATE:0": "O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1O",
                     "METHYL-GLYOXAL:0": "CC(=O)C=O"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == []
    assert products_rxn == ["O=P(O)(O)OCC1OC(O)(COP(=O)(O)O)C(O)C1O","CC(=O)C=O"]

def test_process_products_dict_5():
    products_dict = {"CPD-13092:0": "CC/C=C\\CC(O)C(O)C/C=C\\CCCCCCCC(=O)O"} # the slashes indicate E/Z stereochemistry

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == []
    assert products_rxn == ["CCC=CCC(O)C(O)CC=CCCCCCCCC(=O)O"]

def test_process_products_dict_6():
    products_dict = {"CPDQT-296:0": "CSCCCCCCC/C(=N\\OS(=O)(=O)O)S[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
                     "NADPH:0": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1",
                     "OXYGEN-MOLECULE:0": "O=O"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1",
                                 "O=O"]
    assert products_rxn == ["CSCCCCCCCC(=NOS(=O)(=O)O)SC1OC(CO)C(O)C(O)C1O"]

def test_process_products_dict_7():
    products_dict = {"1-L-MYO-INOSITOL-1-P:0": "O=P(O)(O)O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
                     "UDP-N-ACETYL-D-GLUCOSAMINE:0": "CC(=O)N[C@H]1[C@@H](OP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H](n3ccc(=O)[nH]c3=O)[C@H](O)[C@@H]2O)O[C@H](CO)[C@@H](O)[C@@H]1O"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == []
    assert products_rxn == ["O=P(O)(O)OC1C(O)C(O)C(O)C(O)C1O",
                            "CC(=O)NC1C(OP(=O)(O)OP(=O)(O)OCC2OC(n3ccc(=O)[nH]c3=O)C(O)C2O)OC(CO)C(O)C1O"]

def test_process_products_dict_8():
    products_dict = {"CPD-8666:0": "CC(/C=C/C=C(\\C)C(=O)O[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)C(=O)O[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O",
                     "UDP:0": "O=c1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"]
    assert products_rxn == ["CC(C=CC=C(C)C(=O)OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O)=CC=CC=C(C)C=CC=C(C)C(=O)OC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O"]

def test_process_products_dict_9():
    products_dict = {"Oxidized-Factor-F420:0": "*C(=O)[C@H](C)OP(=O)(O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2nc(=O)[nH]c(=O)c-2cc2ccc(O)cc21",
                      "Secondary-Alcohols:0": "*CC(O)C*"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list =  cofactors_list)

    assert rhs_cofactors_rxn == ["*C(=O)C(C)OP(=O)(O)OCC(O)C(O)C(O)Cn1c2nc(=O)[nH]c(=O)c-2cc2ccc(O)cc21"]
    assert products_rxn == ["*CC(O)C*"]

def test_process_products_dict_10():
    products_dict = {"ACETYL-COA:0": "CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O",
                     "CPD-7221:0": "CCCCCCCC/C=C\\CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O",
                     "HYDROGEN-PEROXIDE:0": "OO"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O",
                                 "OO"]

    assert products_rxn == ["CCCCCCCCC=CCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O"]

def test_process_products_dict_11():
    products_dict = {"(S)-perillaldehyde:0": "C=C(C)C1CC=C(C=O)CC1",
                     "NADH:0": "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]

    assert products_rxn == ["C=C(C)C1CC=C(C=O)CC1"]

def test_process_products_dict_12():
    products_dict = {"NADP+:0": "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1",
                     "acetaldehyde:0": "CC=O",
                     "coenzyme A:0": "CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1",
                                 "CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS"]

    assert products_rxn == ["CC=O"]

def test_process_products_dict_13():
    products_dict = {"2-aminomuconate:0": "NC(=CC=CC(=O)O)C(=O)O",
                     "water:0": "O"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == ["O"]
    assert products_rxn == ["NC(=CC=CC(=O)O)C(=O)O"]

def test_process_products_dict_14():
    products_dict = {"cpd00823:0": "O[C@H]1CSSC[C@@H]1O",
                     "cpd01798:0": "C/C(=C\\CC1(O)C(=O)c2ccccc2C(=O)C1C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C"}

    rhs_cofactors_rxn, products_rxn = preprocessing.process_products_dict(products_dict = products_dict,
                                                                          cofactors_list = cofactors_list)

    assert rhs_cofactors_rxn == []
    assert products_rxn == ["OC1CSSCC1O",
                            "CC(=CCC1(O)C(=O)c2ccccc2C(=O)C1C)CCCC(C)CCCC(C)CCCC(C)C"]

# --------- Tests for 'reorder_smiles_by_mol_weight' in preprocessing package ---------
def test_reorder_smiles_by_mol_weight_1():
    smiles_list = ['O','N','CC(=O)C']
    descending = False
    assert preprocessing.reorder_smiles_by_mol_weight(smiles_list = smiles_list, descending = descending) == \
           ['N','O','CC(=O)C']

def test_reorder_smiles_by_mol_weight_2():
    smiles_list = ['O','N','CC(=O)C']
    descending = True
    assert preprocessing.reorder_smiles_by_mol_weight(smiles_list = smiles_list, descending = descending) == \
           ['CC(=O)C','O','N']

def test_reorder_smiles_by_mol_weight_3():
    smiles_list = ['c1ccccc1','C=C','CC(=O)C']
    descending = False
    assert preprocessing.reorder_smiles_by_mol_weight(smiles_list = smiles_list, descending = descending) == \
           ['C=C','CC(=O)C','c1ccccc1']

def test_reorder_smiles_by_mol_weight_4():
    smiles_list = ['c1ccccc1','C=C','CC(=O)C']
    descending = True
    assert preprocessing.reorder_smiles_by_mol_weight(smiles_list = smiles_list, descending = descending) == \
           ['c1ccccc1','CC(=O)C','C=C']

def test_reorder_smiles_by_mol_weight_5():
    # [Taxol, vitamin D3, adenine]
    smiles_list = ['CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C',
                   'CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C',
                   'C1=NC2=NC=NC(=C2N1)N']
    descending = False
    assert preprocessing.reorder_smiles_by_mol_weight(smiles_list = smiles_list, descending = descending) == \
           ['C1=NC2=NC=NC(=C2N1)N',
            'CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C',
            'CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C']

def test_reorder_smiles_by_mol_weight_6():
    smiles_list = [
        'CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C',
        'CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C',
        'C1=NC2=NC=NC(=C2N1)N']
    descending = True
    assert preprocessing.reorder_smiles_by_mol_weight(smiles_list=smiles_list, descending=descending) == \
           ['CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C',
            'CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C',
            'C1=NC2=NC=NC(=C2N1)N']

# check that an error is raised when SMILES strings with stereochemistry information still present are used
def test_reorder_smiles_by_mol_weight_7():
    stereo_smiles_list = ['CC[C@H](C)Br', 'C[C@@H](O)CC', 'CC[C@@H](C)[C@H](C)O']
    descending = False
    with pytest.raises(Exception) as exc_info:
        preprocessing.reorder_smiles_by_mol_weight(smiles_list = stereo_smiles_list, descending = descending)
    assert exc_info.type == ValueError
    assert "Stereochemistry encountered. Ensure SMILES are canonicalized and have no stereochemistry" in str(
                                                                                                        exc_info.value)

# --------- Tests for 'construct_rxn_str' in preprocessing package ---------
def test_construct_rxn_str_1():
    substrates = ["*C1CCC2C3C(=O)CC4CCCCC4(C)C3CCC12C"]
    lhs_cofactors = ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1"]
    products = ["*C1CCC2C3C(O)CC4CCCCC4(C)C3CCC12C"]
    rhs_cofactors = ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                            lhs_cofactors = lhs_cofactors,
                                            products = products,
                                            rhs_cofactors = rhs_cofactors)

    assert rxn_eq == \
           "*C1CCC2C3C(=O)CC4CCCCC4(C)C3CCC12C + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1 = " \
           "*C1CCC2C3C(O)CC4CCCCC4(C)C3CCC12C + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1"

def test_construct_rxn_str_2():
    substrates = ["CCC=CCC=CCC=CCC=CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O"]
    lhs_cofactors = ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]
    products = ["CCC=CCC=CCC=CCC=CCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O"]
    rhs_cofactors = ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "CCC=CCC=CCC=CCC=CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = CCC=CCC=CCC=CCC=CCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"

def test_construct_rxn_str_3():
    substrates = ["O=CC=Cc1ccc(C(=O)O)oc1=O"]
    lhs_cofactors = ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]
    products = ["O=C(O)c1ccc(C=CCO)c(=O)o1"]
    rhs_cofactors = ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "O=CC=Cc1ccc(C(=O)O)oc1=O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = O=C(O)c1ccc(C=CCO)c(=O)o1 + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"

def test_construct_rxn_str_4():
    substrates = ["CSCCCCCC(C(=O)O)C(=O)C(=O)O"]
    lhs_cofactors = ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]
    products = ["CSCCCCCC(C(=O)O)C(O)C(=O)O"]
    rhs_cofactors = ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "CSCCCCCC(C(=O)O)C(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = CSCCCCCC(C(=O)O)C(O)C(=O)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"

def test_construct_rxn_str_5():
    substrates = ["CSCCCCCC(C(=O)O)C(=O)C(=O)O"]
    lhs_cofactors = ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]
    products = ["CSCCCCCC(C(=O)O)C(O)C(=O)O"]
    rhs_cofactors = ["NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "CSCCCCCC(C(=O)O)C(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1 = CSCCCCCC(C(=O)O)C(O)C(=O)O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"

def test_construct_rxn_str_6():
    substrates = ["Oc1ccccc1"]
    lhs_cofactors = ["O=O","NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1"]
    products = ["Oc1ccccc1O"]
    rhs_cofactors = ["O","NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == f"{substrates[0]} + {lhs_cofactors[0]} + {lhs_cofactors[1]} = {products[0]} + {rhs_cofactors[0]} + {rhs_cofactors[1]}"

def test_construct_rxn_str_7():
    substrates = ["OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O"]
    lhs_cofactors = ["O"]
    products = ["OCC1OC(O)C(O)C(O)C1O",
                "OCC1OC(O)C(O)C(O)C1O"]
    rhs_cofactors = []

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O + O = OCC1OC(O)C(O)C(O)C1O + OCC1OC(O)C(O)C(O)C1O"

def test_construct_rxn_str_8():
    substrates = ["*C(=O)OC(CO)COP(=O)(O)OCC[N+](C)(C)C"]
    lhs_cofactors = ["O"]
    products = ["*C(=O)O",
                "C[N+](C)(C)CCOP(=O)(O)OCC(O)CO"]
    rhs_cofactors = []

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "*C(=O)OC(CO)COP(=O)(O)OCC[N+](C)(C)C + O = *C(=O)O + C[N+](C)(C)CCOP(=O)(O)OCC(O)CO"

def test_construct_rxn_str_9():
    substrates = ["Oc1ccc(-c2[o+]c3cc(O)cc(O)c3cc2O)cc1"]
    lhs_cofactors = ["O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1"]
    products = ["OCC1OC(Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)cc2)C(O)C(O)C1O"]
    rhs_cofactors = ["O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "Oc1ccc(-c2[o+]c3cc(O)cc(O)c3cc2O)cc1 + O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1 = OCC1OC(Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)cc2)C(O)C(O)C1O + O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"

def test_construct_rxn_str_10():
    substrates = ["NCCCCCC(=O)NCCCCCC(=O)O"]
    lhs_cofactors = ["O"]
    products = ["NCCCCCC(=O)O",
                "NCCCCCC(=O)O"]
    rhs_cofactors = []

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "NCCCCCC(=O)NCCCCCC(=O)O + O = NCCCCCC(=O)O + NCCCCCC(=O)O"

def test_construct_rxn_str_11():
    substrates = ["Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O"]
    lhs_cofactors = ["O"]
    products = ["Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O"]
    rhs_cofactors = ["O=P(O)(O)O"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O + O = Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O + O=P(O)(O)O"

def test_construct_rxn_str_12():
    substrates = ["C=C(N)C(=O)O"]
    lhs_cofactors = ["O"]
    products = ["NC(CO)C(=O)O"]
    rhs_cofactors = []

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "C=C(N)C(=O)O + O = NC(CO)C(=O)O"

def test_construct_rxn_str_13():
    substrates = ["CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCSC(=O)CCC(=O)O",
                  "NC(CCO)C(=O)O"]
    lhs_cofactors = []
    products = ["NC(CCOC(=O)CCC(=O)O)C(=O)O"]
    rhs_cofactors = ["CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O)C(O)C(=O)NCCC(=O)NCCS"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == f"{substrates[0]} + {substrates[1]} = {products[0]} + {rhs_cofactors[0]}"

def test_construct_rxn_str_14():
    substrates = ["CC(=O)C=O"]
    lhs_cofactors = ["O",
                     "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"]
    products = ["CC(=O)C(=O)O"]
    rhs_cofactors = ["NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"]

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "CC(=O)C=O + O + NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1 = CC(=O)C(=O)O + NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

def test_construct_rxn_str_15():
    substrates = ["Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC(CO)CO)C(O)C2O)c(=O)n1"]
    lhs_cofactors = ["O=P(O)(O)OP(=O)(O)O"]
    products = ["O=P(O)(O)OC(CO)CO",
                "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)n1"]
    rhs_cofactors = []

    rxn_eq = preprocessing.construct_rxn_str(substrates = substrates,
                                             lhs_cofactors = lhs_cofactors,
                                             products = products,
                                             rhs_cofactors = rhs_cofactors)

    assert rxn_eq == "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC(CO)CO)C(O)C2O)c(=O)n1 + O=P(O)(O)OP(=O)(O)O = O=P(O)(O)OC(CO)CO + Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)n1"

# --------- Tests for 'remove_processed_rxns_with_asterisked_cpd_smiles' in preprocessing package ---------
def test_remove_processed_rxns_with_asterisked_cpd_smiles():
    processed_rxns = [{'Reaction ID': 'RXN-123', 'Substrates': ['C*C'], 'LHS_cofactors': ['O2'],
             'Products': ['CC'], 'RHS_cofactors': ['H2O']},
            {'Reaction ID': 'RXN-124', 'Substrates': ['CC'], 'LHS_cofactors': ['O2'],
             'Products': ['C*C'], 'RHS_cofactors': ['H2O']},
            {'Reaction ID': 'RXN-124', 'Substrates': ['CC'], 'LHS_cofactors': ['O2'],
                       'Products': ['CC'], 'RHS_cofactors': ['H2O']}
                      ]
    assert preprocessing.count_processed_rxns_with_asterisked_cpd_smiles(processed_rxns) == 2

# --------- Tests for 'remove_processed_rxns_with_asterisked_cpd_smiles' in preprocessing package ---------
def remove_processed_rxns_with_asterisked_cpd_smiles():
    processed_rxns = [
        {'Reaction ID': 'RXN-123', 'Substrates': ['C*C'], 'LHS_cofactors': ['O2'],
         'Products': ['CC'], 'RHS_cofactors': ['H2O']},
        {'Reaction ID': 'RXN-124', 'Substrates': ['CC'], 'LHS_cofactors': ['O2'],
         'Products': ['CC'], 'RHS_cofactors': ['H2O']}
    ]

    cleaned_rxns = preprocessing.remove_processed_rxns_with_asterisked_cpd_smiles(processed_rxns)
    assert len(cleaned_rxns) == 1
    assert cleaned_rxns == [{'Reaction ID': 'RXN-124', 'Substrates': ['CC'], 'LHS_cofactors': ['O2'],
         'Products': ['CC'], 'RHS_cofactors': ['H2O']}]