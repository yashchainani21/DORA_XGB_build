import numpy as np

from ML_utils import ML_utils

def test_initialize_rxn_fp_for_ecfp_by_ascending_MW():
    fp_type = 'ecfp4'
    max_species = 4
    cofactor_positioning_method = 'by_ascending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_ecfp_by_descending_MW():
    fp_type = 'ecfp4'
    max_species = 4
    cofactor_positioning_method = 'by_descending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_ecfp_by_add_concat():
    fp_type = 'ecfp4'
    max_species = 4
    cofactor_positioning_method = 'add_concat'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2*2048

def test_initialize_rxn_fp_for_ecfp_add_subtract():
    fp_type = 'ecfp4'
    max_species = 4
    cofactor_positioning_method = 'add_subtract'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2048

def test_initialize_rxn_fp_for_ecfp_half_random():
    fp_type = 'ecfp4'
    max_species = 4
    cofactor_positioning_method = 'half_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_ecfp_full_random():
    fp_type = 'ecfp4'
    max_species = 4
    cofactor_positioning_method = 'full_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_atom_pair_by_ascending_MW():
    fp_type = 'atom_pair'
    max_species = 4
    cofactor_positioning_method = 'by_ascending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_atom_pair_by_descending_MW():
    fp_type = 'atom_pair'
    max_species = 4
    cofactor_positioning_method = 'by_descending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_atom_pair_by_add_concat():
    fp_type = 'atom_pair'
    max_species = 4
    cofactor_positioning_method = 'add_concat'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2*2048

def test_initialize_rxn_fp_for_atom_pair_add_subtract():
    fp_type = 'atom_pair'
    max_species = 4
    cofactor_positioning_method = 'add_subtract'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2048

def test_initialize_rxn_fp_for_atom_pair_half_random():
    fp_type = 'atom_pair'
    max_species = 4
    cofactor_positioning_method = 'half_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_atom_pair_full_random():
    fp_type = 'ecfp4'
    max_species = 4
    cofactor_positioning_method = 'full_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MACCS_by_ascending_MW():
    fp_type = 'MACCS'
    max_species = 4
    cofactor_positioning_method = 'by_ascending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*167

def test_initialize_rxn_fp_for_MACCS_by_descending_MW():
    fp_type = 'MACCS'
    max_species = 4
    cofactor_positioning_method = 'by_descending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*167

def test_initialize_rxn_fp_for_MACCS_by_add_concat():
    fp_type = 'MACCS'
    max_species = 4
    cofactor_positioning_method = 'add_concat'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2*167

def test_initialize_rxn_fp_for_MACCS_add_subtract():
    fp_type = 'MACCS'
    max_species = 4
    cofactor_positioning_method = 'add_subtract'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 167

def test_initialize_rxn_fp_for_MACCS_half_random():
    fp_type = 'MACCS'
    max_species = 4
    cofactor_positioning_method = 'half_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*167

def test_initialize_rxn_fp_for_MACCS_full_random():
    fp_type = 'MACCS'
    max_species = 4
    cofactor_positioning_method = 'full_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*167

def test_initialize_rxn_fp_for_mordred_by_ascending_MW():
    fp_type = 'mordred'
    max_species = 4
    cofactor_positioning_method = 'by_ascending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*1613

def test_initialize_rxn_fp_for_mordred_by_descending_MW():
    fp_type = 'mordred'
    max_species = 4
    cofactor_positioning_method = 'by_descending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*1613

def test_initialize_rxn_fp_for_mordred_by_add_concat():
    fp_type = 'mordred'
    max_species = 4
    cofactor_positioning_method = 'add_concat'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2*1613

def test_initialize_rxn_fp_for_mordred_add_subtract():
    fp_type = 'mordred'
    max_species = 4
    cofactor_positioning_method = 'add_subtract'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 1613

def test_initialize_rxn_fp_for_mordred_half_random():
    fp_type = 'mordred'
    max_species = 4
    cofactor_positioning_method = 'half_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*1613

def test_initialize_rxn_fp_for_mordred_full_random():
    fp_type = 'mordred'
    max_species = 4
    cofactor_positioning_method = 'full_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*1613

def test_initialize_rxn_fp_for_MAP4_by_ascending_MW():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'by_ascending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MAP4_by_descending_MW():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'by_descending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MAP4_by_add_concat():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'add_concat'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2*2048

def test_initialize_rxn_fp_for_MAP4_add_subtract():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'add_subtract'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2048

def test_initialize_rxn_fp_for_MAP4_half_random():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'half_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MAP4_full_random():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'full_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MHFP4_by_ascending_MW():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'by_ascending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MHFP4_by_descending_MW():
    fp_type = 'MAP4'
    max_species = 4
    cofactor_positioning_method = 'by_descending_MW'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MHFP4_by_add_concat():
    fp_type = 'MHFP4'
    max_species = 4
    cofactor_positioning_method = 'add_concat'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2*2048

def test_initialize_rxn_fp_for_MHFP4_add_subtract():
    fp_type = 'MHFP4'
    max_species = 4
    cofactor_positioning_method = 'add_subtract'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 2048

def test_initialize_rxn_fp_for_MHFP4_half_random():
    fp_type = 'MHFP4'
    max_species = 4
    cofactor_positioning_method = 'half_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048

def test_initialize_rxn_fp_for_MHFP4_full_random():
    fp_type = 'MHFP4'
    max_species = 4
    cofactor_positioning_method = 'full_random'

    initialized_fp, fp_length = ML_utils.initialize_rxn_fp(cofactor_positioning_method = cofactor_positioning_method,
                                                           fp_type = fp_type,
                                                           max_species = max_species)

    assert type(initialized_fp) == np.ndarray
    assert sum(initialized_fp) == 0
    assert fp_length == 8*2048