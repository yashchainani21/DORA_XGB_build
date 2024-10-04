"""
This script creates a table of cofactors and their corresponding SMILES strings without any stereochemistry.
An initial set of cofactors is obtained from JN1224MIN and processed here. NADH and NAD(+) are also added.
These cofactors allow DORA_XGB models to first parse an input reaction string into substrates, products, and cofactors.
Subsequently, molecular fingerprints of participating species can be organized into reaction fingerprints.
This organization is based on the configuration chosen with which to arrange species in creating reaction fingerprints.
"""
import pandas as pd
from featurizations import featurizations

raw_cofactors_list = '../data/raw/all_cofactors.tsv'
cofactors_df = pd.read_csv(raw_cofactors_list, delimiter='\t')
cofactors_smiles_list = list(cofactors_df['SMILES'])
cofactors_smiles_list_without_stereo = []

for cofactor_smiles in cofactors_smiles_list:
    cpd_object = featurizations.compound(cofactor_smiles)
    smi = cpd_object.remove_stereo() # will canonicalize SMILES string first then remove stereo
    cofactors_smiles_list_without_stereo.append(smi)

# replace SMILES list
cofactors_df['SMILES'] = cofactors_smiles_list_without_stereo

# also add NAD+ and NADH in our cofactors list underneath (only NADP+ and NADPH present currently)
nadplus_correct = ("NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1")

nadh_correct = ("NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1")

NAD_plus_entry = {'#ID': 'NAD', 'Name': 'NAD', 'SMILES': nadplus_correct}
NADH_entry = {'#ID': 'NADH', 'Name': 'NADH', 'SMILES': nadh_correct}

cofactors_df.loc[cofactors_df.shape[0]] = NAD_plus_entry # insert in NAD+
cofactors_df.loc[cofactors_df.shape[0] + 1] = NADH_entry # insert in NADH

cofactors_df.to_csv('../data/processed/expanded_cofactors_no_stereochem.tsv')