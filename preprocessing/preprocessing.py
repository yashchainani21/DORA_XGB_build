import json
import pandas as pd
from typing import Tuple, Optional
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from typing import List, Dict

try:
    from featurizations import featurizations
except ImportError:
    import sys
    sys.path.append('../featurizations')
    import featurizations

# Silence non-critical RDKit warnings to minimize unnecessary outputs
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# All checks are manual
# # All unit tests are in '../tests/02_preprocessing_test/all_preprocessing_tests.py

# Checked and tested
def load_cofactors(cofactors_filepath: str) -> Tuple[pd.DataFrame, List[str]]:
    cofactors_df = pd.read_csv(cofactors_filepath, delimiter=',')
    cofactors_list = list(pd.read_csv(cofactors_filepath, delimiter=',')['SMILES'])
    return cofactors_df, cofactors_list

# Checked and tested
def load_dataset(db: str) -> dict:
    """
    Load enzymatic reaction datasets from specified databases into a Python dictionary.

    Parameters:
    - db (str): The type of enzymatic database to load. Supported database options include:
        - "MetaCyc": For loading the MetaCyc database.
        - "BRENDA": For loading the BRENDA database.
        - "KEGG": For loading the KEGG database.
        - "MetaCyc_V24": For loading version 24 of the MetaCyc database.

    Returns:
    - raw_data (dict): A dictionary containing the reactions within the specified enzymatic database.
      The structure of the dictionary depends on the JSON data structure within the source files for each database.

    Raises:
    - FileNotFoundError: If the specified database file does not exist at the expected location.
    - json.JSONDecodeError: If there is an issue with decoding the JSON file.

    Note:
    The JSON files for each database are expected to be located in a `../data/raw/` directory relative to the
    script executing this function.
    """
    if db == "metacyc":
        try:
            metacyc_filepath = "../data/raw/metacyc.json"
            with open(metacyc_filepath) as f:
                metacyc_raw_data = json.load(f)

        except FileNotFoundError:
            metacyc_filepath = "../../data/raw/metacyc.json"
            with open(metacyc_filepath) as f:
                metacyc_raw_data = json.load(f)

        return metacyc_raw_data

    elif db == "brenda":
        try:
            brenda_filepath = "../data/raw/brenda.json"
            with open(brenda_filepath) as f:
                brenda_raw_data = json.load(f)

        except FileNotFoundError: # need to go up an additional directory if running a test
            brenda_filepath = "../../data/raw/brenda.json"
            with open(brenda_filepath) as f:
                brenda_raw_data = json.load(f)

        return brenda_raw_data

    elif db == "kegg":
        try:
            kegg_filepath = "../data/raw/kegg.json"
            with open(kegg_filepath) as f:
                kegg_raw_data = json.load(f)

        except FileNotFoundError: # need to go up an additional directory if running a test
            kegg_filepath = "../../data/raw/kegg.json"
            with open(kegg_filepath) as f:
                kegg_raw_data = json.load(f)

        return kegg_raw_data


    elif db == "metacyc_V24":
        try:
            metacyc_v24_filepath = "../data/raw/metacyc24.json"
            with open(metacyc_v24_filepath) as f:
                metacyc_v24_raw_data = json.load(f)
        except:
            metacyc_v24_filepath = "../../data/raw/metacyc24.json"
            with open(metacyc_v24_filepath) as f:
                metacyc_v24_raw_data = json.load(f)

        return metacyc_v24_raw_data

    else:
        raise ValueError(f"Database '{db}' is not supported.")

# Checked and tested
def get_rule_info(rules_df: pd.DataFrame, query_rule: str) -> Tuple[int, int, List[str], List[str]]:
    """
    Extracts and returns detailed information about a specified reaction rule from a reaction rules DataFrame.

    Parameters
    - rules_df : pd.DataFrame
        - A pandas DataFrame containing reaction rules. This DataFrame should have columns named "Name",
        - "Reactants", and "Products", where "Reactants" and "Products" are semicolon-separated strings
        - representing the number and types of reactants and products in each rule.
    - query_rule : str
        - The name of the rule for which information is to be retrieved. It should match a value in the
        - "Name" column of `rules_df`.

    Returns
    - num_substrates : int
        - The number of substrates (indicated by "Any" in the "Reactants" column) in the specified rule.
    - num_products : int
        - The number of products (indicated by "Any" in the "Products" column) in the specified rule.
    - lhs_cofactors : list of str
        - A list of cofactors on the reactants' side (left-hand side) of the reaction for the specified rule.
        - For example, ['NAD_CoF'] for a rule with NAD_CoF as a reactant cofactor.
    - rhs_cofactors : list of str
        - A list of cofactors on the products' side (right-hand side) of the reaction for the specified rule.
        - For example, ['NADH_CoF'] for a rule with NADH_CoF as a product cofactor.

    Notes
    -----
    This function assumes that substrates and products specified as "Any" in the "Reactants" and "Products"
    columns are placeholders for actual substrates and products, while other entries are considered as cofactors.
    The function splits the "Reactants" and "Products" strings on semicolons to separate individual components.
    """
    selected_rule = rules_df[rules_df["Name"] == query_rule]
    reactants = list(selected_rule["Reactants"])[0].split(";")
    products = list(selected_rule["Products"])[0].split(";")

    # Count number of substrates on LHS then remove them, leaving only cofactors
    num_substrates = reactants.count("Any")
    for i in range(num_substrates):
        reactants.remove("Any")
    lhs_cofactors = reactants

    # Count number of products on RHS then remove them, leaving only cofactors
    num_products = products.count("Any")
    for i in range(num_products):
        products.remove("Any")
    rhs_cofactors = products

    return num_substrates, num_products, lhs_cofactors, rhs_cofactors

# Checked and tested
def lookup_cofactor_by_ID(ID: str, cofactors_df: pd.DataFrame) -> str:
    return list(cofactors_df[cofactors_df['#ID'] == ID]['SMILES'])[0]

# Checked and tested
def cofactor_ids_to_smiles(cofactor_ids_list: list, cofactors_df: pd.DataFrame) -> list:
    """
    Converts a list of cofactor IDs into their corresponding SMILES strings using the provided cofactors DataFrame.

    This function iterates through each cofactor ID in the given list and retrieves its associated SMILES
    string from the DataFrame. For special cofactors like NAD+ and NADH, the function adds both their commonly
    encountered forms (NAD+, NADH) and their official Pickaxe V2.0 forms (NADP+, NADPH) to the output list.

    Parameters:
        cofactor_ids_list (list): A list of strings representing the cofactor IDs.
        cofactors_df (pd.DataFrame): A pandas DataFrame containing at least two columns: one for the cofactor IDs
                                     and another for the corresponding SMILES strings. The DataFrame is expected
                                     to have columns named 'ID' and 'SMILES', where 'ID' matches the cofactor IDs
                                     and 'SMILES' contains the corresponding SMILES strings.

    Returns:
        list: A list of SMILES strings corresponding to the input cofactor IDs. For NAD+ and NADH cofactors,
              the list includes both their common and official Pickaxe V2.0 forms.

    Notes:
        The function relies on the assumption that the DataFrame provided contains accurate and up-to-date
        information regarding the cofactors and their corresponding SMILES strings. The user is responsible for
        ensuring that the 'ID' and 'SMILES' columns in the DataFrame are correctly formatted and populated.
    """

    cofactor_smiles_list = []
    for cofactor_id in cofactor_ids_list:

        # if a rule contains an NAD+ cofactor
        if cofactor_id == 'NAD_CoF':
            # add in both NAD+ (often encountered in mined reactions) and NADP+ (official Pickaxe V2.0 cofactor)
            nad_plus_smiles = lookup_cofactor_by_ID(ID = 'NAD', cofactors_df = cofactors_df)
            cofactor_smiles_list.append(nad_plus_smiles)

            nadp_plus_smiles = lookup_cofactor_by_ID(ID = 'NAD_CoF', cofactors_df = cofactors_df)
            cofactor_smiles_list.append(nadp_plus_smiles)

        # if a rule contains an NADH cofactor
        if cofactor_id == 'NADH_CoF':
            # add in both NADH (often encountered in mined reactions) and NADPH (official Pickaxe V2.0 cofactor)
            nadh_smiles = lookup_cofactor_by_ID(ID = 'NADH', cofactors_df = cofactors_df)
            cofactor_smiles_list.append(nadh_smiles)

            nadph_smiles = lookup_cofactor_by_ID(ID = 'NADH_CoF', cofactors_df = cofactors_df)
            cofactor_smiles_list.append(nadph_smiles)

        if cofactor_id != 'NAD_CoF' and cofactor_id != 'NADH_CoF':
            cofactor_smiles = lookup_cofactor_by_ID(ID = cofactor_id, cofactors_df = cofactors_df)
            cofactor_smiles_list.append(cofactor_smiles)

    return cofactor_smiles_list

# Checked and tested
def fix_cofactors(mined_cofactor_smiles: str, cofactors_list: list) -> Optional[str]:
    """
    Often, the SMILES strings of cofactors in mined versions of Brenda, KEGG, and MetaCyc have missing atoms
    Such missing atoms are typically indicated with an asterisk '*'
    If such SMILES strings are copied and pasted into Chemdraw, the asterisks will often appear as question marks
    Sometimes, the asterisks will also correspond to incorrect atoms (such as an 'A' in the incorrect NAD+ SMILES)
    This function accepts a cofactor SMILES string and if incorrect, fixes it by swapping the incorrect SMILES
    with the correct version from our pre-compiled list of cofactors without any stereochemical information
    This list can be found in '../data/processed/cofactors_processed.tsv'

    Even if a cofactor may not contain missing atoms as denoted by an asterisk, it may still be an uncommon cofactor
    In such cases, we still swap out its SMILES for that present in our cofactors list
    For instance, GDP, CDP, NDP are all common glucosyl donors
    But we only use UDP in our cofactors list
    As such, this function will swap out GDP, CDP, or NDP or UDP

    Parameters
    ----------
    mined_cofactor_smiles : str
        SMILES string of cofactor that was mined from either BRENDA, KEGG, or MetaCyc

    cofactors_list: list
        List of cofactors established within Pickaxe V2.0

    Returns
    -------
    correct_cofactor_smiles : str
    """
    # --------- NAD+ correction ---------
    nad_plus_incorrect = "*OC1C(O)C(COP(=O)(O)OP(=O)(O)OCC2OC([n+]3cccc(C(N)=O)c3)C(O)C2O)OC1n1cnc2c(N)ncnc21"
    nad_plus_correct = "NC(=O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)c1"

    if mined_cofactor_smiles == nad_plus_incorrect:
        return nad_plus_correct

    # --------- NADH correction ---------
    nadh_incorrect = "*OC1C(O)C(COP(=O)(O)OP(=O)(O)OCC2OC(N3C=CCC(C(N)=O)=C3)C(O)C2O)OC1n1cnc2c(N)ncnc21"
    nadh_correct = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(O)C3O)C(O)C2O)C=CC1"

    if mined_cofactor_smiles == nadh_incorrect:
        return nadh_correct

    # --------- glucosyl acceptors that sometimes appears in BRENDA (MetaCyc and KEGG should be fine) ---------
    GDP_smiles = "Nc1nc(=O)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)O)C(O)C3O)c2[nH]1"
    CDP_smiles = "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)n1"
    NDP_smiles = "CC1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O"
    dTDP_smiles = "Cc1cn(C2CC(O)C(COP(=O)(O)OP(=O)(O)O)O2)c(=O)[nH]c1=O"

    # the glucosyl acceptor officially in our cofactors list
    UDP_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"

    if mined_cofactor_smiles in [GDP_smiles, CDP_smiles, NDP_smiles, dTDP_smiles]:
        return UDP_smiles

    # --------- glucosyl donors that sometimes appears in BRENDA (MetaCyc and KEGG should be fine) ---------
    GDP_mannose_smiles = "Nc1nc(=O)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)OC4OC(CO)C(O)C(O)C4O)C(O)C3O)c2[nH]1"
    CDP_alpha_D_glucose_smiles = "Nc1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)n1"
    NDP_glucose_smiles = "CC1OC(COP(=O)(O)OP(=O)(O)OC2OC(CO)C(O)C(O)C2O)C(O)C1O"
    dTDP_alpha_D_glucose_smiles = "Cc1cn(C2CC(O)C(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)O2)c(=O)[nH]c1=O"

    # glucosyl donor that is officially in our cofactors list
    CPD_12575_smiles = "O=c1ccn(C2OC(COP(=O)(O)OP(=O)(O)OC3OC(CO)C(O)C(O)C3O)C(O)C2O)c(=O)[nH]1"

    if mined_cofactor_smiles in [GDP_mannose_smiles, CDP_alpha_D_glucose_smiles, NDP_glucose_smiles,
                                 dTDP_alpha_D_glucose_smiles]:
        return CPD_12575_smiles

    # --------- phosphate acceptors that sometimes appear in either BRENDA, MetaCyc, or KEGG ---------
    phosphate_acceptor_incorrect = "*C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O"

    # phosphate acceptor that is officially in our cofactors list
    ADP_smiles = "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O"

    if mined_cofactor_smiles == phosphate_acceptor_incorrect:
        return ADP_smiles

    # --------- phosphate donors that sometimes appear in either BRENDA, MetaCyc, or KEGG ---------
    phosphate_donor_incorrect_1 = "*C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O" # missing 1 atom
    phosphate_donor_incorrect_2 = "*C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1*" # missing 2 atoms

    # phosphate donor that is officially in our cofactors list
    ATP_smiles = "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O"

    if mined_cofactor_smiles in [phosphate_donor_incorrect_1, phosphate_donor_incorrect_2]:
        return ATP_smiles

    if mined_cofactor_smiles in cofactors_list:
        return mined_cofactor_smiles

    else:
        return None

# Checked and tested
def process_reactants_dict(reactants_dict: Dict[str,str], cofactors_list: list) -> Tuple[List[str], List[str]]:

    lhs_cofactors_list = []
    substrates_list = []

    # iterate over all reactants on reaction LHS
    for key in reactants_dict.keys():
        reactant_smiles = reactants_dict[key]
        reactant_object = featurizations.compound(reactant_smiles)
        reactant_smiles_no_stereo = reactant_object.remove_stereo()

        # after removing stereochemical info, check if reactant is a cofactor within our cofactors list
        is_reactant_cofactor = fix_cofactors(mined_cofactor_smiles = reactant_smiles_no_stereo,
                                             cofactors_list = cofactors_list)

        if is_reactant_cofactor:
            lhs_cofactors_list.append(is_reactant_cofactor)

        # if reactant is not in our cofactors list, it counts as a substrate
        else:
            substrates_list.append(reactant_smiles_no_stereo)

    return lhs_cofactors_list,substrates_list

# Checked and tested
def process_products_dict(products_dict: Dict[str,str], cofactors_list: List[str]) -> Tuple[List[str], List[str]]:

    rhs_cofactors_list = []
    products_list = []

    # iterate over all products on reaction RHS
    for key in products_dict.keys():
        product_smiles = products_dict[key]
        product_object = featurizations.compound(product_smiles)
        product_smiles_no_stereo = product_object.remove_stereo()

        # after removing stereochemical info, check if product is a cofactor within our cofactors list
        is_product_cofactor = fix_cofactors(mined_cofactor_smiles = product_smiles_no_stereo,
                                            cofactors_list = cofactors_list)

        if is_product_cofactor:
            rhs_cofactors_list.append(is_product_cofactor)

        # if product is not in our cofactors list, it counts as a product of the reaction itself
        else:
            products_list.append(product_smiles_no_stereo)

    return rhs_cofactors_list, products_list

# Checked and tested
def reorder_smiles_by_mol_weight(smiles_list: List[str], descending: bool) -> List[str]:
    """
    Reorders a list of SMILES (Simplified Molecular Input Line Entry System) strings
    based on their molecular weights.

    Parameters:
    - smiles_list (list of str): A list of SMILES strings representing molecules.
    - descending (bool, optional): If True, the returned list is sorted in descending
      order of molecular weight. If False (default), the list is sorted in ascending order.

    Returns:
    - list of str: A new list of SMILES strings sorted by molecular weight in the specified order.

    The function calculates the molecular weight for each SMILES string using RDKit and sorts
    the list accordingly. The RDKit library must be installed and importable in the environment
    where this function is used.
    """
    # First, ensure all stereochemistry information has been removed from the SMILES strings
    try:
        assert all(['@' not in smiles for smiles in smiles_list])
    except:
        raise ValueError("Stereochemistry encountered. Ensure SMILES are canonicalized and have no stereochemistry.")

    # Calculate molecular weights and associate them with SMILES strings
    smiles_with_weights = [(smiles, MolWt(Chem.MolFromSmiles(smiles))) for smiles in smiles_list]

    # Sort the list by molecular weight
    # The key parameter is set to a lambda function that returns the molecular weight
    sorted_smiles = sorted(smiles_with_weights, key=lambda x: x[1], reverse=descending)

    # Return only the sorted SMILES strings, excluding weights
    return [smiles for smiles, weight in sorted_smiles]

# Checked and tested
def construct_rxn_str(substrates: List[str], lhs_cofactors: List[str],
                      products: List[str], rhs_cofactors: List[str]) -> str:
    """
    Construct a reaction string of the form "substrate_smiles + cofactor_smiles = product_smiles + cofactor_smiles"
    Standardizing this arrangement will be helpful for featurizing reactions for machine learning
    This order can definitely be switch but it is important that it be standardized across all reactions
    This reaction string representation is also helpful for calculating reaction thermodynamics

    :param substrates_rxn (list): List of substrates on the LHS of the reaction
    :param lhs_cofactors_rxn (list): List of cofactors on the LHS of the reaction
    :param products_rxn (list): List of products on the RHS of the reaction
    :param rhs_cofactors_rxn (list): List of cofactors on the RHS of the reaction
    :return rxn_str (rxn): string of the form "substrate_smiles + cofactor_smiles = product_smiles + cofactor_smiles"
    """
    rxn_str = ""

    # Reactants first on LHS
    for substrate_smi in substrates:
        rxn_str += substrate_smi
        rxn_str += " + "

    # then cofactors on LHS
    for lhs_cofactor_smi in lhs_cofactors:
        rxn_str += lhs_cofactor_smi
        rxn_str += " + "

    # Remove extra ' + ' from LHS
    rxn_str = rxn_str.rstrip(" + ")

    # Insert ' = ' before switching over to RHS
    rxn_str += " = "

    # Products first on RHS
    for product_smi in products:
        rxn_str += product_smi
        rxn_str += " + "

    # then cofactors on RHS
    for rhs_cofactor_smi in rhs_cofactors:
        rxn_str += rhs_cofactor_smi
        rxn_str += " + "

    # Remove extra ' + ' from RHS
    rxn_str = rxn_str.rstrip(" + ")

    return rxn_str

# Checked and tested
def count_processed_rxns_with_asterisked_cpd_smiles(processed_rxns: List[Dict[str, str]]) -> int:
    """
    Counts the number of unique chemical reactions within a given list where at least one of the compounds
    (substrates, left-hand side cofactors, products, right-hand side cofactors) in the reaction contains
    an asterisk (*) in its SMILES representation. This asterisk often indicates a site of reactivity
    or a point of attachment in chemical structures, which may be undesirable for certain analyses.

    Parameters:
    - processed_rxns (List[Dict[str, str]]): A list of dictionaries, each representing a chemical reaction.
      These dictionaries should contain keys for 'Substrates', 'LHS_cofactors', 'Products',
      and 'RHS_cofactors', each associated with a list of SMILES strings.

    Returns:
    - int: The count of unique reaction IDs where at least one involved compound's SMILES string contains an asterisk.
    """
    rxns_ids_w_asterisk_cpds = []
    for rxn in processed_rxns:
        # Check each type of compound in the reaction
        for key in ['Substrates', 'LHS_cofactors', 'Products', 'RHS_cofactors']:
            cpds_list = rxn[key]
            # If any compound contains an asterisk, record the reaction ID
            if any('*' in cpd_smiles for cpd_smiles in cpds_list):
                print(rxn)
                rxns_ids_w_asterisk_cpds.append(rxn['Reaction ID'])
                break  # Stop checking this reaction if any compound has an asterisk
    return len(set(rxns_ids_w_asterisk_cpds))

# Checked and testing now
def remove_processed_rxns_with_asterisked_cpd_smiles(processed_rxns: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """
    Removes chemical reactions from a given list where any of the compounds
    (substrates, left-hand side cofactors, products, right-hand side cofactors)
    in the reaction contains an asterisk (*) in its SMILES representation.
    This function is useful for filtering out reactions that involve reactive
    sites or unspecified attachments, as indicated by asterisks in their
    chemical structures.

    Parameters:
    - processed_rxns (List[Dict[str, str]]): A list of dictionaries, each representing a chemical reaction.
      These dictionaries should contain keys for 'Substrates', 'LHS_cofactors', 'Products',
      and 'RHS_cofactors', each associated with a list of SMILES strings.

    Returns:
    - List[Dict[str, str]]: A list of dictionaries, each representing a chemical reaction,
      excluding those reactions where any involved compound's SMILES string contains an asterisk.
    """
    # Create a new list for reactions to keep
    rxns_to_keep = []
    for rxn in processed_rxns:
        # Flag to indicate whether a reaction should be removed
        remove_rxn = False
        for key in ['Substrates', 'LHS_cofactors', 'Products', 'RHS_cofactors']:
            cpds_list = rxn[key]
            if any('*' in cpd_smiles for cpd_smiles in cpds_list):
                remove_rxn = True
                break  # No need to check further if any compound has an asterisk
        if not remove_rxn:
            rxns_to_keep.append(rxn)  # Add to keep list if no asterisk was found
    return rxns_to_keep