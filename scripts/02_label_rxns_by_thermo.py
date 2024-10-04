import pickle
import copy

# critical driving force in KJ/ mol - rxns with mdf below this value are negative (label 0)
crit_DF = 10

# pick from either 'all_BKM_rxns', 'all_AdH_rxns', 'all_Monox_rxns'
dataset = "all_MetaCyc_V27_rxns"
input_filepath = None
output_filepath = None

if dataset == 'all_AdH_rxns':
    input_filepath = '../data/processed/all_AdH_rxns_with_thermo.pkl'
    output_filepath = f'../data/processed/all_AdH_rxns_with_thermo_labelled_crit_DF_{crit_DF}.pkl'

if dataset == 'all_Monox_rxns':
    input_filepath = '../data/processed/all_Monox_rxns_with_thermo.pkl'
    output_filepath = f'../data/processed/all_Monox_rxns_with_thermo_labelled_crit_DF_{crit_DF}.pkl'

if dataset == 'all_BKM_rxns':
    input_filepath = '../data/processed/all_processed_rxns_with_thermo.pkl'
    output_filepath = f'../data/processed/all_processed_rxns_with_thermo_labelled_crit_DF_{crit_DF}.pkl'

if dataset == 'all_MetaCyc_V24_rxns':
    input_filepath = '../data/processed/all_MetaCyc_V24_rxns_with_thermo.pkl'
    output_filepath = f'../data/processed/all_MetaCyc_V24_rxns_with_thermo_labelled_crit_DF_{crit_DF}.pkl'

if dataset == 'all_MetaCyc_V26_rxns':
    input_filepath = '../data/processed/all_MetaCyc_V26_rxns_with_thermo.pkl'
    output_filepath = f'../data/processed/all_MetaCyc_V26_rxns_with_thermo_labelled_crit_DF_{crit_DF}.pkl'

if dataset == 'all_MetaCyc_V27_rxns':
    input_filepath = '../data/processed/all_MetaCyc_V27_rxns_with_thermo.pkl'
    output_filepath = f'../data/processed/all_MetaCYc_V27_rxns_with_thermo_labelled_crit_DF_{crit_DF}.pkl'

with open(input_filepath,'rb') as in_file:
    input_file = pickle.load(in_file)

rxns_wo_thermo_count = 0 # count reactions for which we did not even calculate thermodynamics
rxns_w_high_uncertainty_thermo_count = 0 # count reactions for which thermo values exist but with high error
feasible_rxns_count = 0 # count thermodynamically feasible reactions per the specified MDF threshold
infeasible_rxns_count = 0 # count thermodynamically infeasible reactions per the specified MDF threshold

all_rxns = copy.deepcopy(input_file)

# for each processed reaction
for rxn in all_rxns:

    # if this reaction was not even sent to eQuilibrator, we cannot label it for training
    # for such reactions that are missing thermodynamic values, one of the following fields will be None
    if any( [rxn['Reaction eq in accession IDs'] is None,
             rxn['Physiological dG'] is None,
             rxn['Physiological dG error'] is None,
             rxn['Standard dG'] is None,
             rxn['Standard dG error'] is None] ):

        feasibility_label = None
        remark = 'No thermodynamic values calculated'
        rxn.update( {'feasibility_label': feasibility_label, 'remark': remark} )
        rxns_wo_thermo_count += 1

    # but if this reaction was sent to eQuilibrator
    else:
        # and it was returned with larger errors, we cannot label it either
        if rxn['Physiological dG error'] > 9999 and rxn['Standard dG error'] > 9999:
            feasibility_label = None
            remark = 'dG uncertainty too high'
            rxn.update({'feasibility_label': feasibility_label, 'remark': remark})
            rxns_w_high_uncertainty_thermo_count += 1

        # but if the errors in calculating thermodynamics are reasonable
        else:
            # then assign feasibility labels by the reaction driving forces
            min_DF = rxn['MDF']

            # if mdf more than or equal to critical value, reaction is feasible (can be made downhill at some conc)
            if min_DF >= crit_DF:
                feasibility_label = 1
                remark = 'pos rxn by thermo'
                rxn.update({'feasibility_label': feasibility_label, 'remark': remark})
                feasible_rxns_count += 1

            # if false, reaction is infeasible (i.e. will always be thermodynamically uphill)
            else:
                feasibility_label = 0
                remark = 'neg rxn by thermo'
                rxn.update({'feasibility_label': feasibility_label, 'remark': remark})
                infeasible_rxns_count += 1

print(f'\nReactions that could not be labelled due to the absence of thermodynamic values: {rxns_wo_thermo_count}')
print(f'\nReactions that could not be labelled because of large thermodynamic errors: {rxns_w_high_uncertainty_thermo_count}')
print(f'\nReactions that were found to be thermodynamically feasible: {feasible_rxns_count}')
print(f'\nReactions that were found to be thermodynamically infeasible: {infeasible_rxns_count}')

print('------------')
print(f'\nTotal number of reactions labelled: {feasible_rxns_count + infeasible_rxns_count}')

with open(output_filepath,'wb') as out_file:
    pickle.dump(all_rxns, out_file)