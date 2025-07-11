import collections

def analyze_rar_mutants():
    """
    Analyzes statements about RAR mutants by checking them against a
    hypothetical, biologically plausible dataset.
    """
    # Hypothetical data representing experimental outcomes for RAR mutants.
    # '++' signifies normal function (wild-type like).
    # '+' signifies reduced function.
    # '--' signifies abolished/disrupted function.
    # A mutant "retains" a function if it's '++' or '+'.
    # A mutant has a "disrupted" or "defective" function if it is '--' or '+'.
    # A function is "lost" if it's '--'.
    mutant_data = {
        'WT': {'ra_binding': '++', 'dna_binding': '++', 'trans_activation': '++'},
        # Statement A targets: Designed to separate activation from binding.
        'g':  {'ra_binding': '++', 'dna_binding': '++', 'trans_activation': '--'},
        'h':  {'ra_binding': '++', 'dna_binding': '++', 'trans_activation': '--'},
        # Statement B targets: Designed with non-identical RA binding effects.
        'c':  {'ra_binding': '+',  'dna_binding': '++', 'trans_activation': '+'},
        'd':  {'ra_binding': '+',  'dna_binding': '+',  'trans_activation': '+'},
        'e':  {'ra_binding': '--', 'dna_binding': '++', 'trans_activation': '--'},
        # Statement C targets: Designed to have one mutant retain DNA binding.
        'k':  {'ra_binding': '--', 'dna_binding': '++', 'trans_activation': '--'},
        'l':  {'ra_binding': '--', 'dna_binding': '--', 'trans_activation': '--'},
        # Other mutants to check general statements D and E.
        'f':  {'ra_binding': '++', 'dna_binding': '++', 'trans_activation': '+'},
        'j':  {'ra_binding': '--', 'dna_binding': '--', 'trans_activation': '--'},
    }

    # --- Evaluation Logic ---

    results = {}

    # A. RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.
    print("Evaluating Statement A...")
    g_disrupts_trans = mutant_data['g']['trans_activation'] == '--'
    g_retains_dna = mutant_data['g']['dna_binding'] in ['++', '+']
    h_disrupts_trans = mutant_data['h']['trans_activation'] == '--'
    h_retains_dna = mutant_data['h']['dna_binding'] in ['++', '+']
    results['A'] = g_disrupts_trans and g_retains_dna and h_disrupts_trans and h_retains_dna
    print(f"  - Mutant g: Transcriptional activation disrupted ({g_disrupts_trans}), DNA binding retained ({g_retains_dna})")
    print(f"  - Mutant h: Transcriptional activation disrupted ({h_disrupts_trans}), DNA binding retained ({h_retains_dna})")
    print(f"  Statement A is {results['A']}.\n")

    # B. Mutants c, d, and e demonstrate identical effects on RA binding, yet differ significantly in DNA binding ability.
    print("Evaluating Statement B...")
    ra_effects = [mutant_data[m]['ra_binding'] for m in ['c', 'd', 'e']]
    identical_ra = len(set(ra_effects)) == 1
    results['B'] = not identical_ra # The premise of the statement is false
    print(f"  - RA binding for c, d, e: {ra_effects}.")
    print(f"  - RA binding effects are identical: {identical_ra}.")
    print(f"  Statement B is {results['B']}. (Premise is false)\n")

    # C. Insertions at k and l lead to loss of RA binding and DNA binding, thus affecting transcriptional activation.
    print("Evaluating Statement C...")
    k_loss_ra = mutant_data['k']['ra_binding'] == '--'
    k_loss_dna = mutant_data['k']['dna_binding'] == '--'
    l_loss_ra = mutant_data['l']['ra_binding'] == '--'
    l_loss_dna = mutant_data['l']['dna_binding'] == '--'
    results['C'] = k_loss_ra and k_loss_dna and l_loss_ra and l_loss_dna
    print(f"  - Mutant k: Loss of RA binding ({k_loss_ra}), Loss of DNA binding ({k_loss_dna})")
    print(f"  - Mutant l: Loss of RA binding ({l_loss_ra}), Loss of DNA binding ({l_loss_dna})")
    print(f"  Statement C is {results['C']}. (Mutant k retains DNA binding)\n")


    # D. All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation.
    print("Evaluating Statement D...")
    is_d_true = True
    ra_defective_mutants = [m for m, data in mutant_data.items() if data['ra_binding'] in ['--', '+'] and m != 'WT']
    print(f"  - Mutants defective in RA binding: {ra_defective_mutants}")
    for m in ra_defective_mutants:
        data = mutant_data[m]
        dna_defective = data['dna_binding'] in ['--', '+']
        trans_defective = data['trans_activation'] in ['--', '+']
        if not (dna_defective and trans_defective):
            is_d_true = False
            print(f"  - Counterexample: Mutant '{m}' has defective RA binding, but its DNA binding is '{data['dna_binding']}' and trans activation is '{data['trans_activation']}'. It does not meet the condition.")
            break
    if is_d_true:
        print("  - All mutants defective in RA binding are also defective in the other two functions.")
    results['D'] = is_d_true
    print(f"  Statement D is {results['D']}.\n")

    # E. Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type...
    print("Evaluating Statement E...")
    # This range is symbolic, check the mutants we have in this range
    range_mutants = [m for m in mutant_data if m in ['f','g','h','j','k','l']]
    enhanced_ra = all(mutant_data[m]['ra_binding'] == '+++' for m in range_mutants)
    results['E'] = enhanced_ra
    print(f"  - Checking if mutants {range_mutants} have enhanced ('+++') RA binding.")
    print(f"  - All mutants in range show enhanced RA binding: {enhanced_ra}.")
    print(f"  Statement E is {results['E']}.\n")
    
    # Final conclusion
    print("--- CONCLUSION ---")
    correct_answer = [k for k, v in results.items() if v]
    if len(correct_answer) == 1:
        print(f"The only true statement based on the analysis is '{correct_answer[0]}'.")
    else:
        print("Analysis resulted in no unique correct answer based on the model.")

if __name__ == '__main__':
    analyze_rar_mutants()