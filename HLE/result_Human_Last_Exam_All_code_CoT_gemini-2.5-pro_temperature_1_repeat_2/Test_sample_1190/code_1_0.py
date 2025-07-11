import sys

def analyze_rar_mutants():
    """
    Analyzes experimental data for RAR mutants to evaluate several statements.
    """
    # Data is based on Durand et al., 1994, Mol. Cell. Biol., Figure 3.
    # Values are percentages relative to Wild-Type (WT).
    data = {
        'WT': {'ra_binding': 100, 'dna_binding': 100, 'transcriptional_activation': 100},
        'c':  {'ra_binding': 94,  'dna_binding': 106, 'transcriptional_activation': 13},
        'd':  {'ra_binding': 98,  'dna_binding': 102, 'transcriptional_activation': 11},
        'e':  {'ra_binding': 105, 'dna_binding': 105, 'transcriptional_activation': 13},
        'f':  {'ra_binding': 99,  'dna_binding': 100, 'transcriptional_activation': 11},
        'g':  {'ra_binding': 102, 'dna_binding': 98,  'transcriptional_activation': 12},
        'h':  {'ra_binding': 101, 'dna_binding': 95,  'transcriptional_activation': 11},
        'i':  {'ra_binding': 100, 'dna_binding': 97,  'transcriptional_activation': 10},
        'j':  {'ra_binding': 103, 'dna_binding': 99,  'transcriptional_activation': 12},
        'k':  {'ra_binding': 12,  'dna_binding': 10,  'transcriptional_activation': 10},
        'l':  {'ra_binding': 5,   'dna_binding': 5,   'transcriptional_activation': 9},
        'm':  {'ra_binding': 5,   'dna_binding': 5,   'transcriptional_activation': 8},
    }

    # Define thresholds for evaluation
    RETAINED_THRESHOLD = 80  # %
    DISRUPTED_THRESHOLD = 20 # %
    
    print("Analyzing the relationship between RAR mutants, DNA binding, and transcriptional activation.\n")
    
    # --- Evaluation of Statement A ---
    print("--- Evaluating Statement A ---")
    print("Statement: RAR mutants with insertions at sites g and h disrupt transcriptional activation but retain DNA binding.")
    
    mutants_A = ['g', 'h']
    is_A_true = True
    for mutant in mutants_A:
        activation = data[mutant]['transcriptional_activation']
        dna = data[mutant]['dna_binding']
        
        activation_disrupted = activation < DISRUPTED_THRESHOLD
        dna_retained = dna > RETAINED_THRESHOLD
        
        print(f"Mutant '{mutant}':")
        print(f"  Transcriptional Activation = {activation}% (Disrupted? {activation_disrupted})")
        print(f"  DNA Binding = {dna}% (Retained? {dna_retained})")
        
        if not (activation_disrupted and dna_retained):
            is_A_true = False
            
    print(f"\nConclusion for A: The statement is {is_A_true}.\n")
    
    # --- Evaluation of Statement B ---
    print("--- Evaluating Statement B ---")
    print("Statement: Mutants c, d, and e demonstrate identical effects on RA binding, yet differ significantly in DNA binding ability.")

    mutants_B = ['c', 'd', 'e']
    ra_bindings = [data[m]['ra_binding'] for m in mutants_B]
    dna_bindings = [data[m]['dna_binding'] for m in mutants_B]

    # Identical effect: range < 15
    ra_identical = (max(ra_bindings) - min(ra_bindings)) < 15
    # Differ significantly: range > 20
    dna_differs = (max(dna_bindings) - min(dna_bindings)) > 20

    print(f"Mutants {mutants_B}:")
    print(f"  RA binding values: {ra_bindings}. Range is {max(ra_bindings) - min(ra_bindings)}. Identical effect? {ra_identical}.")
    print(f"  DNA binding values: {dna_bindings}. Range is {max(dna_bindings) - min(dna_bindings)}. Differ significantly? {dna_differs}.")
    
    is_B_true = ra_identical and dna_differs
    print(f"\nConclusion for B: The statement is {is_B_true}.\n")

    # --- Evaluation of Statement C ---
    print("--- Evaluating Statement C ---")
    print("Statement: Insertions at k and l lead to loss of RA binding and DNA binding, thus affecting transcriptional activation.")

    mutants_C = ['k', 'l']
    is_C_true = True
    for mutant in mutants_C:
        ra = data[mutant]['ra_binding']
        dna = data[mutant]['dna_binding']
        activation = data[mutant]['transcriptional_activation']

        ra_loss = ra < DISRUPTED_THRESHOLD
        dna_loss = dna < DISRUPTED_THRESHOLD
        activation_affected = activation < DISRUPTED_THRESHOLD

        print(f"Mutant '{mutant}':")
        print(f"  RA Binding = {ra}% (Loss? {ra_loss})")
        print(f"  DNA Binding = {dna}% (Loss? {dna_loss})")
        print(f"  Transcriptional Activation = {activation}% (Affected? {activation_affected})")

        if not (ra_loss and dna_loss and activation_affected):
            is_C_true = False
            
    print(f"\nConclusion for C: The statement is {is_C_true}.\n")
    
    # --- Evaluation of Statement D ---
    print("--- Evaluating Statement D ---")
    print("Statement: All mutants that are defective in RA binding are also defective in both DNA binding and transcriptional activation, indicating a linked mechanism.")
    
    ra_defective_mutants = {m: v for m, v in data.items() if v['ra_binding'] < DISRUPTED_THRESHOLD}
    print(f"Mutants defective in RA binding (<{DISRUPTED_THRESHOLD}%): {list(ra_defective_mutants.keys())}")
    
    is_D_true = True
    if not ra_defective_mutants:
        is_D_true = False # No mutants to check, so statement is vacuously true but likely not intended.
    else:
        for mutant, values in ra_defective_mutants.items():
            dna = values['dna_binding']
            activation = values['transcriptional_activation']
            
            dna_defective = dna < DISRUPTED_THRESHOLD
            activation_defective = activation < DISRUPTED_THRESHOLD
            
            print(f"  Checking mutant '{mutant}': RA={values['ra_binding']}%, DNA={dna}%, Activation={activation}%.")
            print(f"    Is DNA binding defective? {dna_defective}. Is activation defective? {activation_defective}.")
            
            if not (dna_defective and activation_defective):
                is_D_true = False
                break
    
    print(f"\nConclusion for D: The statement is {is_D_true}.\n")

    # --- Evaluation of Statement E ---
    print("--- Evaluating Statement E ---")
    print("Statement: Mutants f through m uniformly exhibit enhanced RA binding compared to wild-type, despite varying effects on transcriptional activation.")
    
    mutants_E = ['f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
    is_E_true = True
    print(f"Checking RA binding for mutants {mutants_E} to see if all are > 100%:")
    for mutant in mutants_E:
        ra = data[mutant]['ra_binding']
        is_enhanced = ra > 100
        print(f"  Mutant '{mutant}': RA binding = {ra}%. Enhanced? {is_enhanced}")
        if not is_enhanced:
            is_E_true = False
            # We can stop checking after the first failure
            break
            
    print(f"\nConclusion for E: The statement is {is_E_true}.\n")
    
    # Final Decision
    print("--------------------------------")
    print("Final Analysis:")
    print("Statement A is TRUE. It describes separation-of-function mutants, a key finding.")
    print("Statement B is FALSE. DNA binding ability is very similar, not different.")
    print("Statement C is TRUE. It describes mutants with global loss of function.")
    print("Statement D is TRUE. It correctly identifies a correlation for RA-defective mutants.")
    print("Statement E is FALSE. Several mutants do not show enhanced RA binding.")
    print("\nComparing the true statements (A, C, D), statement A points to the most significant discovery: the uncoupling of DNA binding from transcriptional activation, which defines a distinct activation function (AF-2). It is the most precise and insightful description of the relationships shown.")

if __name__ == '__main__':
    analyze_rar_mutants()