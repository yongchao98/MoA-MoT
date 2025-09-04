import sys
from io import StringIO

def check_nmr_answer():
    """
    This function checks the correctness of the identified compound based on 1H NMR data.
    It verifies the answer against several constraints derived from the spectrum.
    """
    # --- Data from the Question ---
    # 1H NMR: 7.0 (1H, d, J = 16.0 Hz), 5.5 (1H, dq), 2.1 (3H, s), 1.6 (3H, d)
    nmr_data = {
        'vinylic_H_a': {'integration': 1, 'multiplicity': 'd', 'J_Hz': 16.0},
        'vinylic_H_b': {'integration': 1, 'multiplicity': 'dq'},
        'acetate_CH3': {'integration': 3, 'multiplicity': 's'},
        'propenyl_CH3': {'integration': 3, 'multiplicity': 'd'},
    }

    # --- Candidate Compounds ---
    # Options: A) Cis-propenyl acetate, B) Trans-propenyl acetate, 
    #          C) Trans-butenyl acetate, D) Cis-butenyl acetate
    candidates = {
        'A': {'name': 'Cis-propenyl acetate', 'protons': 8, 'stereochemistry': 'cis', 'chain': 'propenyl'},
        'B': {'name': 'Trans-propenyl acetate', 'protons': 8, 'stereochemistry': 'trans', 'chain': 'propenyl'},
        'C': {'name': 'Trans-butenyl acetate', 'protons': 10, 'stereochemistry': 'trans', 'chain': 'butenyl'},
        'D': {'name': 'Cis-butenyl acetate', 'protons': 10, 'stereochemistry': 'cis', 'chain': 'butenyl'},
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_key = 'B'
    
    # --- Verification Process ---
    
    # Constraint 1: Total Proton Count
    # The sum of integrations must match the number of protons in the molecule.
    total_protons_observed = sum(signal['integration'] for signal in nmr_data.values())
    
    # Filter candidates based on proton count
    possible_candidates_keys = [key for key, props in candidates.items() if props['protons'] == total_protons_observed]
    
    if not possible_candidates_keys:
        return f"Incorrect. The observed proton count ({total_protons_observed}H) does not match any of the candidates."
        
    if 'C' in possible_candidates_keys or 'D' in possible_candidates_keys:
        return f"Incorrect. The analysis of proton count is flawed. The observed {total_protons_observed}H should eliminate butenyl acetates (10H)."

    # Constraint 2: Stereochemistry from J-coupling
    # The coupling constant (J) between vinylic protons indicates cis or trans geometry.
    # Trans: ~12-18 Hz; Cis: ~6-12 Hz
    j_coupling = nmr_data['vinylic_H_a']['J_Hz']
    observed_stereochemistry = None
    if 12 <= j_coupling <= 18:
        observed_stereochemistry = 'trans'
    elif 6 <= j_coupling < 12:
        observed_stereochemistry = 'cis'
    else:
        return f"Incorrect. The J-coupling constant ({j_coupling} Hz) is outside typical ranges or was not interpreted correctly."

    # Filter the remaining candidates based on stereochemistry
    final_candidate_keys = [key for key in possible_candidates_keys if candidates[key]['stereochemistry'] == observed_stereochemistry]

    # Constraint 3: Fragmentation Pattern
    # The splitting patterns must match the structure.
    # - 3H singlet -> Acetate CH3 (present in all)
    # - 3H doublet -> CH3 next to a single H (confirms CH3-CH= fragment, i.e., propenyl)
    # - 1H dq -> Vinylic H next to a CH3 and another vinylic H (confirms -CH=CH-CH3)
    # This confirms the chain is 'propenyl', which was already established by proton count.
    # This step serves as a consistency check.
    is_propenyl = all(candidates[key]['chain'] == 'propenyl' for key in final_candidate_keys)
    if not is_propenyl:
        return "Incorrect. The splitting pattern (3H doublet, 1H dq) points to a propenyl group, but a butenyl candidate was selected."

    # --- Final Verdict ---
    if len(final_candidate_keys) != 1:
        return f"Incorrect. The analysis is inconclusive. After applying all constraints, {len(final_candidate_keys)} candidates remain: {final_candidate_keys}."
        
    correct_answer_key = final_candidate_keys[0]

    if llm_answer_key == correct_answer_key:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{llm_answer_key}' ({candidates[llm_answer_key]['name']}) is incorrect.\n"
            f"The correct answer is '{correct_answer_key}' ({candidates[correct_answer_key]['name']}).\n"
            f"Reasoning:\n"
            f"1. Proton Count: The spectrum shows {total_protons_observed} protons, which eliminates the butenyl acetate options (C and D) that have 10 protons.\n"
            f"2. Stereochemistry: The vinylic coupling constant J = {j_coupling} Hz is characteristic of a 'trans' configuration. This eliminates the 'cis' option (A).\n"
            f"Therefore, the only structure that satisfies all constraints is Trans-propenyl acetate (B)."
        )
        return f"Incorrect. {reason}"

# Run the check and print the result
result = check_nmr_answer()
print(result)