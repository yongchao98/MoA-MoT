def check_spectroscopy_answer():
    """
    This function programmatically checks the correctness of the answer to a spectroscopy problem.
    It analyzes the given spectral data and evaluates each candidate structure against the evidence.
    """
    # --- 1. Deconstruct the Spectral Data into Constraints ---
    
    # IR Data:
    # - 3420 & 3325 cm-1: Two N-H stretches indicate a primary amine (-NH2).
    # - 1720 cm-1: Strong C=O stretch, characteristic of a conjugated ester. An amide C=O is typically lower (~1680 cm-1).
    ir_constraints = {
        'amine_type': 'primary',
        'carbonyl_type': 'ester'
    }

    # 1H NMR Data:
    # - 1.20 ppm (t, 3H) + 4.5 ppm (q, 2H): Classic signature of an ethyl group attached to an oxygen (-O-CH2-CH3).
    #   The high shift of the quartet (4.5 ppm) is key. A propionate (-CO-CH2-CH3) quartet would be ~2.5 ppm.
    # - 4.0 ppm (bs, 2H): Broad singlet for 2 protons, confirms the primary amine (-NH2).
    # - 7.0 ppm (d, 2H) + 8.0 ppm (d, 2H): Two doublets in the aromatic region indicate a 1,4- (para) disubstituted benzene ring.
    nmr_constraints = {
        'alkyl_group': 'ethyl_ester',  # Specifically -O-CH2-CH3
        'aromatic_pattern': 'para'
    }

    # --- 2. Define Properties of Each Candidate ---
    candidates = {
        'A': {
            'name': 'ethyl 4-aminobenzoate',
            'amine_type': 'primary',
            'carbonyl_type': 'ester',
            'alkyl_group': 'ethyl_ester',
            'aromatic_pattern': 'para'
        },
        'B': {
            'name': '3-ethoxybenzamide',
            'amine_type': 'primary_amide',
            'carbonyl_type': 'amide',
            'alkyl_group': 'ethoxy', # Different from ethyl ester
            'aromatic_pattern': 'meta'
        },
        'C': {
            'name': 'N-(4-ethoxyphenyl)formamide',
            'amine_type': 'secondary_amide',
            'carbonyl_type': 'amide',
            'alkyl_group': 'ethoxy',
            'aromatic_pattern': 'para'
        },
        'D': {
            'name': '4-aminophenyl propionate',
            'amine_type': 'primary',
            'carbonyl_type': 'ester',
            'alkyl_group': 'propionate_ester', # -O-CO-CH2-CH3
            'aromatic_pattern': 'para'
        }
    }

    # --- 3. Systematic Elimination ---
    correct_candidate = None
    reasons_for_failure = {}

    for key, props in candidates.items():
        failures = []
        # Check IR constraints
        if props['amine_type'] not in ['primary', 'primary_amide']:
            failures.append(f"it is a {props['amine_type']} but the IR data (two N-H bands) indicates a primary amine.")
        if props['carbonyl_type'] != ir_constraints['carbonyl_type']:
            failures.append(f"it is an {props['carbonyl_type']} but the IR C=O stretch at 1720 cm-1 indicates an ester.")
        
        # Check NMR constraints
        if props['aromatic_pattern'] != nmr_constraints['aromatic_pattern']:
            failures.append(f"it is {props['aromatic_pattern']}-substituted but the NMR shows a para pattern (two doublets).")
        if props['alkyl_group'] != nmr_constraints['alkyl_group']:
            failures.append(f"its alkyl group is '{props['alkyl_group']}', which is inconsistent with the NMR signal of an ethyl ester (quartet at 4.5 ppm).")

        if not failures:
            if correct_candidate is None:
                correct_candidate = key
            else:
                # This case should not happen if the problem is well-posed
                return "Error: Multiple candidates fit the data."
        else:
            reasons_for_failure[key] = " and ".join(failures)

    # --- 4. Final Verification ---
    llm_answer = 'A' # The answer provided by the LLM being checked

    if correct_candidate == llm_answer:
        return "Correct"
    else:
        if llm_answer in reasons_for_failure:
            reason = reasons_for_failure[llm_answer]
            return f"Incorrect. The provided answer {llm_answer} ({candidates[llm_answer]['name']}) is wrong because {reason}. The correct answer is {correct_candidate} ({candidates[correct_candidate]['name']})."
        else:
            return f"Incorrect. The correct answer is {correct_candidate}, but the provided answer was {llm_answer}."

# Execute the check
result = check_spectroscopy_answer()
print(result)