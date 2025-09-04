def check_answer():
    """
    This function checks the correctness of the provided answer by systematically applying
    the constraints from the spectral data to the candidate molecules.
    """
    # Define the properties of each candidate molecule based on the question's options
    candidates = {
        'A': {
            'name': '3-Chloro-2-hydroxybenzaldehyde',
            'mw_35cl': 156,
            'has_one_cl': True,
            'functional_group': 'aldehyde/phenol',
            'nmr_pattern': '1,2,3-trisubstituted'
        },
        'B': {
            'name': '2-chlorobenzoic acid',
            'mw_35cl': 156,
            'has_one_cl': True,
            'functional_group': 'carboxylic acid',
            'nmr_pattern': 'ortho'
        },
        'C': {
            'name': '4-chlorobenzoic acid',
            'mw_35cl': 156,
            'has_one_cl': True,
            'functional_group': 'carboxylic acid',
            'nmr_pattern': 'para'
        },
        'D': {
            'name': 'Phenyl chloroformate',
            'mw_35cl': 156,
            'has_one_cl': True,
            'functional_group': 'chloroformate',
            'nmr_pattern': 'monosubstituted'
        }
    }
    
    provided_answer = 'C'
    
    # --- Step 1: Mass Spectrometry Check ---
    # Data: m/z = 156 (M+), m/z = 158 (M+2, 32%)
    # Interpretation: MW with 35Cl is 156, and one Cl atom is present.
    survivors_ms = {key: props for key, props in candidates.items() if props['mw_35cl'] == 156 and props['has_one_cl']}
    if len(survivors_ms) != 4:
        return "MS Check Failed: Not all candidates match the basic mass spec data, which is unexpected."

    # --- Step 2: IR Spectroscopy Check ---
    # Data: broad peak 3500-2700 cm^-1, strong sharp peak at 1720 cm-1
    # Interpretation: Presence of a carboxylic acid (-COOH) group.
    survivors_ir = {key: props for key, props in survivors_ms.items() if props['functional_group'] == 'carboxylic acid'}
    
    # --- Step 3: 1H NMR Spectroscopy Check ---
    # Data: 11.0 ppm (s, 1H), 8.02 ppm (d, 2H), 7.72 (d, 2H)
    # Interpretation: Carboxylic acid proton + para-disubstituted aromatic ring.
    survivors_nmr = {key: props for key, props in survivors_ir.items() if props['nmr_pattern'] == 'para'}

    # --- Final Verification ---
    if len(survivors_nmr) == 1:
        correct_key = list(survivors_nmr.keys())[0]
        if correct_key == provided_answer:
            return "Correct"
        else:
            return (f"Incorrect. The provided answer is '{provided_answer}', but the analysis shows the correct answer is '{correct_key}'.\n"
                    f"Reasoning: Only candidate {correct_key} ({candidates[correct_key]['name']}) is a carboxylic acid (satisfies IR) "
                    f"and has a para-substitution pattern (satisfies NMR).")
    elif len(survivors_nmr) == 0:
        return "Incorrect. The analysis shows that no candidate satisfies all the spectral data."
    else:
        return "Incorrect. The analysis is inconclusive as multiple candidates satisfy all conditions."

# Run the check and print the result
result = check_answer()
print(result)