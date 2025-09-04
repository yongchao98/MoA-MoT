import math

def check_answer():
    """
    This function checks the correctness of the provided LLM answer by systematically applying
    the spectral data constraints to the given candidate molecules.
    """

    # Define the properties of each candidate molecule based on chemical knowledge.
    # Using integer masses for nominal mass spectrometry (m/z).
    candidates = {
        'A': {
            'name': '3-Chloro-2-hydroxybenzaldehyde',
            'mw_cl35': 7*12 + 5*1 + 35 + 2*16,  # 156
            'has_cl': True,
            'is_carboxylic_acid': False,
            'substitution_pattern': '1,2,3-trisubstituted',
            'expected_aromatic_nmr': 'complex multiplet'
        },
        'B': {
            'name': 'Phenyl chloroformate',
            'mw_cl35': 7*12 + 5*1 + 35 + 2*16,  # 156
            'has_cl': True,
            'is_carboxylic_acid': False,
            'substitution_pattern': 'monosubstituted',
            'expected_aromatic_nmr': 'complex multiplet'
        },
        'C': {
            'name': '2-chlorobenzoic acid',
            'mw_cl35': 7*12 + 5*1 + 35 + 2*16,  # 156
            'has_cl': True,
            'is_carboxylic_acid': True,
            'substitution_pattern': 'ortho',
            'expected_aromatic_nmr': 'complex multiplet'
        },
        'D': {
            'name': '4-chlorobenzoic acid',
            'mw_cl35': 7*12 + 5*1 + 35 + 2*16,  # 156
            'has_cl': True,
            'is_carboxylic_acid': True,
            'substitution_pattern': 'para',
            'expected_aromatic_nmr': 'two doublets, 2H each'
        }
    }

    # The answer provided by the LLM to be checked.
    llm_answer = 'D'

    # --- Constraint 1: Mass Spectrometry ---
    # Molecular ion at m/z = 156 and M+2 peak at m/z=158 (~32%) indicates MW=156 and one Cl atom.
    passing_ms = []
    for key, props in candidates.items():
        if props['mw_cl35'] == 156 and props['has_cl']:
            passing_ms.append(key)
    
    # All candidates pass this constraint, so we just check if the LLM's answer is in the list.
    if llm_answer not in passing_ms:
        return f"The answer '{llm_answer}' is incorrect. It fails the Mass Spec constraint. Its molecular weight (using ³⁵Cl) is not 156 or it does not contain a chlorine atom."

    # --- Constraint 2: IR Spectroscopy ---
    # Broad peak 3500-2700 cm⁻¹ and strong sharp peak at 1720 cm⁻¹ strongly indicate a carboxylic acid.
    passing_ir = []
    for key in passing_ms:
        if candidates[key]['is_carboxylic_acid']:
            passing_ir.append(key)

    if llm_answer not in passing_ir:
        return f"The answer '{llm_answer}' is incorrect. It fails the IR Spectroscopy constraint. The IR data indicates a carboxylic acid, but candidate {llm_answer} is not a carboxylic acid."

    # --- Constraint 3: 1H NMR Spectroscopy ---
    # - 11.0 ppm (s, 1H) confirms the carboxylic acid proton.
    # - 8.02 ppm (d, 2H) and 7.72 (d, 2H) is a classic pattern for a para-substituted (1,4-disubstituted) benzene ring.
    passing_nmr = []
    for key in passing_ir:
        # Check for both carboxylic acid presence (already filtered by IR) and the correct substitution pattern.
        if candidates[key]['substitution_pattern'] == 'para':
            passing_nmr.append(key)

    if llm_answer not in passing_nmr:
        props = candidates[llm_answer]
        return f"The answer '{llm_answer}' is incorrect. It fails the 1H NMR constraint. Its '{props['substitution_pattern']}' substitution pattern would produce a complex multiplet, not the two doublets (2H each) observed, which is characteristic of a 'para' substitution."

    # --- Final Conclusion ---
    # Check if the LLM's answer is the single candidate that passed all filters.
    if len(passing_nmr) == 1 and passing_nmr[0] == llm_answer:
        return "Correct"
    elif len(passing_nmr) == 0:
        return f"The answer '{llm_answer}' is incorrect, and no other candidate fits all criteria."
    else:
        correct_answer = passing_nmr[0]
        return f"The answer '{llm_answer}' is incorrect. The only candidate that satisfies all spectral data is '{correct_answer}' ({candidates[correct_answer]['name']})."

# Execute the check and print the result.
result = check_answer()
print(result)