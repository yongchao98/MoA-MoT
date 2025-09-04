def check_correctness():
    """
    This function checks the correctness of the provided answer based on the spectral data.
    The question asks to identify a compound from MS, IR, and 1H NMR data.
    The provided final answer from the LLM is 'D', which corresponds to 4-chlorobenzoic acid.
    """

    # --- Define Constraints from Spectral Data ---
    # MS: MW=156 (nominal), 1 Cl atom
    # IR: Carboxylic acid group (-COOH)
    # NMR: Carboxylic acid proton, 4 aromatic protons, para-substitution pattern

    # --- Define Properties of Candidate Molecules ---
    # Note: The lettering (A, B, C, D) is based on the question's list.
    candidates = {
        'A': {
            'name': 'Phenyl chloroformate',
            'nominal_mass': 156,
            'has_cl': True,
            'has_carboxylic_acid': False,
            'aromatic_protons': 5,
            'substitution': 'mono'
        },
        'B': {
            'name': '3-Chloro-2-hydroxybenzaldehyde',
            'nominal_mass': 156,
            'has_cl': True,
            'has_carboxylic_acid': False,  # It has phenol and aldehyde groups
            'aromatic_protons': 3,
            'substitution': '1,2,3-tri'
        },
        'C': {
            'name': '2-chlorobenzoic acid',
            'nominal_mass': 156,
            'has_cl': True,
            'has_carboxylic_acid': True,
            'aromatic_protons': 4,
            'substitution': 'ortho'  # 1,2-disubstituted
        },
        'D': {
            'name': '4-chlorobenzoic acid',
            'nominal_mass': 156,
            'has_cl': True,
            'has_carboxylic_acid': True,
            'aromatic_protons': 4,
            'substitution': 'para'  # 1,4-disubstituted
        }
    }

    # The final answer provided in the prompt to be checked
    answer_to_check = 'D'

    # Get the properties of the chosen answer
    chosen_molecule = candidates[answer_to_check]

    # --- Check against constraints ---

    # 1. Check MS data constraints
    if not chosen_molecule['has_cl'] or chosen_molecule['nominal_mass'] != 156:
        return f"MS constraint not satisfied: The MS data (m/z=156, M+2 peak) indicates a molecule with nominal mass 156 containing one chlorine atom. This is not met by {chosen_molecule['name']}."

    # 2. Check IR data constraints
    if not chosen_molecule['has_carboxylic_acid']:
        return f"IR constraint not satisfied: The IR spectrum (broad peak 3500-2700 cm-1, sharp peak 1720 cm-1) indicates a carboxylic acid functional group, which {chosen_molecule['name']} lacks."

    # 3. Check 1H NMR data constraints
    if chosen_molecule['aromatic_protons'] != 4:
        return f"NMR constraint not satisfied: The NMR data shows 4 aromatic protons (a 2H doublet and another 2H doublet). {chosen_molecule['name']} has {chosen_molecule['aromatic_protons']} aromatic protons."
    
    if chosen_molecule['substitution'] != 'para':
        return f"NMR constraint not satisfied: The NMR shows a pattern of two doublets (2H each), which is characteristic of a para-substituted (1,4) benzene ring. {chosen_molecule['name']} has '{chosen_molecule['substitution']}' substitution, which would produce a different splitting pattern."

    # If all constraints are satisfied
    return "Correct"

# Execute the check and print the result
print(check_correctness())