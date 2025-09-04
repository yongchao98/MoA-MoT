def check_chemistry_answer():
    """
    This function checks the correctness of the answer for the given analytical chemistry problem.
    It evaluates the proposed answer (4-chlorobenzoic acid) against the spectral data.
    """

    # --- 1. Define Spectral Data Constraints from the Question ---
    spectral_data = {
        'ms': {
            'm_z': 156,
            'm_plus_2_ratio': 0.32,  # Indicates one Cl atom
        },
        'ir': {
            'is_carboxylic_acid': True, # Based on broad 3500-2700 cm-1 and C=O at 1720 cm-1
        },
        'nmr': {
            'has_acid_proton': True, # Based on 11.0 ppm peak
            'aromatic_protons': 4,
            'aromatic_pattern': 'para' # Based on two doublets of 2H each
        }
    }

    # --- 2. Define Properties of Candidate Molecules ---
    molecules = {
        'A': {
            'name': '2-chlorobenzoic acid',
            'mw_35cl': 156,
            'has_cl': 1,
            'is_carboxylic_acid': True,
            'aromatic_protons': 4,
            'aromatic_pattern': 'ortho' # 1,2-disubstituted
        },
        'B': {
            'name': '4-chlorobenzoic acid',
            'mw_35cl': 156,
            'has_cl': 1,
            'is_carboxylic_acid': True,
            'aromatic_protons': 4,
            'aromatic_pattern': 'para' # 1,4-disubstituted
        },
        'C': {
            'name': '3-Chloro-2-hydroxybenzaldehyde',
            'mw_35cl': 156,
            'has_cl': 1,
            'is_carboxylic_acid': False, # It's a phenol and aldehyde
            'aromatic_protons': 3,
            'aromatic_pattern': 'trisubstituted' # 1,2,3-trisubstituted
        },
        'D': {
            'name': 'Phenyl chloroformate',
            'mw_35cl': 156,
            'has_cl': 1,
            'is_carboxylic_acid': False, # It's an ester/acid chloride
            'aromatic_protons': 5,
            'aromatic_pattern': 'monosubstituted'
        }
    }

    # --- 3. Check the Provided Answer (B) ---
    answer_key = 'B'
    candidate = molecules[answer_key]

    # Check MS data
    # All candidates have MW 156 and 1 Cl, so this is a basic check that all pass.
    if candidate['mw_35cl'] != spectral_data['ms']['m_z'] or candidate['has_cl'] != 1:
        return f"Incorrect: The candidate {candidate['name']} does not match the basic MS data (MW and Cl presence)."

    # Check IR data (presence of carboxylic acid)
    if candidate['is_carboxylic_acid'] != spectral_data['ir']['is_carboxylic_acid']:
        return f"Incorrect: The IR spectrum indicates a carboxylic acid, but {candidate['name']} is not a carboxylic acid."

    # Check NMR data (presence of carboxylic acid proton)
    if candidate['is_carboxylic_acid'] != spectral_data['nmr']['has_acid_proton']:
        return f"Incorrect: The NMR spectrum shows a carboxylic acid proton at 11.0 ppm, which is inconsistent with the structure of {candidate['name']}."

    # Check NMR data (number of aromatic protons)
    if candidate['aromatic_protons'] != spectral_data['nmr']['aromatic_protons']:
        return f"Incorrect: The NMR spectrum shows {spectral_data['nmr']['aromatic_protons']} aromatic protons, but {candidate['name']} has {candidate['aromatic_protons']}."

    # Check NMR data (aromatic substitution pattern)
    if candidate['aromatic_pattern'] != spectral_data['nmr']['aromatic_pattern']:
        return f"Incorrect: The NMR spectrum's pattern of two doublets (2H each) indicates a '{spectral_data['nmr']['aromatic_pattern']}' substitution, but {candidate['name']} has an '{candidate['aromatic_pattern']}' substitution pattern."

    # If all checks for the candidate pass, the answer is correct.
    # A more robust check would also verify that all other options fail at least one test.
    # Let's quickly verify that:
    # A fails the NMR pattern check.
    # C fails the IR, NMR acid proton, and NMR aromatic proton count checks.
    # D fails the IR, NMR acid proton, and NMR aromatic proton count checks.
    # Since only B passes all checks, the answer is uniquely determined and correct.
    
    return "Correct"

# Run the check and print the result
result = check_chemistry_answer()
print(result)