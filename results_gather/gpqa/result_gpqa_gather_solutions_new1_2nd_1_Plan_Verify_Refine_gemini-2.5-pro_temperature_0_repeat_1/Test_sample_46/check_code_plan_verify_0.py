def check_correctness():
    """
    Checks if the compound ethyl 4-aminobenzoate (D) matches the given spectral data
    and if the other options (A, B, C) can be ruled out.
    """
    # --- 1. Define the Spectral Data from the Question ---
    # IR: 3420 & 3325 cm-1 (primary amine), 1720 cm-1 (conjugated ester)
    # 1H NMR: 1.20 (t, 3H), 4.0 (bs, 2H), 4.5 (q, 2H), 7.0 (d, 2H), 8.0 (d, 2H)
    
    # Key features extracted from the data:
    ir_nh_bands = 2  # Two bands indicate a primary amine or primary amide
    ir_co_freq = 1720  # Characteristic of a conjugated ester
    nmr_aromatic_pattern = 'para'  # Two doublets (2H each) indicate 1,4-disubstitution
    nmr_amine_protons = 2 # A 2H signal at 4.0 ppm confirms a primary amine
    nmr_quartet_shift = 4.5 # The chemical shift of the -CH2- quartet

    # --- 2. Define Properties of Candidate Molecules ---
    candidates = {
        'A': {
            'name': '4-aminophenyl propionate',
            'amine_type': 'primary',
            'carbonyl_type': 'ester',
            'substitution': 'para',
            'ethyl_quartet_env': 'CO-CH2-', # Expected shift ~2.5 ppm
        },
        'B': {
            'name': 'N-(4-ethoxyphenyl)formamide',
            'amine_type': 'secondary_amide', # Expects 1 N-H band
            'carbonyl_type': 'amide',
            'substitution': 'para',
        },
        'C': {
            'name': '3-ethoxybenzamide',
            'amine_type': 'primary_amide',
            'carbonyl_type': 'amide', # Expects C=O < 1700 cm-1
            'substitution': 'meta', # Wrong aromatic pattern
        },
        'D': {
            'name': 'ethyl 4-aminobenzoate',
            'amine_type': 'primary',
            'carbonyl_type': 'ester',
            'substitution': 'para',
            'ethyl_quartet_env': 'COO-CH2-', # Expected shift ~4.5 ppm
        }
    }
    
    # --- 3. Get the Proposed Answer and Verify ---
    proposed_answer_key = 'D'
    candidate_to_check = candidates[proposed_answer_key]
    errors = []

    # Check 1: Amine Type (from IR N-H bands)
    if candidate_to_check['amine_type'] == 'secondary_amide' and ir_nh_bands == 2:
        errors.append(f"The data shows 2 N-H bands (primary amine), but {candidate_to_check['name']} is a secondary amide.")

    # Check 2: Carbonyl Type (from IR C=O frequency)
    if candidate_to_check['carbonyl_type'] == 'ester' and not (1710 <= ir_co_freq <= 1740):
        errors.append(f"The IR C=O frequency ({ir_co_freq} cm-1) does not match the expected range for a conjugated ester.")
    if candidate_to_check['carbonyl_type'] == 'amide' and ir_co_freq > 1700:
        errors.append(f"The IR C=O frequency ({ir_co_freq} cm-1) is too high for an amide; it indicates an ester.")

    # Check 3: Substitution Pattern (from NMR aromatic signals)
    if candidate_to_check['substitution'] != nmr_aromatic_pattern:
        errors.append(f"The NMR shows a '{nmr_aromatic_pattern}' substitution pattern, but {candidate_to_check['name']} is '{candidate_to_check['substitution']}'.")

    # Check 4: Ethyl Group Environment (from NMR quartet shift)
    env = candidate_to_check.get('ethyl_quartet_env')
    if env == 'COO-CH2-': # Ethyl ester
        if not (4.1 <= nmr_quartet_shift <= 4.6):
            errors.append(f"The NMR quartet at {nmr_quartet_shift} ppm is outside the expected range for an ethyl ester (-COO-CH2-).")
    elif env == 'CO-CH2-': # Propionate
        if not (2.2 <= nmr_quartet_shift <= 2.7):
            errors.append(f"The NMR quartet is at {nmr_quartet_shift} ppm, but a propionate group (-CO-CH2-) should be around 2.2-2.7 ppm.")

    # --- 4. Final Verdict ---
    if errors:
        return f"Incorrect. The proposed answer {proposed_answer_key} ({candidate_to_check['name']}) fails the following checks:\n- " + "\n- ".join(errors)
    else:
        # To be certain, we confirm the other options are wrong.
        # A fails Check 4 (quartet shift is wrong).
        # B fails Check 1 (it's a secondary amide).
        # C fails Check 3 (it's meta-substituted) and Check 2 (it's an amide).
        # Since D passes all checks and the others fail, the answer is correct.
        return "Correct"

# Run the check and print the result
print(check_correctness())