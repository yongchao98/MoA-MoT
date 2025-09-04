def check_spectroscopy_answer():
    """
    This function checks if the proposed compound, ethyl 4-aminobenzoate,
    matches the provided IR and 1H NMR spectral data.
    """
    
    # --- Provided Data from the Question ---
    # IR data
    ir_nh_stretches = [3420, 3325]
    ir_co_stretch = 1720
    
    # 1H NMR data
    nmr_signals = {
        1.20: {'multiplicity': 't', 'integration': 3, 'assignment': 'ethyl_CH3'},
        4.0: {'multiplicity': 'bs', 'integration': 2, 'assignment': 'amine_NH2'},
        4.5: {'multiplicity': 'q', 'integration': 2, 'assignment': 'ethyl_ester_OCH2'},
        7.0: {'multiplicity': 'd', 'integration': 2, 'assignment': 'aromatic_H_ortho_to_NH2'},
        8.0: {'multiplicity': 'd', 'integration': 2, 'assignment': 'aromatic_H_ortho_to_COOR'}
    }
    
    # --- Properties of the Proposed Answer: C) ethyl 4-aminobenzoate ---
    proposed_compound = {
        'name': 'ethyl 4-aminobenzoate',
        'amine_type': 'primary',  # Has -NH2 group
        'carbonyl_type': 'conjugated_ester', # C=O stretch ~1715-1730 cm-1
        'aromatic_substitution': 'para', # 1,4-disubstituted
        'nmr_fragments': ['ethyl_ester', 'primary_amine', 'para_aromatic_ring']
    }
    
    # --- List to store reasons for incorrectness ---
    errors = []

    # --- Verification Steps ---

    # 1. Check IR Data
    # The two bands at 3420 and 3325 cm-1 indicate a primary amine (-NH2).
    if proposed_compound['amine_type'] != 'primary':
        errors.append(f"IR Mismatch: The data shows two N-H bands ({ir_nh_stretches} cm-1), characteristic of a primary amine. The proposed compound is not a primary amine.")
        
    # The strong band at 1720 cm-1 indicates a conjugated ester. Amides are typically < 1700 cm-1.
    if proposed_compound['carbonyl_type'] != 'conjugated_ester':
        errors.append(f"IR Mismatch: The C=O band at {ir_co_stretch} cm-1 is characteristic of a conjugated ester. The proposed compound does not fit this description.")

    # 2. Check 1H NMR Data
    # Check for ethyl ester pattern: triplet (3H) ~1.2 ppm and quartet (2H) ~4.5 ppm.
    has_ethyl_ch3 = any(s['multiplicity'] == 't' and s['integration'] == 3 for s in nmr_signals.values())
    has_ester_och2 = any(s['multiplicity'] == 'q' and s['integration'] == 2 and k > 4.2 for k, s in nmr_signals.items()) # High shift is key
    if not (has_ethyl_ch3 and has_ester_och2 and 'ethyl_ester' in proposed_compound['nmr_fragments']):
        errors.append("NMR Mismatch: The data's triplet at 1.20 ppm and quartet at 4.5 ppm indicate an ethyl ester group (-COOCH2CH3), which is inconsistent with the proposed structure.")

    # Check for primary amine protons: broad singlet (2H) ~4.0 ppm.
    has_amine_protons = any(s['multiplicity'] == 'bs' and s['integration'] == 2 for s in nmr_signals.values())
    if not (has_amine_protons and 'primary_amine' in proposed_compound['nmr_fragments']):
        errors.append("NMR Mismatch: The data's broad singlet at 4.0 ppm (2H) indicates a primary amine, which is inconsistent with the proposed structure.")

    # Check for para-substitution pattern: two doublets, each integrating to 2H.
    doublet_count = sum(1 for s in nmr_signals.values() if s['multiplicity'] == 'd' and s['integration'] == 2)
    if not (doublet_count == 2 and proposed_compound['aromatic_substitution'] == 'para'):
        errors.append("NMR Mismatch: The data's two doublets in the aromatic region indicate a para-substituted ring, which is inconsistent with the proposed structure.")

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. The proposed answer does not satisfy the following constraints:\n- " + "\n- ".join(errors)

# Execute the check and print the result
result = check_spectroscopy_answer()
print(result)