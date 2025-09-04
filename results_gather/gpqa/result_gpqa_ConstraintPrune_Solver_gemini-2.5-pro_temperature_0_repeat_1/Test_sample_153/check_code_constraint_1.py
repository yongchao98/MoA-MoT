def check_answer():
    """
    This function checks if the proposed structure, 4-chlorobenzoic acid,
    is consistent with the provided spectral data.
    """
    
    # --- 1. Define Spectral Data from the Question ---
    MS_DATA = {'m_plus': 156, 'm_plus_2': 158, 'm_plus_2_intensity_percent': 32}
    IR_DATA = {'broad_peak_range': (2700, 3500), 'sharp_peak': 1720}
    NMR_DATA = [
        {'ppm': 11.0, 'type': 's', 'integration': 1},
        {'ppm': 8.02, 'type': 'd', 'integration': 2},
        {'ppm': 7.72, 'type': 'd', 'integration': 2}
    ]
    
    # --- 2. Define Properties of the Proposed Answer (D) 4-chlorobenzoic acid ---
    PROPOSED_ANSWER = {
        'name': "4-chlorobenzoic acid",
        'formula': "C7H5ClO2",
        'functional_groups': ["carboxylic_acid", "aryl_halide"],
        'nmr_symmetry': "para-disubstituted"
    }

    # --- 3. Perform Checks ---

    # Check 3a: Mass Spectrometry
    # Calculate expected MW for M+ (with 35Cl) and M+2 (with 37Cl)
    expected_m_plus = (7 * 12) + (5 * 1) + 35 + (2 * 16)
    expected_m_plus_2 = (7 * 12) + (5 * 1) + 37 + (2 * 16)
    
    if expected_m_plus != MS_DATA['m_plus']:
        return f"Incorrect: Mass Spec check failed. The calculated M+ peak for 4-chlorobenzoic acid is {expected_m_plus}, but the data shows {MS_DATA['m_plus']}."
    
    if expected_m_plus_2 != MS_DATA['m_plus_2']:
        return f"Incorrect: Mass Spec check failed. The calculated M+2 peak for 4-chlorobenzoic acid is {expected_m_plus_2}, but the data shows {MS_DATA['m_plus_2']}."

    # Check M+2 isotope ratio for one chlorine atom (approx. 3:1 or 33%)
    intensity_ratio = MS_DATA['m_plus_2_intensity_percent']
    if not (28 <= intensity_ratio <= 38):
        return f"Incorrect: Mass Spec check failed. The M+2 peak intensity ({intensity_ratio}%) is not consistent with the presence of one chlorine atom (expected ~33%)."

    # Check 3b: Infrared (IR) Spectroscopy
    # The combination of a very broad O-H (3500-2700 cm-1) and a C=O (~1720 cm-1) is classic for a carboxylic acid.
    if "carboxylic_acid" not in PROPOSED_ANSWER['functional_groups']:
        return "Incorrect: IR check failed. The data strongly indicates a carboxylic acid, but the proposed molecule is not one."
    
    c_double_o_peak = IR_DATA['sharp_peak']
    if not (1700 <= c_double_o_peak <= 1725):
        # This range is typical for aromatic carboxylic acids. 1720 is a perfect fit.
        return f"Incorrect: IR check failed. The C=O peak at {c_double_o_peak} cm-1 is not in the expected range for an aromatic carboxylic acid."

    # Check 3c: 1H NMR Spectroscopy
    # Check for carboxylic acid proton
    acidic_proton_signal = next((s for s in NMR_DATA if 10 <= s['ppm'] <= 13 and s['integration'] == 1 and s['type'] == 's'), None)
    if not acidic_proton_signal:
        return "Incorrect: NMR check failed. The proposed molecule is a carboxylic acid, but the characteristic acidic proton signal (10-13 ppm, 1H, singlet) is not found in the data."

    # Check aromatic region for para-substitution pattern
    aromatic_signals = [s for s in NMR_DATA if 6.5 <= s['ppm'] <= 8.5]
    if PROPOSED_ANSWER['nmr_symmetry'] == "para-disubstituted":
        # Expect two doublets, each integrating to 2H
        doublets_2H = [s for s in aromatic_signals if s['type'] == 'd' and s['integration'] == 2]
        if len(doublets_2H) != 2:
            return "Incorrect: NMR check failed. The aromatic region pattern does not match the expected two doublets (each 2H) for a para-substituted ring."
    else:
        return "Incorrect: NMR check failed. The proposed structure's symmetry does not match the observed NMR pattern."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)