import math

def check_correctness_of_answer():
    """
    This function programmatically checks if the proposed answer, 4-chlorobenzoic acid,
    is consistent with all the spectral data provided in the question.
    """
    
    # Part 1: Define the constraints from the spectral data provided in the question.
    
    # Mass Spectrometry Data
    ms_mz_m_ion = 156
    ms_mz_m_plus_2 = 158
    ms_ratio = 100 / 32.0

    # IR Data
    ir_carbonyl_peak = 1720  # cm^-1
    # The presence of a broad peak from 3500-2700 cm^-1 is a qualitative feature.

    # 1H NMR Data
    nmr_signals = [
        {'shift': 11.0, 'multiplicity': 's', 'integration': 1},
        {'shift': 8.02, 'multiplicity': 'd', 'integration': 2},
        {'shift': 7.72, 'multiplicity': 'd', 'integration': 2}
    ]

    # Part 2: Define the expected properties of the proposed answer, 4-chlorobenzoic acid.
    
    # Molecular Formula: C7H5ClO2
    # Using integer atomic masses for m/z calculation
    expected_mw_m_ion = (7 * 12) + (5 * 1) + (1 * 35) + (2 * 16)
    expected_mw_m_plus_2 = (7 * 12) + (5 * 1) + (1 * 37) + (2 * 16)
    
    # The natural abundance ratio of 35Cl to 37Cl is ~3:1.
    # The given ratio of 100:32 (3.125) is a perfect match.
    expected_ms_ratio_range = (3.0, 3.3) 

    # Expected IR features for a conjugated carboxylic acid
    expected_ir_carbonyl_range = (1680, 1740)
    # A broad O-H stretch from ~3300-2500 cm^-1 is expected.

    # Expected 1H NMR features for 4-chlorobenzoic acid (a para-substituted ring)
    # - One acidic proton (COOH), singlet, chemical shift > 10 ppm.
    # - Two sets of aromatic protons due to symmetry.
    # - Each set contains 2 protons and splits the other into a doublet.
    
    # Part 3: Perform the checks and log any errors.
    
    error_log = []

    # Check 1: Mass Spectrometry
    if expected_mw_m_ion != ms_mz_m_ion:
        error_log.append(f"Mass Spec M+ peak mismatch: Calculated MW for C7H5(35)ClO2 is {expected_mw_m_ion}, but observed m/z is {ms_mz_m_ion}.")
    if expected_mw_m_plus_2 != ms_mz_m_plus_2:
        error_log.append(f"Mass Spec M+2 peak mismatch: Calculated MW for C7H5(37)ClO2 is {expected_mw_m_plus_2}, but observed m/z is {ms_mz_m_plus_2}.")
    if not (expected_ms_ratio_range[0] <= ms_ratio <= expected_ms_ratio_range[1]):
        error_log.append(f"Mass Spec isotope ratio mismatch: Observed ratio is {ms_ratio:.2f}, which is outside the expected range {expected_ms_ratio_range} for one chlorine atom.")

    # Check 2: IR Spectroscopy
    if not (expected_ir_carbonyl_range[0] <= ir_carbonyl_peak <= expected_ir_carbonyl_range[1]):
        error_log.append(f"IR mismatch: Carbonyl peak at {ir_carbonyl_peak} cm-1 is outside the expected range {expected_ir_carbonyl_range} for a conjugated carboxylic acid.")
    # The broad peak from 3500-2700 cm-1 is qualitatively consistent with a carboxylic acid O-H stretch. This check passes.

    # Check 3: 1H NMR Spectroscopy
    # Check for the acidic proton
    acid_proton_found = any(s['shift'] >= 10 and s['multiplicity'] == 's' and s['integration'] == 1 for s in nmr_signals)
    if not acid_proton_found:
        error_log.append("NMR mismatch: The characteristic carboxylic acid proton (1H, singlet, >10 ppm) was not found in the provided data.")
    
    # Check for the para-substitution pattern in the aromatic region
    aromatic_signals = [s for s in nmr_signals if 6.5 < s['shift'] < 9.0]
    if len(aromatic_signals) != 2:
        error_log.append(f"NMR mismatch: Expected 2 aromatic signals for a para-substituted ring, but found {len(aromatic_signals)}.")
    else:
        are_doublets = all(s['multiplicity'] == 'd' for s in aromatic_signals)
        are_2H_integrals = all(s['integration'] == 2 for s in aromatic_signals)
        if not (are_doublets and are_2H_integrals):
            error_log.append("NMR mismatch: The aromatic signals do not match the expected pattern of two doublets, each integrating to 2H, which is characteristic of para-substitution.")

    # Part 4: Return the final verdict.
    if not error_log:
        return "Correct"
    else:
        error_header = "The provided answer is incorrect for the following reason(s):"
        return error_header + "\n- " + "\n- ".join(error_log)

# Run the check and print the result.
result = check_correctness_of_answer()
print(result)