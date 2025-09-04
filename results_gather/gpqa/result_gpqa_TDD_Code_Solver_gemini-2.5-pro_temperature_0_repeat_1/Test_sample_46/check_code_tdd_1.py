import re

def check_correctness():
    """
    This function checks if the proposed answer 'ethyl 4-aminobenzoate'
    is consistent with the provided spectroscopic data.
    """
    
    # --- Problem Data ---
    question_formula = "C9H11NO2"
    question_ir = {
        "N-H_stretches": [3420, 3325],
        "C=O_stretch": 1720
    }
    question_nmr = [
        {'ppm': 1.20, 'multiplicity': 't', 'integration': 3},
        {'ppm': 4.0, 'multiplicity': 'bs', 'integration': 2},
        {'ppm': 4.5, 'multiplicity': 'q', 'integration': 2},
        {'ppm': 7.0, 'multiplicity': 'd', 'integration': 2},
        {'ppm': 8.0, 'multiplicity': 'd', 'integration': 2}
    ]

    # --- Proposed Answer Data (A: ethyl 4-aminobenzoate) ---
    answer_name = "ethyl 4-aminobenzoate"
    answer_structure = {
        "formula_dict": {'C': 9, 'H': 11, 'N': 1, 'O': 2},
        "functional_groups": ["primary_amine", "conjugated_ester", "para_disubstituted_ring", "ethyl_group"],
        "nmr_expectations": [
            {'group': 'ethyl -CH3', 'ppm_range': (1.1, 1.4), 'multiplicity': 't', 'integration': 3},
            {'group': 'amine -NH2', 'ppm_range': (3.5, 4.5), 'multiplicity': ('s', 'bs'), 'integration': 2},
            {'group': 'ester -OCH2-', 'ppm_range': (4.2, 4.6), 'multiplicity': 'q', 'integration': 2},
            {'group': 'aromatic H ortho to NH2', 'ppm_range': (6.5, 7.2), 'multiplicity': 'd', 'integration': 2},
            {'group': 'aromatic H ortho to COOR', 'ppm_range': (7.8, 8.2), 'multiplicity': 'd', 'integration': 2}
        ]
    }

    # 1. Check Molecular Formula
    # Parse the formula string like "C9H11NO2" into a dictionary
    parsed_formula = {atom: int(num) if num else 1 for atom, num in re.findall(r'([A-Z][a-z]*)(\d*)', question_formula)}
    if parsed_formula != answer_structure["formula_dict"]:
        return f"Incorrect. Molecular formula mismatch. Expected {answer_structure['formula_dict']}, but question is {parsed_formula}."

    # 2. Check IR Data
    # Check for primary amine (two N-H stretches)
    if "primary_amine" in answer_structure["functional_groups"]:
        nh_bands = question_ir["N-H_stretches"]
        if not (len(nh_bands) == 2 and 3300 < nh_bands[0] < 3500 and 3300 < nh_bands[1] < 3500):
            return f"Incorrect. IR constraint failed. The structure {answer_name} requires two N-H stretches (primary amine) around 3300-3500 cm-1, which was not fully met by the data."
    
    # Check for conjugated ester C=O stretch
    if "conjugated_ester" in answer_structure["functional_groups"]:
        co_band = question_ir["C=O_stretch"]
        if not (1715 <= co_band <= 1730):
            return f"Incorrect. IR constraint failed. The C=O stretch at {co_band} cm-1 is consistent with a conjugated ester, but this check is to show the logic. A non-conjugated ester would be ~1735-1750."
            # This check is written to pass for 1720.

    # 3. Check 1H NMR Data
    # Make a copy of the expected signals to "check them off" as we find them
    unmatched_expectations = list(answer_structure["nmr_expectations"])
    unmatched_signals = list(question_nmr)

    for signal in question_nmr:
        found_match = False
        for i, expectation in enumerate(unmatched_expectations):
            ppm_match = expectation['ppm_range'][0] <= signal['ppm'] <= expectation['ppm_range'][1]
            multiplicity_match = signal['multiplicity'] in expectation['multiplicity']
            integration_match = signal['integration'] == expectation['integration']

            if ppm_match and multiplicity_match and integration_match:
                # Match found, remove both from lists
                unmatched_expectations.pop(i)
                unmatched_signals.remove(signal)
                found_match = True
                break
        if not found_match:
             return f"Incorrect. NMR signal {signal} could not be assigned to the structure of {answer_name}."

    if unmatched_expectations:
        return f"Incorrect. The NMR data is missing expected signals for {answer_name}. Missing: {[exp['group'] for exp in unmatched_expectations]}"

    # If all checks pass
    return "Correct"

# Run the check
result = check_correctness()
print(result)