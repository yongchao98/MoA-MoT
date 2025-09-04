import re

def check_nmr_answer():
    """
    Checks the correctness of the LLM's answer for the NMR spectroscopy question.
    """
    question_constraints = {
        "ring_type": "di-substituted 6-membered aromatic",
        "total_carbons": 8,
        "functional_groups": ["carbonyl", "aromatic_halogen"]
    }

    # The options as presented in the prompt for the final answer being checked.
    options = {
        "A": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "B": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "C": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "D": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)"
    }
    
    llm_answer = "A"

    # --- Step 1: Deduce the expected structure and NMR features from constraints ---
    # Aromatic ring (C6) + 2 substituents = 8 carbons total.
    # Substituents must have 2 carbons.
    # One substituent is a halogen (-X), which has 0 carbons.
    # The other must have 2 carbons and a carbonyl (C=O). This is an acetyl group (-COCH3).
    # The structure is a haloacetophenone (X-C6H4-COCH3).
    # Total protons = 4 on the ring + 3 on the methyl group = 7.

    # Expected NMR for the most likely isomer (para-haloacetophenone):
    expected_features = {
        "aromatic_protons": {
            "count": 4,
            "pattern": "two doublets",
            "shift_range": (6.5, 8.5)
        },
        "acetyl_protons": {
            "count": 3,
            "pattern": "singlet",
            "shift_range": (2.0, 2.7)
        }
    }

    # --- Step 2: Analyze the NMR data selected by the LLM ---
    selected_data_str = options.get(llm_answer)
    if not selected_data_str:
        return f"Invalid answer key '{llm_answer}'. The key must be one of {list(options.keys())}."

    # Use regex to parse the NMR data string
    # Pattern to find: shift (protonsH, multiplicity)
    signals = re.findall(r"(\d+\.\d+)\s*\((\d+)H,\s*([a-z]+)\)", selected_data_str)
    
    if not signals:
        return f"Could not parse the NMR data for option {llm_answer}: '{selected_data_str}'"

    parsed_signals = [{'shift': float(s[0]), 'protons': int(s[1]), 'multiplicity': s[2]} for s in signals]

    # --- Step 3: Check if the selected data satisfies all constraints ---
    
    # Check 1: Acetyl group signal
    acetyl_signal = [s for s in parsed_signals if s['protons'] == 3 and s['multiplicity'] == 's']
    if not acetyl_signal:
        return f"Incorrect. The selected data for answer {llm_answer} does not contain a singlet for 3 protons (3H, s), which is required for the acetyl group (-COCH3)."
    
    acetyl_shift = acetyl_signal[0]['shift']
    if not (expected_features["acetyl_protons"]["shift_range"][0] <= acetyl_shift <= expected_features["acetyl_protons"]["shift_range"][1]):
        return f"Incorrect. The methyl singlet at {acetyl_shift} ppm is outside the expected range of {expected_features['acetyl_protons']['shift_range']} ppm for an acetyl group."

    # Check 2: Aromatic signals
    aromatic_signals = [s for s in parsed_signals if s['shift'] >= expected_features["aromatic_protons"]["shift_range"][0]]
    total_aromatic_protons = sum(s['protons'] for s in aromatic_signals)
    
    if total_aromatic_protons != 4:
        return f"Incorrect. The selected data for answer {llm_answer} has {total_aromatic_protons} aromatic protons, but a di-substituted benzene ring must have 4."

    # Check 3: Para-substitution pattern
    doublet_signals = [s for s in aromatic_signals if s['multiplicity'] == 'd' and s['protons'] == 2]
    if len(doublet_signals) != 2:
        return f"Incorrect. The aromatic signals in answer {llm_answer} do not match the expected pattern of two doublets for 2H each, which is characteristic of a para-substituted ring."

    # If all checks pass for the selected answer, it is correct.
    return "Correct"

# Run the check
result = check_nmr_answer()
print(result)