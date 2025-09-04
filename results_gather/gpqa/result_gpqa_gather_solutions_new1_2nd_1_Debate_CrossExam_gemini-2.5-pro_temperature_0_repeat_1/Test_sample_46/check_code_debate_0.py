import re

def check_spectroscopy_answer():
    """
    Checks the correctness of the answer for a spectroscopy problem.

    The function defines the problem's data and the properties of each candidate molecule.
    It then systematically checks each candidate against the data (Molecular Formula, IR, NMR)
    to find the one that matches all constraints. Finally, it compares this result
    with the provided answer.
    """
    # --- Problem Definition ---
    question_data = {
        "molecular_formula": "C9H11NO2",
        "ir_peaks": {
            "N-H": [3420, 3325],  # Two peaks indicate primary amine
            "C=O": 1720           # Ester range
        },
        "nmr_signals": [
            {"ppm": 1.20, "type": "t", "H": 3}, # -CH3 of ethyl
            {"ppm": 4.0, "type": "bs", "H": 2}, # -NH2
            {"ppm": 4.5, "type": "q", "H": 2}, # -O-CH2- of ethyl ester
            {"ppm": 7.0, "type": "d", "H": 2}, # Aromatic
            {"ppm": 8.0, "type": "d", "H": 2}  # Aromatic
        ]
    }
    
    # --- Candidate Definitions ---
    # Define expected properties for each candidate molecule.
    candidates = {
        "A": {
            "name": "3-ethoxybenzamide",
            "formula": "C9H11NO2",
            "features": {
                "amine_type": "primary_amide", # Has -NH2, but is an amide
                "carbonyl_type": "amide",      # C=O stretch ~1680 cm-1
                "substitution": "meta",        # Complex aromatic NMR, not 2 doublets
                "ethyl_type": "ethoxy"         # -O-CH2- quartet ~4.0 ppm
            }
        },
        "B": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "formula": "C9H11NO2",
            "features": {
                "amine_type": "secondary_amide", # Only one N-H peak in IR
                "carbonyl_type": "amide",
                "substitution": "para",
                "ethyl_type": "ethoxy"
            }
        },
        "C": {
            "name": "4-aminophenyl propionate",
            "formula": "C9H11NO2",
            "features": {
                "amine_type": "primary_amine",
                "carbonyl_type": "ester",
                "substitution": "para",
                "ethyl_type": "propionate" # -CO-CH2- quartet ~2.5 ppm
            }
        },
        "D": {
            "name": "ethyl 4-aminobenzoate",
            "formula": "C9H11NO2",
            "features": {
                "amine_type": "primary_amine", # Two N-H peaks in IR
                "carbonyl_type": "ester",      # C=O stretch ~1720 cm-1
                "substitution": "para",        # Two aromatic doublets
                "ethyl_type": "ethyl_ester"    # -O-CH2- quartet ~4.5 ppm
            }
        }
    }
    
    provided_answer = "D"

    # --- Analysis Functions ---
    def get_formula_parts(formula_str):
        parts = re.findall(r'([A-Z][a-z]*)(\d*)', formula_str)
        return {element: int(count) if count else 1 for element, count in parts}

    def check_dou(formula_parts):
        # DoU = C + 1 - H/2 + N/2
        c = formula_parts.get('C', 0)
        h = formula_parts.get('H', 0)
        n = formula_parts.get('N', 0)
        return c + 1 - (h / 2) + (n / 2)

    # --- Verification Logic ---
    correct_candidate_key = None
    all_reasons = {}

    for key, candidate in candidates.items():
        reasons_for_failure = []
        
        # 1. Check Molecular Formula and DoU
        candidate_formula_parts = get_formula_parts(candidate["formula"])
        question_formula_parts = get_formula_parts(question_data["molecular_formula"])
        if candidate_formula_parts != question_formula_parts:
            reasons_for_failure.append(f"Incorrect molecular formula. Expected {question_data['molecular_formula']}, but candidate is {candidate['formula']}.")
        
        if check_dou(candidate_formula_parts) != 5:
             reasons_for_failure.append(f"Incorrect Degree of Unsaturation. Expected 5, but candidate has {check_dou(candidate_formula_parts)}.")

        # 2. Check IR Data
        features = candidate["features"]
        if features["amine_type"] == "primary_amine" or features["amine_type"] == "primary_amide":
            if len(question_data["ir_peaks"]["N-H"]) != 2:
                reasons_for_failure.append("IR mismatch: Expected 2 N-H peaks for a primary amine/amide, but data shows otherwise.")
        elif features["amine_type"] == "secondary_amide":
            if len(question_data["ir_peaks"]["N-H"]) != 1:
                 reasons_for_failure.append("IR mismatch: Data shows 2 N-H peaks, which is inconsistent with a secondary amide.")
        
        if features["carbonyl_type"] == "ester":
            # Conjugated ester C=O is ~1715-1730 cm-1. 1720 is a perfect match.
            if not (1715 <= question_data["ir_peaks"]["C=O"] <= 1730):
                reasons_for_failure.append(f"IR mismatch: C=O peak at {question_data['ir_peaks']['C=O']} is not typical for a conjugated ester.")
        elif features["carbonyl_type"] == "amide":
            # Conjugated amide C=O is ~1650-1690 cm-1. 1720 is too high.
            if not (1650 <= question_data["ir_peaks"]["C=O"] <= 1690):
                reasons_for_failure.append(f"IR mismatch: C=O peak at {question_data['ir_peaks']['C=O']} is too high for an amide.")

        # 3. Check NMR Data
        # Check aromatic substitution
        aromatic_doublets = [s for s in question_data["nmr_signals"] if s["type"] == "d" and s["H"] == 2]
        if features["substitution"] == "para":
            if len(aromatic_doublets) != 2:
                reasons_for_failure.append("NMR mismatch: Data does not show the two doublets expected for a para-substituted ring.")
        elif features["substitution"] == "meta":
            if len(aromatic_doublets) == 2:
                reasons_for_failure.append("NMR mismatch: Data shows two clean doublets, which is inconsistent with meta-substitution.")

        # Check ethyl group quartet chemical shift
        quartet_signal = next((s for s in question_data["nmr_signals"] if s["type"] == "q"), None)
        if quartet_signal:
            quartet_ppm = quartet_signal["ppm"]
            if features["ethyl_type"] == "ethyl_ester": # -O-CH2-
                if not (4.2 <= quartet_ppm <= 4.6):
                    reasons_for_failure.append(f"NMR mismatch: Quartet at {quartet_ppm} ppm is outside the expected range for an ethyl ester (-O-CH2-).")
            elif features["ethyl_type"] == "propionate": # -CO-CH2-
                if not (2.2 <= quartet_ppm <= 2.8):
                    reasons_for_failure.append(f"NMR mismatch: Quartet at {quartet_ppm} ppm is incorrect for a propionate (-CO-CH2-), which should be ~2.5 ppm.")
            elif features["ethyl_type"] == "ethoxy": # Ar-O-CH2-
                 if not (3.9 <= quartet_ppm <= 4.2):
                    reasons_for_failure.append(f"NMR mismatch: Quartet at {quartet_ppm} ppm is slightly high for a typical ethoxy group.")

        if not reasons_for_failure:
            correct_candidate_key = key
        else:
            all_reasons[key] = reasons_for_failure

    # --- Final Verdict ---
    if correct_candidate_key == provided_answer:
        return "Correct"
    elif correct_candidate_key is None:
        return "The provided answer is incorrect. No candidate perfectly matches the data based on the checks."
    else:
        correct_name = candidates[correct_candidate_key]['name']
        provided_name = candidates[provided_answer]['name']
        failure_reasons = all_reasons.get(provided_answer, ["Unknown reason."])
        return (f"Incorrect. The provided answer is {provided_answer} ({provided_name}), but the data points to {correct_candidate_key} ({correct_name}).\n"
                f"The provided answer {provided_answer} is wrong because:\n- " + "\n- ".join(failure_reasons))

# Run the check
result = check_spectroscopy_answer()
print(result)