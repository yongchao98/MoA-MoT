def check_answer():
    """
    Checks the correctness of the identified compound based on spectral data.
    """
    # --- Define Constraints from Spectral Data ---
    constraints = {
        "formula": "C9H11NO2",
        "amine_type": "primary",  # From IR (2 bands) and NMR (bs, 2H)
        "carbonyl_type": "ester", # From IR (1720 cm-1)
        "nmr_ethyl_group_type": "ester_O_CH2CH3", # From NMR (quartet at 4.5 ppm)
        "substitution_pattern": "para" # From NMR (2 doublets, 2H each)
    }

    # --- Define Properties of Candidate Molecules ---
    candidates = {
        "A": {
            "name": "ethyl 4-aminobenzoate",
            "formula": "C9H11NO2",
            "amine_type": "primary",
            "carbonyl_type": "ester",
            "nmr_ethyl_group_type": "ester_O_CH2CH3", # -O-CH2CH3 shift is high (~4.3 ppm)
            "substitution_pattern": "para"
        },
        "B": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "formula": "C9H11NO2",
            "amine_type": "secondary_amide", # It's a secondary amide (-NH-), not a primary amine.
            "carbonyl_type": "amide",
            "nmr_ethyl_group_type": "ether_O_CH2CH3", # -O-CH2CH3 shift is ~4.0 ppm
            "substitution_pattern": "para"
        },
        "C": {
            "name": "4-aminophenyl propionate",
            "formula": "C9H11NO2",
            "amine_type": "primary",
            "carbonyl_type": "ester",
            "nmr_ethyl_group_type": "ester_CO_CH2CH3", # -CO-CH2CH3 shift is low (~2.5 ppm)
            "substitution_pattern": "para"
        },
        "D": {
            "name": "3-ethoxybenzamide",
            "formula": "C9H11NO2",
            "amine_type": "primary_amide", # It's a primary amide (-CONH2), not an amine on the ring.
            "carbonyl_type": "amide",
            "nmr_ethyl_group_type": "ether_O_CH2CH3", # -O-CH2CH3 shift is ~4.0 ppm
            "substitution_pattern": "meta" # 1,3-disubstituted
        }
    }

    llm_answer = "A" # The final answer provided by the LLM.

    correct_candidate = None
    errors = {}

    for key, properties in candidates.items():
        is_match = True
        error_log = []

        # Check formula
        if properties["formula"] != constraints["formula"]:
            is_match = False
            error_log.append(f"Formula mismatch: Expected {constraints['formula']}, got {properties['formula']}")

        # Check amine type
        if properties["amine_type"] != constraints["amine_type"]:
            is_match = False
            error_log.append(f"Amine type mismatch: Expected '{constraints['amine_type']}', but structure is '{properties['amine_type']}' (IR shows two N-H bands).")

        # Check carbonyl type
        if properties["carbonyl_type"] != constraints["carbonyl_type"]:
            is_match = False
            error_log.append(f"Carbonyl type mismatch: Expected '{constraints['carbonyl_type']}', but structure is '{properties['carbonyl_type']}' (IR shows C=O at 1720 cm-1).")
        
        # Check NMR ethyl group
        if properties["nmr_ethyl_group_type"] != constraints["nmr_ethyl_group_type"]:
            is_match = False
            error_log.append(f"NMR ethyl group mismatch: Expected '{constraints['nmr_ethyl_group_type']}', but structure is '{properties['nmr_ethyl_group_type']}' (NMR quartet is at a high shift of 4.5 ppm).")

        # Check substitution pattern
        if properties["substitution_pattern"] != constraints["substitution_pattern"]:
            is_match = False
            error_log.append(f"Substitution pattern mismatch: Expected '{constraints['substitution_pattern']}', but structure is '{properties['substitution_pattern']}' (NMR shows two doublets in aromatic region).")

        if is_match:
            if correct_candidate is None:
                correct_candidate = key
            else:
                # This case would mean the question is ambiguous
                return "Error: Multiple candidates match the spectral data."
        else:
            errors[key] = error_log

    if correct_candidate is None:
        return "Error: No candidate matches all the spectral data."

    if llm_answer == correct_candidate:
        return "Correct"
    else:
        reason = f"Incorrect. The correct answer is {correct_candidate} ({candidates[correct_candidate]['name']}).\n"
        reason += f"The provided answer was {llm_answer}, which is incorrect because:\n"
        for err in errors[llm_answer]:
            reason += f"- {err}\n"
        return reason

# Run the check
result = check_answer()
print(result)