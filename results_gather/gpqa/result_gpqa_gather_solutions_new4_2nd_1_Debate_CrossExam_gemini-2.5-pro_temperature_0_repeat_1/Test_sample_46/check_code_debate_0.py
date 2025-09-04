import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer for a spectroscopy problem by systematically
    evaluating each possible compound against the provided spectral data.
    """
    # --- Data from the Question ---
    question_data = {
        "formula": "C9H11NO2",
        "ir": {
            "N-H_stretches": 2,  # Two peaks at 3420, 3325 cm-1 indicate primary amine/amide
            "C=O_stretch": 1720  # cm-1
        },
        "nmr": {
            "ethyl_quartet_ppm": 4.5,
            "aromatic_pattern": "two_doublets", # Two doublets (2H each) indicate para-substitution
            "amine_protons": True # bs, 2H at 4.0 ppm
        }
    }

    # --- Candidate Compounds and their expected spectral features ---
    # This dictionary defines the expected properties for each option based on its structure.
    options = {
        "A": {
            "name": "ethyl 4-aminobenzoate",
            "formula": "C9H11NO2",
            "N-H_type": "primary_amine",    # Expects 2 N-H IR peaks
            "C=O_type": "conjugated_ester", # Expects C=O ~1715-1730 cm-1
            "aromatic_type": "para",        # Expects two doublets in NMR
            "ethyl_quartet_env": "O-CH2"    # Expects quartet ~4.1-4.6 ppm
        },
        "B": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "formula": "C9H11NO2",
            "N-H_type": "secondary_amide",  # Expects 1 N-H IR peak
            "C=O_type": "amide",            # Expects C=O ~1650-1690 cm-1
            "aromatic_type": "para",
            "ethyl_quartet_env": "ArO-CH2"   # Expects quartet ~3.9-4.1 ppm
        },
        "C": {
            "name": "3-ethoxybenzamide",
            "formula": "C9H11NO2",
            "N-H_type": "primary_amide",    # Expects 2 N-H IR peaks
            "C=O_type": "amide",            # Expects C=O ~1650-1690 cm-1
            "aromatic_type": "meta",        # Expects complex NMR pattern, not two doublets
            "ethyl_quartet_env": "ArO-CH2"   # Expects quartet ~3.9-4.1 ppm
        },
        "D": {
            "name": "4-aminophenyl propionate",
            "formula": "C9H11NO2",
            "N-H_type": "primary_amine",    # Expects 2 N-H IR peaks
            "C=O_type": "ester",            # Expects C=O ~1730-1750 cm-1
            "aromatic_type": "para",
            "ethyl_quartet_env": "CO-CH2"    # Expects quartet ~2.2-2.6 ppm
        }
    }

    # The final answer provided by the LLM
    llm_answer_choice = "A"

    # --- Analysis ---
    # Find the candidate that matches all spectral data
    identified_compound_key = None
    error_log = {}

    for key, props in options.items():
        errors = []
        
        # Check 1: IR N-H stretch
        if props["N-H_type"] in ["primary_amine", "primary_amide"]:
            if question_data["ir"]["N-H_stretches"] != 2:
                errors.append("IR mismatch: Expected 2 N-H peaks for primary amine/amide.")
        elif props["N-H_type"] == "secondary_amide":
            if question_data["ir"]["N-H_stretches"] != 1:
                errors.append("IR mismatch: Expected 1 N-H peak for secondary amide, but data shows 2.")

        # Check 2: IR C=O stretch
        co_stretch = question_data["ir"]["C=O_stretch"]
        if "ester" in props["C=O_type"]:
            if not (1715 <= co_stretch <= 1750):
                errors.append(f"IR mismatch: C=O at {co_stretch} is outside the expected ester range.")
        elif props["C=O_type"] == "amide":
            if not (1650 <= co_stretch <= 1690):
                errors.append(f"IR mismatch: C=O at {co_stretch} is too high for an amide (expected < 1700).")

        # Check 3: NMR Aromatic Pattern
        if props["aromatic_type"] == "para":
            if question_data["nmr"]["aromatic_pattern"] != "two_doublets":
                errors.append("NMR mismatch: Expected para-pattern (two doublets).")
        elif props["aromatic_type"] == "meta":
            if question_data["nmr"]["aromatic_pattern"] == "two_doublets":
                errors.append("NMR mismatch: Expected complex meta-pattern, not simple para-pattern.")

        # Check 4: NMR Ethyl Quartet Shift
        q_ppm = question_data["nmr"]["ethyl_quartet_ppm"]
        env = props["ethyl_quartet_env"]
        if env == "O-CH2": # e.g., ethyl ester
            if not (4.1 <= q_ppm <= 4.6):
                errors.append(f"NMR mismatch: Quartet at {q_ppm} ppm is outside the 4.1-4.6 ppm range for an -O-CH2- group.")
        elif env == "CO-CH2": # e.g., propionate
            if not (2.2 <= q_ppm <= 2.6):
                errors.append(f"NMR mismatch: Quartet at {q_ppm} ppm is not in the 2.2-2.6 ppm range for a -CO-CH2- group.")
        elif env == "ArO-CH2": # e.g., ethoxy
             if not (3.9 <= q_ppm <= 4.1):
                errors.append(f"NMR mismatch: Quartet at {q_ppm} ppm is not in the 3.9-4.1 ppm range for an Ar-O-CH2- group.")

        if not errors:
            identified_compound_key = key
        else:
            error_log[key] = errors

    # --- Verdict ---
    if identified_compound_key == llm_answer_choice:
        return "Correct"
    elif identified_compound_key is None:
        return "Analysis Error: No single option correctly fits all the provided spectral data."
    else:
        correct_name = options[identified_compound_key]['name']
        llm_name = options[llm_answer_choice]['name']
        reason = f"Incorrect. The provided answer is {llm_answer_choice} ({llm_name}), but the spectral data points to {identified_compound_key} ({correct_name}).\n"
        
        # Provide specific reasons why the LLM's choice is wrong
        if llm_answer_choice in error_log:
            reason += f"The chosen answer {llm_answer_choice} is incorrect because:\n"
            for err in error_log[llm_answer_choice]:
                reason += f"- {err}\n"
        else:
             reason += "An unknown error occurred during analysis."
        return reason.strip()

# Execute the check and print the result
result = check_correctness()
print(result)