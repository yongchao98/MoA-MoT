import sys
from io import StringIO

def check_compound_identification():
    """
    This function checks the correctness of the identified compound based on spectral data.
    It simulates the step-by-step analysis a chemist would perform.
    """
    
    # --- Define Spectral Data Constraints from the Question ---
    constraints = {
        "molecular_formula": "C9H11NO2",
        "ir_amine": "primary",  # Two bands at 3420, 3325 cm-1 indicate -NH2
        "ir_carbonyl": "ester", # Strong band at 1720 cm-1 indicates ester C=O
        "nmr_aromatic_pattern": "para", # Two doublets (2H each) indicate 1,4-disubstitution
        "nmr_ethyl_group_type": "ethyl_ester" # Quartet at high shift (4.5 ppm) indicates -O-CH2CH3
    }

    # --- Define Properties of Candidate Compounds ---
    candidates = {
        "A": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "formula": "C9H11NO2",
            "amine_type": "secondary_amide", # -NH-CHO group
            "carbonyl_type": "amide",
            "aromatic_pattern": "para",
            "ethyl_group_type": "ethoxy" # Ar-O-CH2CH3
        },
        "B": {
            "name": "ethyl 4-aminobenzoate",
            "formula": "C9H11NO2",
            "amine_type": "primary", # -NH2 group
            "carbonyl_type": "ester",
            "aromatic_pattern": "para",
            "ethyl_group_type": "ethyl_ester" # Ar-COO-CH2CH3
        },
        "C": {
            "name": "4-aminophenyl propionate",
            "formula": "C9H11NO2",
            "amine_type": "primary", # -NH2 group
            "carbonyl_type": "ester",
            "aromatic_pattern": "para",
            "ethyl_group_type": "propionate" # Ar-O-CO-CH2CH3
        },
        "D": {
            "name": "3-ethoxybenzamide",
            "formula": "C9H11NO2",
            "amine_type": "primary_amide", # -CONH2 group
            "carbonyl_type": "amide",
            "aromatic_pattern": "meta",
            "ethyl_group_type": "ethoxy" # Ar-O-CH2CH3
        }
    }

    # --- The Answer to Check ---
    llm_answer = "B"

    # --- Analysis Logic ---
    correct_candidates = []
    reasons = {}

    for key, props in candidates.items():
        is_correct = True
        reason_list = []

        # Constraint 1: IR Amine Type (Primary Amine)
        if props["amine_type"] != constraints["ir_amine"]:
            is_correct = False
            reason_list.append(f"Fails IR amine check: Has a {props['amine_type']}, but data (3420, 3325 cm-1) indicates a primary amine.")
        
        # Constraint 2: IR Carbonyl Type (Ester)
        if props["carbonyl_type"] != constraints["ir_carbonyl"]:
            is_correct = False
            reason_list.append(f"Fails IR carbonyl check: Has an {props['carbonyl_type']}, but data (1720 cm-1) indicates an ester.")

        # Constraint 3: NMR Aromatic Pattern (Para)
        if props["aromatic_pattern"] != constraints["nmr_aromatic_pattern"]:
            is_correct = False
            reason_list.append(f"Fails NMR aromatic check: Is {props['aromatic_pattern']}-substituted, but data (two doublets) indicates para-substitution.")

        # Constraint 4: NMR Ethyl Group Type (Ethyl Ester)
        # This is the key differentiator between B and C.
        # An ethyl ester (-COO-CH2CH3) has the -CH2- quartet at ~4.1-4.6 ppm.
        # A propionate (-CO-CH2CH3) has the -CH2- quartet at ~2.2-2.6 ppm.
        if props["ethyl_group_type"] != constraints["nmr_ethyl_group_type"]:
            is_correct = False
            reason_list.append(f"Fails NMR ethyl group check: Has a {props['ethyl_group_type']} group. The data (quartet at 4.5 ppm) specifically indicates an ethyl ester (-O-CH2-), not a propionate (-CO-CH2-) or simple ethoxy group.")

        if is_correct:
            correct_candidates.append(key)
        
        reasons[key] = reason_list

    # --- Final Verdict ---
    if not correct_candidates:
        return "Incorrect. No candidate satisfies all the spectral data constraints."
    
    if len(correct_candidates) > 1:
        return f"Incorrect. Multiple candidates ({', '.join(correct_candidates)}) satisfy the constraints. The data may be ambiguous."

    determined_correct_answer = correct_candidates[0]

    if llm_answer == determined_correct_answer:
        return "Correct"
    else:
        error_message = f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {determined_correct_answer}.\n"
        error_message += f"Reasoning:\n"
        error_message += f"The provided answer {llm_answer} ({candidates[llm_answer]['name']}) is incorrect because:\n"
        for reason in reasons[llm_answer]:
            error_message += f"- {reason}\n"
        
        error_message += f"\nCandidate {determined_correct_answer} ({candidates[determined_correct_answer]['name']}) is the only one that satisfies all constraints:\n"
        error_message += f"- It has a primary amine (-NH2), matching the IR bands at 3420 and 3325 cm-1.\n"
        error_message += f"- It has a conjugated ester, matching the IR band at 1720 cm-1.\n"
        error_message += f"- It is para-substituted, matching the two doublets in the NMR aromatic region.\n"
        error_message += f"- It has an ethyl ester group (-COO-CH2CH3), matching the NMR triplet at 1.20 ppm and the highly deshielded quartet at 4.5 ppm."
        
        return error_message

# Execute the check and print the result
result = check_compound_identification()
print(result)