def check_answer():
    """
    Checks the correctness of the answer by analyzing the spectral data for C9H11NO2.
    """
    # Define the properties of each candidate molecule
    candidates = {
        'A': {
            "name": "ethyl 4-aminobenzoate",
            "amine_type": "primary",  # -NH2
            "carbonyl_type": "ester", # -COOR
            "substitution": "para",   # 1,4-disubstituted
            "ethyl_attachment": "ester_oxygen" # -COOCH2CH3
        },
        'B': {
            "name": "N-(4-ethoxyphenyl)formamide",
            "amine_type": "secondary_amide", # -NH-
            "carbonyl_type": "amide", # -CONH-
            "substitution": "para",
            "ethyl_attachment": "ether_oxygen" # -OCH2CH3
        },
        'C': {
            "name": "3-ethoxybenzamide",
            "amine_type": "primary_amide", # -CONH2
            "carbonyl_type": "amide",
            "substitution": "meta", # 1,3-disubstituted
            "ethyl_attachment": "ether_oxygen"
        },
        'D': {
            "name": "4-aminophenyl propionate",
            "amine_type": "primary",
            "carbonyl_type": "ester",
            "substitution": "para",
            "ethyl_attachment": "carbonyl_carbon" # -COCH2CH3
        }
    }

    # The final answer provided by the LLM
    llm_answer_key = 'A'
    
    # Select the candidate structure based on the LLM's answer
    chosen_candidate = candidates.get(llm_answer_key)

    if not chosen_candidate:
        return f"Invalid answer key '{llm_answer_key}'. The key must be one of {list(candidates.keys())}."

    # --- Start of Checks based on Spectral Data ---

    # Check 1: IR bands at 3420 & 3325 cm-1 indicate a primary amine or primary amide.
    # This means it must have an -NH2 group.
    if chosen_candidate["amine_type"] not in ["primary", "primary_amide"]:
        return (f"Incorrect. The IR spectrum shows two bands at 3420 and 3325 cm-1, "
                f"which is characteristic of a primary amine or primary amide (-NH2). "
                f"Candidate {llm_answer_key} ({chosen_candidate['name']}) is a {chosen_candidate['amine_type']}, which would show only one N-H band.")

    # Check 2: Strong IR band at 1720 cm-1 indicates a conjugated ester.
    # Amides typically absorb at a lower frequency (< 1700 cm-1).
    if chosen_candidate["carbonyl_type"] != "ester":
        return (f"Incorrect. The IR spectrum shows a strong band at 1720 cm-1, "
                f"which is characteristic of a conjugated ester. "
                f"Candidate {llm_answer_key} ({chosen_candidate['name']}) is an {chosen_candidate['carbonyl_type']}, "
                f"which would absorb at a lower frequency (typically < 1700 cm-1).")

    # Check 3: 1H NMR shows two doublets in the aromatic region (7.0 & 8.0 ppm).
    # This is the classic pattern for a 1,4- (para) disubstituted ring.
    if chosen_candidate["substitution"] != "para":
        return (f"Incorrect. The 1H NMR spectrum shows two doublets in the aromatic region, "
                f"which indicates a 1,4- (para) disubstituted ring. "
                f"Candidate {llm_answer_key} ({chosen_candidate['name']}) is {chosen_candidate['substitution']}-substituted.")

    # Check 4: 1H NMR shows a quartet at 4.5 ppm.
    # This indicates a -CH2- group attached to an oxygen, not a carbonyl.
    # This distinguishes an ethyl ester (-COOCH2CH3) from a propionate (-COCH2CH3).
    if chosen_candidate["ethyl_attachment"] not in ["ester_oxygen", "ether_oxygen"]:
        return (f"Incorrect. The 1H NMR spectrum shows a quartet at a high chemical shift of 4.5 ppm, "
                f"indicating a -CH2- group is attached to an oxygen. "
                f"In candidate {llm_answer_key} ({chosen_candidate['name']}), the ethyl group's -CH2- is attached to a "
                f"{chosen_candidate['ethyl_attachment']}, which would result in a quartet at a much lower chemical shift (~2.5 ppm).")

    # Check 5: 1H NMR shows a broad singlet for 2H at 4.0 ppm.
    # This confirms the presence of a primary amine (-NH2).
    if chosen_candidate["amine_type"] != "primary":
        # This check is slightly redundant with check 1 & 2 but confirms the amine vs amide distinction.
        return (f"Incorrect. The 1H NMR shows a broad singlet for 2H at 4.0 ppm, confirming a primary amine (-NH2). "
                f"Candidate {llm_answer_key} ({chosen_candidate['name']}) is an {chosen_candidate['amine_type']}, which is inconsistent.")

    # If all checks pass for the chosen candidate
    return "Correct"

# Run the check
result = check_answer()
print(result)