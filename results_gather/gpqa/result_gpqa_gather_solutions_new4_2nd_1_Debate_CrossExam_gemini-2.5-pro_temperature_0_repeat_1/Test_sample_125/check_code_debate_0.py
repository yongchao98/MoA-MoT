def check_organic_chemistry_answer():
    """
    This function checks the correctness of the answer to the organic chemistry HPLC problem.
    It simulates the chemical analysis step-by-step to determine the correct number of peaks
    and compares it with the provided answer.
    """

    # --- Step 1: Define the problem's options and the LLM's chosen answer ---
    options = {
        "A": {"chiral": 5, "normal": 4},
        "B": {"chiral": 3, "normal": 2},
        "C": {"chiral": 3, "normal": 3},
        "D": {"chiral": 4, "normal": 2}
    }
    llm_answer_choice = "A"

    # --- Step 2: Perform the chemical analysis programmatically ---

    # Reaction I: (S)-5-methoxyhexan-3-one (chiral) is reduced, creating one new stereocenter.
    # This results in two products that are diastereomers of each other.
    # Diastereomers are separable by both chiral and normal-phase HPLC.
    rxn1_products = {
        "description": "Two diastereomers",
        "total_stereoisomers": 2,
        "achiral_separable_groups": 2  # Each diastereomer is a separate group
    }

    # Reaction II: Pentane-2,4-dione (achiral) is reduced, creating two new stereocenters.
    # This results in a pair of enantiomers and a meso compound.
    # The meso compound is a diastereomer of the enantiomeric pair.
    rxn2_products = {
        "description": "One enantiomeric pair and one meso compound",
        "total_stereoisomers": 3,  # (R,R), (S,S), and meso
        "achiral_separable_groups": 2  # The enantiomeric pair (1 group) + the meso compound (1 group)
    }

    # --- Step 3: Calculate the expected number of peaks for each HPLC type ---

    # Chiral HPLC can separate all unique stereoisomers.
    # We assume products from Rxn I and Rxn II are structurally different and will not co-elute.
    expected_chiral_peaks = rxn1_products["total_stereoisomers"] + rxn2_products["total_stereoisomers"]

    # Normal-phase (achiral) HPLC separates diastereomers but not enantiomers.
    # Enantiomeric pairs elute as a single peak.
    expected_normal_phase_peaks = rxn1_products["achiral_separable_groups"] + rxn2_products["achiral_separable_groups"]

    # --- Step 4: Validate the LLM's answer ---
    if llm_answer_choice not in options:
        return f"The answer choice '{llm_answer_choice}' is not a valid option (A, B, C, or D)."

    llm_answer_values = options[llm_answer_choice]

    # Check if the values in the chosen option match the calculated correct values.
    is_chiral_correct = llm_answer_values["chiral"] == expected_chiral_peaks
    is_normal_correct = llm_answer_values["normal"] == expected_normal_phase_peaks

    if is_chiral_correct and is_normal_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_chiral_correct:
            error_messages.append(
                f"the number of peaks in chiral HPLC is incorrect. Expected {expected_chiral_peaks}, but the answer claims {llm_answer_values['chiral']}."
            )
        if not is_normal_correct:
            error_messages.append(
                f"the number of peaks in normal-phase HPLC is incorrect. Expected {expected_normal_phase_peaks}, but the answer claims {llm_answer_values['normal']}."
            )
        return "The answer is incorrect because " + " and ".join(error_messages)

# Execute the check and print the result
result = check_organic_chemistry_answer()
print(result)