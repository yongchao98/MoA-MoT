import re

def check_peptide_analysis_answer():
    """
    Checks the correctness of the answer to the peptide analysis question.

    This function simulates the logical deduction process based on the provided
    analytical data (NMR, LC-MS) to determine the correct explanation among the
    given options. It then compares this derived correct answer to the provided
    LLM's answer.
    """
    # The final answer from the LLM to be checked.
    llm_final_answer = "<<<D>>>"

    # --- Step 1: Define the constraints from the problem description ---
    # Constraint 1: Both species have the same mass (from MS data).
    # Constraint 2: The two species are separable by standard LC (achiral).
    # Constraint 3: The two species are distinguishable by standard NMR (achiral).
    observations = {
        "is_isomer": True,
        "is_separable_by_achiral_lc": True,
        "is_distinguishable_by_achiral_nmr": True
    }

    # --- Step 2: Define the properties of each possible answer choice ---
    # Based on the original question's lettering.
    options_properties = {
        "A": {  # Precursor
            "name": "Contaminated with a precursor",
            "is_isomer": False,  # Different molecule, different mass
            "is_separable_by_achiral_lc": True,
            "is_distinguishable_by_achiral_nmr": True
        },
        "B": {  # Enantiomers
            "name": "Mixture of enantiomers",
            "is_isomer": True,
            "is_separable_by_achiral_lc": False,  # Not separable by achiral LC
            "is_distinguishable_by_achiral_nmr": False  # Not distinguishable by achiral NMR
        },
        "C": {  # Double coupling product
            "name": "'Double coupling' has occurred",
            "is_isomer": False,  # Different molecule, different mass
            "is_separable_by_achiral_lc": True,
            "is_distinguishable_by_achiral_nmr": True
        },
        "D": {  # Diastereoisomers
            "name": "Mixture of diastereoisomers",
            "is_isomer": True,
            "is_separable_by_achiral_lc": True,
            "is_distinguishable_by_achiral_nmr": True
        }
    }

    # --- Step 3: Find the correct option by matching properties to observations ---
    correct_option_key = None
    for option_key, properties in options_properties.items():
        # Check if all observations match the properties of this option
        if all(properties.get(obs_key) == obs_value for obs_key, obs_value in observations.items()):
            correct_option_key = option_key
            break

    # --- Step 4: Validate the LLM's answer ---
    # Extract the letter from the LLM's answer string
    match = re.search(r'<<<([A-D])>>>', llm_final_answer)
    if not match:
        return f"Invalid answer format: '{llm_final_answer}'. Expected format is '<<<X>>>' where X is A, B, C, or D."

    llm_answer_key = match.group(1)

    if llm_answer_key == correct_option_key:
        return "Correct"
    else:
        # Build a detailed explanation for the error
        reason = f"The provided answer '{llm_answer_key}' is incorrect. The correct answer is '{correct_option_key}'.\n\n"
        reason += "Reasoning:\n"
        reason += "1. The MS data shows both peaks have the same mass, meaning they are isomers. This rules out options A (precursor) and C (double coupling product), which would have different masses.\n"
        reason += "2. The LC and NMR data show two distinct, separable peaks using standard (achiral) methods. This rules out option B (enantiomers), as they are indistinguishable by achiral techniques.\n"
        reason += f"3. Only option D ({options_properties['D']['name']}) is consistent with all observations: diastereomers are isomers (same mass) but have different physical/chemical properties, making them separable by LC and distinguishable by NMR."
        return reason

# Execute the check and print the result
result = check_peptide_analysis_answer()
print(result)