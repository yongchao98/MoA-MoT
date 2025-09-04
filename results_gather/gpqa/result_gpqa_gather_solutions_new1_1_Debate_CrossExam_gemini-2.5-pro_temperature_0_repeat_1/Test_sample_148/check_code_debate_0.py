def check_answer_correctness():
    """
    Checks the correctness of the final answer by logically evaluating each option
    against the experimental constraints provided in the question.
    """
    # The final answer from the LLM analysis to be checked.
    final_answer = "D"

    # --- Define Experimental Observations as Constraints ---
    # 1. From LC-MS: "Both peaks have the same mass spectrum"
    # This means the two species are isomers.
    constraint_are_isomers = True
    
    # 2. From LC-MS: "two clearly defined peaks"
    # This means the two species are separable by standard (achiral) chromatography.
    constraint_are_separable_by_lc = True
    
    # 3. From 1H NMR: "two peaks that both correspond to the same alpha-proton"
    # This means the two species are distinguishable by standard (achiral) NMR.
    constraint_are_distinguishable_by_nmr = True

    # --- Define Properties of Each Option ---
    options = {
        "A": {
            "description": "Mixture of enantiomers",
            "is_isomer": True,
            "separable_by_achiral_lc": False,
            "distinguishable_by_achiral_nmr": False
        },
        "B": {
            "description": "Contaminated with a precursor",
            "is_isomer": False,  # A precursor has a different (lower) mass.
            "separable_by_achiral_lc": True,
            "distinguishable_by_achiral_nmr": True
        },
        "C": {
            "description": "'Double coupling' has occurred",
            "is_isomer": False,  # This product would have a different (higher) mass.
            "separable_by_achiral_lc": True,
            "distinguishable_by_achiral_nmr": True
        },
        "D": {
            "description": "Mixture of diastereoisomers",
            "is_isomer": True,
            "separable_by_achiral_lc": True,
            "distinguishable_by_achiral_nmr": True
        }
    }

    # --- Apply Constraints to Find the Correct Answer ---
    valid_options = []
    for option_key, properties in options.items():
        # Check against Constraint 1: Must be isomers (same mass)
        if properties["is_isomer"] != constraint_are_isomers:
            continue
        
        # Check against Constraint 2: Must be separable by LC
        if properties["separable_by_achiral_lc"] != constraint_are_separable_by_lc:
            continue

        # Check against Constraint 3: Must be distinguishable by NMR
        if properties["distinguishable_by_achiral_nmr"] != constraint_are_distinguishable_by_nmr:
            continue
            
        # If all constraints are met, this is a valid option
        valid_options.append(option_key)

    # --- Final Verification ---
    if len(valid_options) == 1:
        derived_correct_answer = valid_options[0]
        if derived_correct_answer == final_answer:
            return "Correct"
        else:
            return (f"Incorrect. The provided answer is '{final_answer}', but the only option that satisfies all "
                    f"experimental constraints is '{derived_correct_answer}' ({options[derived_correct_answer]['description']}).")
    elif len(valid_options) == 0:
        return "Incorrect. No option satisfies all the experimental constraints."
    else:
        return (f"Incorrect. The problem is ambiguous as multiple options ({', '.join(valid_options)}) "
                f"satisfy all experimental constraints.")

# Run the check
result = check_answer_correctness()
print(result)