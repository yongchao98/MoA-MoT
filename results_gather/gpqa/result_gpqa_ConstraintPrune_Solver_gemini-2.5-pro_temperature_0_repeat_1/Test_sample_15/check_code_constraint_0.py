def check_answer():
    """
    This function checks the correctness of the provided answer about the number of optically active compounds.
    It uses a ground truth dictionary where each compound's optical activity is determined by established chemical principles.
    """

    # Ground truth: A dictionary mapping compound names to their optical activity status (True for active, False for inactive).
    ground_truth_activity = {
        "(Z)-1-chloro-2-methylbut-1-ene": False,  # Achiral (plane of symmetry, no chiral centers)
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione": True,  # Chiral (specific, non-meso stereoisomer)
        "(2R,3S)-2,3-dimethylsuccinic acid": False, # Achiral (meso compound)
        "(2R,3R)-2,3-dimethylsuccinic acid": True,   # Chiral (specific enantiomer)
        "(R)-cyclohex-3-en-1-ol": True,   # Chiral (one chiral center, specific enantiomer)
        "(1s,3s,5s)-cyclohexane-1,3,5-triol": False, # Achiral (high symmetry)
        "1-cyclopentyl-3-methylbutan-1-one": False  # Achiral (no chiral centers)
    }

    # The answer provided by the LLM to be checked.
    llm_answer_count = 3
    llm_answer_option = "D"
    
    # The options provided in the question.
    options = {'A': 5, 'B': 2, 'C': 4, 'D': 3}

    # 1. Check if the LLM's option letter corresponds to its stated count.
    if options.get(llm_answer_option) != llm_answer_count:
        return f"Incorrect. The LLM's final option '{llm_answer_option}' corresponds to a count of {options.get(llm_answer_option)}, but its analysis concluded a count of {llm_answer_count}."

    # 2. Calculate the correct count based on the ground truth.
    correct_count = sum(1 for is_active in ground_truth_activity.values() if is_active)

    # 3. Compare the LLM's count with the correct count.
    if llm_answer_count != correct_count:
        active_compounds = [name for name, is_active in ground_truth_activity.items() if is_active]
        return (f"Incorrect. The LLM's answer is {llm_answer_count}, but the correct count of optically active compounds is {correct_count}.\n"
                f"The optically active compounds are:\n" + "\n".join(f"- {c}" for c in active_compounds))

    # 4. Check if the LLM correctly identified the active compounds.
    llm_identified_active = [
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
        "(2R,3R)-2,3-dimethylsuccinic acid",
        "(R)-cyclohex-3-en-1-ol"
    ]
    ground_truth_active = [name for name, is_active in ground_truth_activity.items() if is_active]

    if sorted(llm_identified_active) != sorted(ground_truth_active):
        return (f"Incorrect. The LLM arrived at the correct count ({correct_count}) but for the wrong reasons. "
                f"The list of active compounds is not correct.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)