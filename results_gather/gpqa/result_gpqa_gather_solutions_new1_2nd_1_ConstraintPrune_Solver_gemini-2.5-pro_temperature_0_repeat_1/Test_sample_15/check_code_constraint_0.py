def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing each compound for optical activity.
    """
    # The final answer provided by the LLM to be checked.
    llm_final_answer = "<<<A>>>"

    # Define the mapping from options to the number of compounds.
    # This is based on the question text: A) 3, B) 2, C) 5, D) 4
    options_map = {
        "A": 3,
        "B": 2,
        "C": 5,
        "D": 4
    }

    # A dictionary to store the analysis of each compound.
    # The value is a tuple: (is_optically_active, reason).
    compounds_analysis = {
        "(Z)-1-chloro-2-methylbut-1-ene": 
            (False, "Achiral. The molecule is planar around the C=C bond, creating a plane of symmetry."),
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione": 
            (True, "Chiral. The name specifies a single enantiomer (3aR,7aS) of a complex chiral molecule."),
        "(2R,3S)-2,3-dimethylsuccinic acid": 
            (False, "Achiral. This is a meso compound with an internal plane of symmetry, making it optically inactive despite having two chiral centers."),
        "(2R,3R)-2,3-dimethylsuccinic acid": 
            (True, "Chiral. This is a single enantiomer (2R,3R) of a chiral pair, lacking internal symmetry."),
        "(R)-cyclohex-3-en-1-ol": 
            (True, "Chiral. The name specifies a single enantiomer (R) with a chiral center at C1."),
        "(1s,3s,5s)-cyclohexane-1,3,5-triol": 
            (False, "Achiral. The all-cis isomer is highly symmetric and possesses multiple planes of symmetry."),
        "1-cyclopentyl-3-methylbutan-1-one": 
            (False, "Achiral. The molecule lacks a chiral center (C3 of the butane chain has two identical methyl groups).")
    }

    # Calculate the correct number of optically active compounds from our analysis.
    correct_count = sum(1 for is_active, reason in compounds_analysis.values() if is_active)
    
    # Extract the chosen option from the LLM's answer.
    try:
        chosen_option = llm_final_answer.strip().replace('<', '').replace('>', '')
    except Exception:
        return "Error: The provided answer format is invalid."

    # Check if the chosen option is valid.
    if chosen_option not in options_map:
        return f"Error: The chosen option '{chosen_option}' is not one of the valid options (A, B, C, D)."

    # Get the count corresponding to the LLM's chosen option.
    llm_count = options_map[chosen_option]

    # Compare the LLM's count with the correct count.
    if llm_count == correct_count:
        return "Correct"
    else:
        optically_active_compounds = [compound for compound, (is_active, _) in compounds_analysis.items() if is_active]
        
        error_message = (
            f"The answer is incorrect.\n"
            f"The provided answer '{chosen_option}' corresponds to {llm_count} optically active compounds.\n"
            f"The correct number of optically active compounds is {correct_count}.\n\n"
            f"The optically active compounds are:\n"
        )
        for compound in optically_active_compounds:
            error_message += f"- {compound}\n"
            
        return error_message

# Execute the check and print the result.
result = check_answer_correctness()
print(result)