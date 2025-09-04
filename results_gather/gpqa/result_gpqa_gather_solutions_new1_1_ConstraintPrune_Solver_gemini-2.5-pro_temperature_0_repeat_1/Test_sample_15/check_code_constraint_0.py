def check_answer():
    """
    Checks the correctness of the provided answer about optically active compounds.
    """
    # A compound is optically active if it is chiral.
    # A compound is achiral if it has an internal plane of symmetry or is a meso compound.
    # We encode the chemical analysis of each compound into a dictionary.
    compounds = [
        {
            "name": "(Z)-1-chloro-2-methylbut-1-ene",
            "is_chiral": False, # Planar molecule, has a plane of symmetry
            "is_meso": False,
        },
        {
            "name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
            "is_chiral": True, # Name specifies a single enantiomer of a complex, asymmetric molecule
            "is_meso": False,
        },
        {
            "name": "(2R,3S)-2,3-dimethylsuccinic acid",
            "is_chiral": False, # This is a meso compound, achiral due to internal symmetry
            "is_meso": True,
        },
        {
            "name": "(2R,3R)-2,3-dimethylsuccinic acid",
            "is_chiral": True, # This is a chiral diastereomer, specified as a single enantiomer
            "is_meso": False,
        },
        {
            "name": "(R)-cyclohex-3-en-1-ol",
            "is_chiral": True, # Name specifies a single enantiomer of a chiral molecule
            "is_meso": False,
        },
        {
            "name": "(1s,3s,5s)-cyclohexane-1,3,5-triol",
            "is_chiral": False, # All-cis configuration is highly symmetric, has planes of symmetry
            "is_meso": False,
        },
        {
            "name": "1-cyclopentyl-3-methylbutan-1-one",
            "is_chiral": False, # No chiral centers, has a plane of symmetry
            "is_meso": False,
        }
    ]

    # The final answer from the LLM
    llm_answer_choice = "D"
    llm_reasoning_active_list = [
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
        "(2R,3R)-2,3-dimethylsuccinic acid",
        "(R)-cyclohex-3-en-1-ol"
    ]

    # Map options to counts
    option_map = {'A': 5, 'B': 2, 'C': 4, 'D': 3}
    expected_count = option_map.get(llm_answer_choice)

    if expected_count is None:
        return f"Invalid answer choice '{llm_answer_choice}'. Options are A, B, C, D."

    # Determine optical activity based on our encoded properties
    code_determined_active_compounds = []
    for c in compounds:
        # A compound is optically active if it is chiral.
        if c["is_chiral"]:
            code_determined_active_compounds.append(c["name"])
    
    actual_count = len(code_determined_active_compounds)

    # Check 1: The count must match the answer choice.
    if actual_count != expected_count:
        return (f"Incorrect. The number of optically active compounds is {actual_count}, "
                f"but the answer choice '{llm_answer_choice}' corresponds to {expected_count}.")

    # Check 2: The list of identified compounds must match the reasoning.
    if set(code_determined_active_compounds) != set(llm_reasoning_active_list):
        return (f"Incorrect. The reasoning identifies a different set of active compounds.\n"
                f"Expected: {sorted(llm_reasoning_active_list)}\n"
                f"Found: {sorted(code_determined_active_compounds)}")

    return "Correct"

# Run the check
result = check_answer()
print(result)