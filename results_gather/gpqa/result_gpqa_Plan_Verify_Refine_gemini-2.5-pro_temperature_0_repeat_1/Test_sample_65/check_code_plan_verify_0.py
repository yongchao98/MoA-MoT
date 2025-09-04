def check_smeft_symmetry_answer():
    """
    Checks the correctness of the answer about required symmetries in SMEFT.

    The function encodes the fundamental principles of SMEFT construction:
    1.  It must be a Lorentz- and Poincare-invariant theory.
    2.  It must obey the CPT theorem as a consequence of being a local, Lorentz-invariant QFT.
    3.  It is designed to parameterize new sources of CP violation, so CP is NOT a required symmetry.
    """

    # Define the symmetries and their corresponding numbers from the question
    symmetries = {
        1: "Lorentz Symmetry",
        2: "Poincare symmetry",
        3: "CP symmetry",
        4: "CPT symmetry"
    }

    # Ground truth: A dictionary stating whether a symmetry MUST be respected by ALL SMEFT operators.
    smeft_requirements = {
        "Lorentz Symmetry": True,
        "Poincare symmetry": True,
        "CP symmetry": False,
        "CPT symmetry": True
    }

    # Determine the set of correct symmetry numbers based on the ground truth
    correct_symmetry_numbers = {num for num, name in symmetries.items() if smeft_requirements[name]}

    # Define the options provided in the question
    options = {
        "A": {3, 4},
        "B": {1, 3, 4},
        "C": {1, 2, 4},
        "D": {1, 2}
    }

    # The answer from the LLM to be checked
    llm_answer = "C"

    # Check if the provided answer key exists in the options
    if llm_answer not in options:
        return f"Invalid Answer: The key '{llm_answer}' is not a valid option (A, B, C, or D)."

    # Get the set of symmetries corresponding to the LLM's answer
    llm_selected_symmetries = options[llm_answer]

    # Compare the LLM's answer with the correct set of symmetries
    if llm_selected_symmetries == correct_symmetry_numbers:
        return "Correct"
    else:
        # Build a detailed reason for the incorrectness
        reason = "Incorrect.\n"
        
        # Check for missing required symmetries
        missing_symmetries = correct_symmetry_numbers - llm_selected_symmetries
        if missing_symmetries:
            missing_names = [f"'{symmetries[num]}' ({num})" for num in sorted(list(missing_symmetries))]
            reason += f"The answer fails to include required symmetries: {', '.join(missing_names)}.\n"

        # Check for incorrectly included symmetries
        incorrectly_included_symmetries = llm_selected_symmetries - correct_symmetry_numbers
        if incorrectly_included_symmetries:
            incorrect_names = [f"'{symmetries[num]}' ({num})" for num in sorted(list(incorrectly_included_symmetries))]
            reason += f"The answer incorrectly includes symmetries that are not required for all SMEFT operators: {', '.join(incorrect_names)}.\n"

        # Provide the correct reasoning
        reason += "\nExplanation:\n"
        reason += "- Lorentz (1) and Poincare (2) symmetries are fundamental spacetime symmetries that any relativistic QFT, including SMEFT, must respect.\n"
        reason += "- CPT symmetry (4) is a direct consequence of the CPT theorem, which applies to all local, Lorentz-invariant QFTs.\n"
        reason += "- CP symmetry (3) is NOT a required symmetry. SMEFT is explicitly constructed to include operators that can violate CP symmetry, parameterizing new physics beyond the Standard Model's known CP violation."
        
        return reason

# Execute the check and print the result
print(check_smeft_symmetry_answer())