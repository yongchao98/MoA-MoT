def check_smeft_symmetry_answer():
    """
    Checks the correctness of an answer about the required symmetries in SMEFT.

    This function encodes the fundamental principles of the Standard Model Effective
    Field Theory (SMEFT) as a ground truth and compares the given answer against it.
    """

    # Step 1: Define the ground truth based on the established principles of SMEFT.
    # A symmetry is 'True' if it MUST be respected by all operators in SMEFT.
    # A symmetry is 'False' if it is NOT a required symmetry for all operators.
    ground_truth = {
        "Lorentz Symmetry": True,
        "Poincare symmetry": True,
        "CP symmetry": False,
        "CPT symmetry": True
    }

    # Step 2: Map the question's numbering to the symmetry names for clarity.
    symmetry_map = {
        1: "Lorentz Symmetry",
        2: "Poincare symmetry",
        3: "CP symmetry",
        4: "CPT symmetry"
    }

    # Step 3: Define the options provided in the multiple-choice question.
    options = {
        'A': {1, 2},
        'B': {1, 2, 4},
        'C': {3, 4},
        'D': {1, 3, 4}
    }

    # The answer from the LLM to be checked.
    llm_answer = 'B'

    # Step 4: Perform the verification.

    # Determine the set of required symmetries based on the ground truth.
    required_symmetries_indices = {
        idx for idx, name in symmetry_map.items() if ground_truth.get(name, False)
    }

    # Get the set of symmetries corresponding to the LLM's answer.
    llm_chosen_indices = options.get(llm_answer)

    if llm_chosen_indices is None:
        return f"Incorrect. The provided answer '{llm_answer}' is not a valid option."

    # Compare the LLM's answer with the ground truth.
    if llm_chosen_indices == required_symmetries_indices:
        return "Correct"
    else:
        # Generate a detailed error report if the answer is incorrect.
        error_messages = []
        
        # Check for required symmetries that were missed by the answer.
        missing_symmetries = required_symmetries_indices - llm_chosen_indices
        if missing_symmetries:
            for idx in sorted(list(missing_symmetries)):
                name = symmetry_map[idx]
                reason = ""
                if name == "Poincare symmetry" or name == "Lorentz Symmetry":
                    reason = "The SMEFT Lagrangian is constructed to be invariant under the Poincaré group (which includes Lorentz transformations) by definition."
                elif name == "CPT symmetry":
                    reason = "The CPT theorem requires any local, unitary, Lorentz-invariant QFT to be CPT invariant. SMEFT fulfills these conditions."
                error_messages.append(f"The answer fails to include the required symmetry: '{name}'. Reason: {reason}")

        # Check for symmetries that were incorrectly included in the answer.
        incorrectly_included_symmetries = llm_chosen_indices - required_symmetries_indices
        if incorrectly_included_symmetries:
            for idx in sorted(list(incorrectly_included_symmetries)):
                name = symmetry_map[idx]
                reason = ""
                if name == "CP symmetry":
                    reason = "CP symmetry is violated in the Standard Model and SMEFT is a framework used to parameterize new sources of CP violation. Therefore, it is not a required symmetry."
                error_messages.append(f"The answer incorrectly includes the non-required symmetry: '{name}'. Reason: {reason}")
        
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result.
result = check_smeft_symmetry_answer()
print(result)