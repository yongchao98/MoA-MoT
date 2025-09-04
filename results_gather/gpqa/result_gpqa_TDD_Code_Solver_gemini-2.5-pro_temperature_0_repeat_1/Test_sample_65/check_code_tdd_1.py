def check_smeft_symmetry_answer():
    """
    Checks the correctness of an answer about required symmetries in SMEFT.

    This function encodes the fundamental principles of SMEFT construction
    to determine which symmetries must be respected by all operators. It then
    compares this ground truth against the provided answer.
    """

    # 1. Define the status of each symmetry based on QFT principles for SMEFT.
    # A symmetry is required if it's a foundational assumption of the theory's construction.
    symmetries_database = {
        1: {
            "name": "Lorentz Symmetry",
            "is_required": True,
            "reason": "SMEFT is a relativistic theory, making Lorentz invariance a fundamental requirement."
        },
        2: {
            "name": "Poincare symmetry",
            "is_required": True,
            "reason": "As a relativistic QFT, SMEFT must be invariant under Poincare transformations (Lorentz transformations + spacetime translations)."
        },
        3: {
            "name": "CP symmetry",
            "is_required": False,
            "reason": "CP is violated in the Standard Model. SMEFT is built to parameterize this and other potential sources of CP violation, so it is not a required symmetry for all operators."
        },
        4: {
            "name": "CPT symmetry",
            "is_required": True,
            "reason": "The CPT theorem applies to any local, Lorentz-invariant QFT like SMEFT, making CPT a mandatory symmetry."
        }
    }

    # 2. Define the question's options and the LLM's answer.
    llm_answer_choice = "D"
    options = {
        "A": {3, 4},
        "B": {1, 3, 4},
        "C": {1, 2},
        "D": {1, 2, 4}
    }

    # 3. Determine the set of correct required symmetries from the database.
    ground_truth_indices = {
        idx for idx, props in symmetries_database.items() if props["is_required"]
    }

    # 4. Get the set of symmetries from the LLM's answer.
    llm_answer_indices = options.get(llm_answer_choice)

    if llm_answer_indices is None:
        return f"Error: The provided answer choice '{llm_answer_choice}' is not a valid option."

    # 5. Compare the LLM's answer with the ground truth and report.
    if llm_answer_indices == ground_truth_indices:
        return "Correct"
    else:
        # Identify which required symmetries were missed.
        missing_symmetries = ground_truth_indices - llm_answer_indices
        error_report = []
        for idx in missing_symmetries:
            name = symmetries_database[idx]['name']
            reason = symmetries_database[idx]['reason']
            error_report.append(
                f"Constraint violated: The answer omits the required symmetry '{name}' (number {idx}). Reason: {reason}"
            )

        # Identify which non-required symmetries were incorrectly included.
        extra_symmetries = llm_answer_indices - ground_truth_indices
        for idx in extra_symmetries:
            name = symmetries_database[idx]['name']
            reason = symmetries_database[idx]['reason']
            error_report.append(
                f"Constraint violated: The answer incorrectly includes '{name}' (number {idx}), which is not a required symmetry. Reason: {reason}"
            )
        
        return "\n".join(error_report)

# Execute the check and print the result.
result = check_smeft_symmetry_answer()
print(result)