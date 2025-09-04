def check_impulse_approximation_answer():
    """
    Checks the correctness of the selected option for the impulse approximation question.
    """
    # Define the statements given in the question
    statements = {
        1: "The interaction current only interacts with individual nucleons.",
        2: "The nucleus is transparent apart from the selected nucleon.",
        3: "The quarks internal to the selected nucleon are non-relativistic.",
        4: "The interaction proceeds as if the selected nucleon experiences no binding forces."
    }

    # Define the options and the statements they correspond to
    options = {
        "A": {2, 3, 4},
        "B": {1, 2, 3},
        "C": {1, 3, 4},
        "D": {1, 2, 4}
    }

    # The provided answer from the LLM
    llm_answer = "D"

    # Define the correct set of assumptions for the impulse approximation based on its definition in physics.
    # 1. Interaction with a single nucleon.
    # 2. No final state interactions (nucleus is transparent).
    # 3. Interaction is fast, so binding forces are ignored during the event.
    correct_assumptions_indices = {1, 2, 4}

    # Retrieve the set of assumptions for the given answer
    answer_indices = options.get(llm_answer)

    if answer_indices is None:
        return f"The provided answer '{llm_answer}' is not a valid option."

    # Compare the answer's assumptions with the correct ones
    if answer_indices == correct_assumptions_indices:
        return "Correct"
    else:
        # Identify what is wrong with the answer
        missing = correct_assumptions_indices - answer_indices
        extra = answer_indices - correct_assumptions_indices

        reason = "The answer is incorrect.\n"
        reason += "The impulse approximation is defined by assumptions 1, 2, and 4.\n"
        
        if missing:
            reason += f"The selected answer is missing the required assumption(s): {sorted(list(missing))}.\n"
            for index in sorted(list(missing)):
                 reason += f"  - Assumption {index}: '{statements[index]}'\n"

        if extra:
            reason += f"The selected answer incorrectly includes assumption(s): {sorted(list(extra))}.\n"
            for index in sorted(list(extra)):
                 reason += f"  - Assumption {index}: '{statements[index]}' is not a defining feature of the impulse approximation. It relates to the internal structure of the nucleon, not its interaction within the nucleus.\n"
        
        return reason.strip()

# Run the check
result = check_impulse_approximation_answer()
print(result)