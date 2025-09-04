def check_impulse_approximation_answer():
    """
    Checks the correctness of the answer for the impulse approximation question.

    The function encodes the physical principles of the impulse approximation
    to determine the correct set of assumptions and compares it against the
    provided answer.
    """
    # Define the four assumptions given in the question
    assumptions_text = {
        1: "The interaction current only interacts with individual nucleons.",
        2: "The nucleus is transparent apart from the selected nucleon.",
        3: "The quarks internal to the selected nucleon are non-relativistic.",
        4: "The interaction proceeds as if the selected nucleon experiences no binding forces."
    }

    # Define the available options
    options = {
        'A': {1, 2, 3},
        'B': {1, 3, 4},
        'C': {2, 3, 4},
        'D': {1, 2, 4}
    }

    # The final answer to be checked
    final_answer = 'D'

    # Based on physics principles, the impulse approximation is defined by
    # assumptions 1, 2, and 4. Assumption 3 is about the internal structure
    # of the nucleon, not its behavior within the nucleus under this approximation.
    correct_assumptions_set = {1, 2, 4}

    # Retrieve the set of assumptions for the given final answer
    if final_answer not in options:
        return f"Incorrect. The answer '{final_answer}' is not a valid option."

    selected_assumptions = options[final_answer]

    # Check if the selected assumptions match the correct set
    if selected_assumptions == correct_assumptions_set:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        missing = correct_assumptions_set - selected_assumptions
        extra = selected_assumptions - correct_assumptions_set
        
        reasons = []
        if missing:
            missing_text = ", ".join([f"'{assumptions_text[num]}'" for num in sorted(list(missing))])
            reasons.append(f"The answer is incorrect because it omits a necessary assumption: {missing_text}.")
        
        if extra:
            extra_text = ", ".join([f"'{assumptions_text[num]}'" for num in sorted(list(extra))])
            reasons.append(f"The answer is incorrect because it includes an unnecessary assumption: {extra_text}.")
            
        return "\n".join(reasons)

# Run the check and print the result
result = check_impulse_approximation_answer()
print(result)