import re

def check_impulse_approximation_answer(llm_response: str) -> str:
    """
    Checks the correctness of the answer for the impulse approximation question.

    Args:
        llm_response: The full response from the LLM, including the final answer
                      in the format <<<...>>>.

    Returns:
        "Correct" if the answer is right, otherwise a string explaining the error.
    """
    # Extract the final answer from the response string
    match = re.search(r'<<<(.+?)>>>', llm_response)
    if not match:
        return "The response does not contain a final answer in the required format <<<...>>>."

    answer = match.group(1).strip().upper()

    # Define the assumptions for each option
    options = {
        'A': {1, 2, 4},
        'B': {1, 3, 4},
        'C': {1, 2, 3},
        'D': {2, 3, 4}
    }

    # The correct set of assumptions for the impulse approximation
    correct_assumptions = {1, 2, 4}
    correct_option = 'A'

    # Check if the provided answer is a valid option
    if answer not in options:
        return f"Invalid option '{answer}'. The answer must be one of A, B, C, or D."

    # Check if the answer is correct
    if answer == correct_option:
        return "Correct"
    else:
        # Analyze why the chosen answer is incorrect
        chosen_assumptions = options[answer]
        
        missing = correct_assumptions - chosen_assumptions
        incorrectly_included = chosen_assumptions - correct_assumptions
        
        reasons = []
        if 3 in incorrectly_included:
            reasons.append("it incorrectly includes assumption 3. The impulse approximation is about the nucleon's interaction within the nucleus, not its internal quark structure.")
        if 1 in missing:
            reasons.append("it is missing assumption 1, which is the core idea of reducing the problem to a single-nucleon interaction.")
        if 2 in missing:
            reasons.append("it is missing assumption 2 (no final-state interactions), which is needed to treat the nucleon as isolated after the interaction.")
        if 4 in missing:
            reasons.append("it is missing assumption 4, which allows ignoring the nuclear binding potential during the brief interaction.")

        return f"The answer '{answer}' is incorrect because {', and '.join(reasons)}. The correct answer is '{correct_option}' (assumptions 1, 2, and 4)."

# The user-provided response to be checked
llm_response_to_check = """
Here are the tests for the impulse approximation problem.

1.  **`test_correct_assumptions`**: The solution must identify the combination of assumptions (1, 2, and 4) as the correct set, which corresponds to option 'A'. The impulse approximation requires treating the interaction as occurring with a single nucleon (1), which is then ejected without further interaction (2), and the interaction itself is so rapid that binding forces are negligible during the event (4).
2.  **`test_excludes_internal_structure_assumption`**: The solution must exclude assumption 3. The internal quark structure of a nucleon (relativistic or not) is a separate physical consideration from how the nucleon as a whole interacts within the nucleus. Options B, C, and D incorrectly include this assumption.
3.  **`test_includes_one_body_current_assumption`**: The solution must include assumption 1. The core idea of the approximation is to reduce a many-body problem to a one-body problem, so the interaction must be with an individual nucleon. Option D incorrectly omits this.
4.  **`test_output_format`**: The function must return a single, uppercase character string, corresponding to one of the options A, B, C, or D.
