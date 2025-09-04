import math
import re

def check_correctness_of_answer():
    """
    This function verifies the solution to the quantum mechanics problem.
    1. It calculates the theoretical value of 'a' based on the normalization principle.
    2. It identifies which of the multiple-choice options is the closest.
    3. It compares this correct option with the final answer provided by the LLM.
    """
    # The final answer provided by the LLM to be checked
    llm_final_answer_text = "<<<B>>>"

    # The options from the question
    options = {'A': 0.6, 'B': 0.85, 'C': 1.1, 'D': 0.35}

    # Step 1: Analytically solve for 'a'
    # The normalization integral ∫[1,3] |ψ(x)|² dx = 1 leads to the equation:
    # a² * ln(2) + 0.5 = 1
    # a = sqrt(0.5 / ln(2))
    try:
        a_calculated = math.sqrt(0.5 / math.log(2))
    except Exception as e:
        return f"Error during the analytical calculation of 'a': {e}"

    # Step 2: Find the closest option to the calculated value
    closest_option_label = None
    min_difference = float('inf')

    for label, value in options.items():
        difference = abs(value - a_calculated)
        if difference < min_difference:
            min_difference = difference
            closest_option_label = label

    if closest_option_label is None:
        return "Could not determine the correct option from the choices."

    # Step 3: Extract the LLM's chosen answer
    match = re.search(r'<<<([A-D])>>>', llm_final_answer_text)
    if not match:
        return f"Could not parse the final answer format from the text: '{llm_final_answer_text}'"
    llm_choice = match.group(1)

    # Step 4: Compare the LLM's choice with the correct choice
    if llm_choice == closest_option_label:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_choice}, but the correct option is {closest_option_label}. "
                f"The theoretical value for 'a' is approximately {a_calculated:.4f}, which is closest to "
                f"the value {options[closest_option_label]} (Option {closest_option_label}).")

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)