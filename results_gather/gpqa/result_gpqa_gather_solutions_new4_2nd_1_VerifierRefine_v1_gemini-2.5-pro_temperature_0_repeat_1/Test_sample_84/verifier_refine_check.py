import math

def check_answer():
    """
    This function checks the correctness of the LLM's answer to the exoplanet temperature ratio problem.
    """
    # Given values from the question
    M_p1 = 7  # Earth masses
    M_p2 = 5  # Earth masses
    delta_lambda1 = 0.03  # Angstroms
    delta_lambda2 = 0.04  # Angstroms

    # The options provided in the question
    options = {
        "A": 0.98,
        "B": 1.05,
        "C": 1.30,
        "D": 0.53
    }

    # The final answer provided by the LLM
    llm_answer_choice = "D"

    # --- Step 1: Perform the correct calculation based on physics principles ---
    # The derivation shows that T_eq1 / T_eq2 = (M_p2 / M_p1) * (delta_lambda1 / delta_lambda2)
    
    # Calculate the ratio of masses
    mass_ratio = M_p2 / M_p1
    
    # Calculate the ratio of Doppler shifts
    doppler_ratio = delta_lambda1 / delta_lambda2
    
    # Calculate the final temperature ratio
    calculated_ratio = mass_ratio * doppler_ratio
    
    # --- Step 2: Check if the LLM's choice is correct ---
    # Get the numerical value corresponding to the LLM's chosen option
    if llm_answer_choice not in options:
        return f"Invalid option '{llm_answer_choice}' selected by the LLM. Valid options are {list(options.keys())}."

    llm_answer_value = options[llm_answer_choice]

    # --- Step 3: Compare the calculated value with the LLM's chosen value ---
    # The options are given with '~', so we check if the calculated value is close to the option value.
    # A tolerance of 0.01 should be sufficient for this problem.
    tolerance = 0.01
    
    if abs(calculated_ratio - llm_answer_value) > tolerance:
        # Find the best matching option for the calculated value
        best_match = min(options.items(), key=lambda item: abs(item[1] - calculated_ratio))
        return (f"Incorrect. The LLM's choice '{llm_answer_choice}' corresponds to a value of {llm_answer_value}. "
                f"The correctly calculated ratio is approximately {calculated_ratio:.4f}. "
                f"This value is closest to option '{best_match[0]}' ({best_match[1]}). "
                f"The LLM selected the wrong option for the calculated value.")

    # --- Step 4: Verify that the chosen option is indeed the best fit ---
    # This confirms that the LLM not only got a close value but picked the *closest* option.
    best_match_choice, best_match_value = min(options.items(), key=lambda item: abs(item[1] - calculated_ratio))
    
    if llm_answer_choice != best_match_choice:
        return (f"Incorrect. While the LLM's choice '{llm_answer_choice}' ({llm_answer_value}) is close to the calculated value ({calculated_ratio:.4f}), "
                f"the best matching option is actually '{best_match_choice}' ({best_match_value}).")

    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)