import math

def check_answer():
    """
    Checks the correctness of the provided answer for the exoplanet orbital period problem.
    """
    # Given values from the question
    delta_lambda_1 = 5  # miliangstrom for planet #1
    delta_lambda_2 = 7  # miliangstrom for planet #2

    # The options provided in the final answer's context
    options = {
        'A': 1.40,
        'B': 0.36,
        'C': 1.96,
        'D': 0.85
    }

    # The final answer choice provided by the LLM
    llm_choice = 'B'

    # --- Step 1: Perform the calculation based on the physics derivation ---
    # The ratio of the orbital periods P2 / P1 is (delta_lambda_1 / delta_lambda_2)^3
    try:
        calculated_ratio = (delta_lambda_1 / delta_lambda_2)**3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Step 2: Check if the calculated value matches the chosen option ---
    if llm_choice not in options:
        return f"The chosen answer '{llm_choice}' is not a valid option."

    expected_value = options[llm_choice]

    # Use a tolerance for floating-point comparison
    if not math.isclose(calculated_ratio, expected_value, rel_tol=1e-2):
        return (f"Incorrect. The calculated ratio P2/P1 is {calculated_ratio:.4f}. "
                f"The chosen option '{llm_choice}' corresponds to a value of {expected_value}, "
                f"which does not match the calculated result.")

    # --- Step 3: Final verification ---
    # The derivation is sound, the calculation is correct, and the result matches the chosen option.
    return "Correct"

# Run the check
result = check_answer()
print(result)