import math

def check_exoplanet_period_ratio():
    """
    This function checks the correctness of the answer to the exoplanet period ratio question.
    It calculates the expected ratio based on the physics of the radial velocity method
    and compares it to the provided options.
    """

    # --- 1. Define problem constraints and the final answer to check ---
    # Given wavelength shifts in miliangstrom
    delta_lambda_1 = 5.0
    delta_lambda_2 = 7.0

    # The multiple-choice options provided in the question
    options = {
        'A': 0.36,
        'B': 0.85,
        'C': 1.40,
        'D': 1.96
    }

    # The final answer from the LLM to be verified
    llm_answer = 'A'

    # --- 2. Perform the calculation based on physics principles ---
    # The ratio of the orbital periods T₂/T₁ is derived from the relationship:
    # T₂ / T₁ = (Δλ₁ / Δλ₂)^3
    try:
        calculated_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- 3. Determine the correct option ---
    # Find the option key whose value is closest to the calculated ratio
    closest_option = None
    min_difference = float('inf')

    for key, value in options.items():
        difference = abs(value - calculated_ratio)
        if difference < min_difference:
            min_difference = difference
            closest_option = key

    # --- 4. Validate the LLM's answer ---
    if closest_option == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. The calculated ratio T₂/T₁ is approximately {calculated_ratio:.4f}. "
                f"This value is closest to option {closest_option} ({options[closest_option]}). "
                f"The provided answer was '{llm_answer}', which is not the correct choice.")

# Execute the check and print the result
result = check_exoplanet_period_ratio()
print(result)