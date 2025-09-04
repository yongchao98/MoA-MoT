import math

def check_exoplanet_period_ratio():
    """
    Checks the correctness of the LLM's answer for the exoplanet period ratio problem.

    The core logic is:
    1. The radial velocity semi-amplitude (K) is proportional to the wavelength shift (Δλ).
       Therefore, K₁/K₂ = Δλ₁/Δλ₂.
    2. For the given conditions (same star/planet masses, same inclination), the orbital period (P)
       is proportional to K⁻³.
    3. This leads to the final relationship: P₂/P₁ = (K₁/K₂)³ = (Δλ₁/Δλ₂)³.
    """
    
    # Given values from the question
    delta_lambda_1 = 5.0  # miliangstrom
    delta_lambda_2 = 7.0  # miliangstrom

    # The options provided in the question text
    options = {
        'A': 0.36,
        'B': 1.40,
        'C': 1.96,
        'D': 0.85
    }

    # The final answer chosen by the LLM
    llm_choice = 'A'

    # Perform the calculation based on the physics
    try:
        # P₂/P₁ = (Δλ₁/Δλ₂)³
        calculated_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # Get the value corresponding to the LLM's chosen option
    chosen_option_value = options.get(llm_choice)

    if chosen_option_value is None:
        return f"The chosen option '{llm_choice}' is invalid. Valid options are {list(options.keys())}."

    # Check if the calculated result is approximately equal to the chosen option's value.
    # The tilde (~) in the options indicates an approximate value, so we use a tolerance.
    # A tolerance of 0.01 is reasonable for this context.
    if math.isclose(calculated_ratio, chosen_option_value, rel_tol=0.02, abs_tol=0.02):
        return "Correct"
    else:
        # Find the correct option based on the calculation
        correct_choice = None
        for option_key, option_value in options.items():
            if math.isclose(calculated_ratio, option_value, rel_tol=0.02, abs_tol=0.02):
                correct_choice = option_key
                break
        
        reason = (
            f"The answer is incorrect. "
            f"The calculation for the period ratio P₂/P₁ is (Δλ₁/Δλ₂)³ = (5/7)³ ≈ {calculated_ratio:.4f}. "
            f"The chosen option '{llm_choice}' corresponds to a value of {chosen_option_value}. "
            f"The calculated value {calculated_ratio:.4f} does not match the value of the chosen option. "
        )
        if correct_choice:
            reason += f"The correct option is '{correct_choice}', which corresponds to a value of ~{options[correct_choice]}."
        else:
            reason += "None of the options match the calculated value."
            
        return reason

# Execute the check and print the result
print(check_exoplanet_period_ratio())