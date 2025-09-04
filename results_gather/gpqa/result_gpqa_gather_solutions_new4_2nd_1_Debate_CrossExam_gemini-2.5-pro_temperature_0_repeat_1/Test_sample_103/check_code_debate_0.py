import math

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer to the exoplanet orbital period question.
    """
    # --- Problem Constraints & Given Data ---
    # Wavelength shift for planet #1's host star
    delta_lambda_1 = 5  # miliangstrom
    # Wavelength shift for planet #2's host star
    delta_lambda_2 = 7  # miliangstrom

    # --- LLM's Final Answer ---
    # The final answer selected by the LLM
    llm_selected_option = "A"
    # The options as listed in the final answer block
    options_map = {
        "A": 0.36,
        "B": 1.96,
        "C": 0.85,
        "D": 1.40
    }

    # --- Verification Calculation ---
    # The relationship between the radial velocity semi-amplitude (K) and wavelength shift (Δλ) is K ∝ Δλ.
    # The relationship between K and orbital period (P) is K ∝ P^(-1/3).
    # From these, we derive the relationship for the ratio of periods: P₂ / P₁ = (K₁ / K₂)^3.
    # Since K is proportional to Δλ, we can substitute: P₂ / P₁ = (Δλ₁ / Δλ₂)^3.
    
    try:
        # Calculate the ratio of the wavelength shifts
        ratio_delta_lambda = delta_lambda_1 / delta_lambda_2
        
        # Calculate the final ratio of the periods
        calculated_ratio_p2_p1 = ratio_delta_lambda ** 3
        
        # Get the value corresponding to the LLM's chosen option
        llm_answer_value = options_map.get(llm_selected_option)
        
        if llm_answer_value is None:
            return f"Incorrect. The chosen option '{llm_selected_option}' is not a valid choice from the provided options."

        # Check if the calculated value is approximately equal to the option's value
        # A tolerance of 0.01 is reasonable given the '~' symbol in the options.
        if math.isclose(calculated_ratio_p2_p1, llm_answer_value, rel_tol=0.0, abs_tol=0.01):
            return "Correct"
        else:
            return (f"Incorrect. The final answer is wrong.\n"
                    f"The correct derivation is P₂/P₁ = (Δλ₁/Δλ₂)^3.\n"
                    f"The calculation gives: P₂/P₁ = (5/7)^3 ≈ {calculated_ratio_p2_p1:.4f}.\n"
                    f"The chosen option '{llm_selected_option}' corresponds to the value {llm_answer_value}.\n"
                    f"The calculated value {calculated_ratio_p2_p1:.4f} does not match the value of the chosen option {llm_answer_value} within a reasonable tolerance.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check and print the result
result = check_answer_correctness()
print(result)