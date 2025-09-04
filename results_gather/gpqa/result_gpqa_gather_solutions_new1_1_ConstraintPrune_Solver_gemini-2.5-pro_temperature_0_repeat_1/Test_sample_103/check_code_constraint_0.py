import math

def check_exoplanet_period_ratio():
    """
    This function verifies the correctness of the LLM's answer by calculating the
    ratio of the orbital periods based on the provided data and physics principles.
    """
    
    # --- Given Data ---
    # Periodic wavelength shift for the star hosting planet #1
    delta_lambda_1 = 5.0  # in miliangstrom
    # Periodic wavelength shift for the star hosting planet #2
    delta_lambda_2 = 7.0  # in miliangstrom

    # --- Options from the question ---
    options = {
        'A': 1.96,
        'B': 0.85,
        'C': 1.40,
        'D': 0.36
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_choice = 'D'

    # --- Calculation based on Physics ---
    # 1. The ratio of radial velocity semi-amplitudes (K) is proportional to the
    #    ratio of wavelength shifts (Δλ).
    #    K₁ / K₂ = Δλ₁ / Δλ₂
    k_ratio_1_over_2 = delta_lambda_1 / delta_lambda_2

    # 2. The relationship between K and orbital period T is K ∝ T⁻¹/³.
    #    This leads to the formula for the period ratio: T₂ / T₁ = (K₁ / K₂)³.
    calculated_period_ratio = k_ratio_1_over_2 ** 3

    # --- Verification ---
    # Find the option that is numerically closest to the calculated result.
    best_choice = None
    min_difference = float('inf')
    
    for choice, value in options.items():
        difference = abs(calculated_period_ratio - value)
        if difference < min_difference:
            min_difference = difference
            best_choice = choice
            
    # Check if the LLM's answer choice matches the calculated best choice.
    if llm_answer_choice == best_choice:
        # The final choice is correct. The reasoning in the provided answer also correctly
        # derives T₂/T₁ = (5/7)³ ≈ 0.3644, which matches option D.
        return "Correct"
    else:
        return (f"Incorrect. The physical derivation shows that the ratio of the orbital periods T₂/T₁ "
                f"is (Δλ₁/Δλ₂)^3. With the given values, this is (5/7)^3 ≈ {calculated_period_ratio:.4f}. "
                f"The closest option is '{best_choice}' ({options[best_choice]}), but the provided answer was "
                f"'{llm_answer_choice}' ({options[llm_answer_choice]}).")

# Run the check
result = check_exoplanet_period_ratio()
print(result)