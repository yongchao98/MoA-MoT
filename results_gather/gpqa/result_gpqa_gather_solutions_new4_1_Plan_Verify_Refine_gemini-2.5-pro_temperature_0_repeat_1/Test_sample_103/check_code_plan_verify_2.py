import math

def check_exoplanet_period_ratio():
    """
    Checks the correctness of the provided answer for the exoplanet period ratio question.
    """
    # --- Problem Setup ---
    # Given values from the question
    delta_lambda_1 = 5.0  # miliangstrom for planet #1
    delta_lambda_2 = 7.0  # miliangstrom for planet #2

    # Options from the question prompt
    # A) ~ 1.40
    # B) ~ 0.36
    # C) ~ 1.96
    # D) ~ 0.85
    options = {'A': 1.40, 'B': 0.36, 'C': 1.96, 'D': 0.85}
    
    # The answer from the LLM to be checked
    llm_answer_choice = 'A'

    # --- Correct Derivation ---
    # 1. The radial velocity semi-amplitude (K) is proportional to the wavelength shift (Δλ).
    #    Therefore, K₂ / K₁ = Δλ₂ / Δλ₁
    # 2. For the given constraints (same star mass, similar planet mass), K is proportional to P^(-1/3).
    #    K ∝ P^(-1/3)
    # 3. From this, we can derive the ratio of the periods:
    #    K₂ / K₁ = (P₁ / P₂)^(1/3)
    #    (K₂ / K₁)^3 = P₁ / P₂
    #    P₂ / P₁ = (K₁ / K₂)^3
    # 4. Since K₁/K₂ = Δλ₁/Δλ₂, the final formula is:
    #    P₂ / P₁ = (Δλ₁ / Δλ₂)^3
    
    correct_period_ratio = (delta_lambda_1 / delta_lambda_2) ** 3

    # --- Verification ---
    llm_answer_value = options.get(llm_answer_choice)
    if llm_answer_value is None:
        return f"Invalid answer choice '{llm_answer_choice}'. Valid options are A, B, C, D."

    # Find which option corresponds to the correct calculation
    tolerance = 0.01
    correct_option = None
    for option_key, option_value in options.items():
        if math.isclose(correct_period_ratio, option_value, rel_tol=tolerance, abs_tol=tolerance):
            correct_option = option_key
            break
    
    if llm_answer_choice == correct_option:
        return "Correct"
    else:
        # Check for a common mistake: providing the velocity ratio K₂/K₁
        velocity_ratio_2_over_1 = delta_lambda_2 / delta_lambda_1
        if math.isclose(llm_answer_value, velocity_ratio_2_over_1, rel_tol=tolerance, abs_tol=tolerance):
            reason = (f"Incorrect. The provided answer '{llm_answer_choice}' corresponds to the value ~{llm_answer_value}, "
                      f"which is the ratio of the velocity amplitudes (K₂/K₁ = 7/5 = 1.4). "
                      f"The question asks for the ratio of the orbital periods (P₂/P₁). "
                      f"The correct calculation is P₂/P₁ = (K₁/K₂)³ = (5/7)³ ≈ 0.36. "
                      f"This corresponds to option B.")
            return reason
        
        # Generic error message if it's not a common mistake
        reason = (f"Incorrect. The provided answer is '{llm_answer_choice}' (~{llm_answer_value}). "
                  f"The correct answer is option '{correct_option}' (~{options[correct_option]}), "
                  f"based on the calculation P₂/P₁ = (5/7)³ ≈ {correct_period_ratio:.2f}.")
        return reason

# Execute the check and print the result
result = check_exoplanet_period_ratio()
print(result)