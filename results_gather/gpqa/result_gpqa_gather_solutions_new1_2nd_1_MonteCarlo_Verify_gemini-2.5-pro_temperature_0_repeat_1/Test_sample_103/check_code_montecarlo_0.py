import math

def check_answer_correctness():
    """
    Checks the correctness of the final answer for the exoplanet period ratio problem.

    The function follows the physical derivation:
    1. The radial velocity semi-amplitude (K) is proportional to the wavelength shift (Δλ).
       Therefore, K₂ / K₁ = Δλ₂ / Δλ₁ = 7 / 5 = 1.4.
    2. For constant stellar and planetary minimum masses (mₚ sin i), the relationship between K and the orbital period (T) is K ∝ T^(-1/3).
    3. This leads to the ratio: K₂ / K₁ = (T₂ / T₁)^(-1/3).
    4. Solving for the desired ratio T₂ / T₁ gives: T₂ / T₁ = (K₂ / K₁)^(-3) = (7 / 5)^(-3) = (5 / 7)³.
    5. The calculated value is then compared against the provided answer.
    """
    # Given values from the question
    delta_lambda_1 = 5.0
    delta_lambda_2 = 7.0

    # The options provided in the question
    options = {
        "A": 1.40,
        "B": 0.36,
        "C": 1.96,
        "D": 0.85
    }
    
    # The final answer from the LLM to be checked
    llm_final_answer_letter = "B"

    # --- Calculation ---
    # Step 1: Calculate the ratio of radial velocity semi-amplitudes (K₂ / K₁)
    try:
        k_ratio = delta_lambda_2 / delta_lambda_1
    except ZeroDivisionError:
        return "Incorrect: The wavelength shift for planet #1 cannot be zero."

    # Step 2: Calculate the period ratio (T₂ / T₁) using the physical relationship
    # T₂ / T₁ = (K₂ / K₁)^(-3)
    calculated_period_ratio = k_ratio ** (-3)

    # --- Verification ---
    # Check if the provided answer letter is a valid option
    if llm_final_answer_letter not in options:
        return f"Incorrect: The provided answer '{llm_final_answer_letter}' is not one of the valid options (A, B, C, D)."

    # Get the numerical value of the LLM's answer
    llm_answer_value = options[llm_final_answer_letter]

    # Check if the calculated value is close to the LLM's chosen option value.
    # A relative tolerance of 5% is reasonable for approximate answers.
    if not math.isclose(calculated_period_ratio, llm_answer_value, rel_tol=0.05):
        return (f"Incorrect: The calculated period ratio T₂/T₁ is approximately {calculated_period_ratio:.4f}. "
                f"The provided answer is '{llm_final_answer_letter}' which corresponds to a value of {llm_answer_value}. "
                f"These values do not match.")

    # As a final check, ensure that the chosen option is indeed the best fit among all options.
    best_fit_option = None
    min_error = float('inf')

    for option_letter, option_value in options.items():
        error = abs(calculated_period_ratio - option_value)
        if error < min_error:
            min_error = error
            best_fit_option = option_letter
            
    if best_fit_option != llm_final_answer_letter:
        return (f"Incorrect: The provided answer is '{llm_final_answer_letter}' ({llm_answer_value}). "
                f"While this value is close to the calculated ratio of {calculated_period_ratio:.4f}, "
                f"option '{best_fit_option}' ({options[best_fit_option]}) is a better numerical match.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer_correctness()
print(result)