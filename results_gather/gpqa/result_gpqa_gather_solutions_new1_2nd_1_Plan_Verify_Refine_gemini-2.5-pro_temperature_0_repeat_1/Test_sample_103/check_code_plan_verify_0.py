import math

def check_answer():
    """
    This function verifies the correctness of the answer to the exoplanet period problem.

    The core physics principles are:
    1. The radial velocity semi-amplitude (K) is directly proportional to the observed wavelength shift (Δλ).
       Therefore, K₂ / K₁ = Δλ₂ / Δλ₁.
    2. For a circular orbit with constant stellar and planetary minimum masses, K is proportional to T^(-1/3), where T is the orbital period.
       Therefore, K₂ / K₁ = (T₂ / T₁)^(-1/3) = (T₁ / T₂)^(1/3).

    Combining these gives: Δλ₂ / Δλ₁ = (T₁ / T₂)^(1/3).
    We can solve for the required ratio T₂ / T₁.
    """

    # Given values from the problem
    delta_lambda_1 = 5.0  # miliangstrom
    delta_lambda_2 = 7.0  # miliangstrom

    # The multiple-choice options provided in the question
    options = {
        'A': 1.40,
        'B': 1.96,
        'C': 0.36,
        'D': 0.85
    }

    # The final answer to be checked
    final_answer_letter = 'C'

    # Step 1: Calculate the ratio of the radial velocity semi-amplitudes (K)
    k_ratio_2_over_1 = delta_lambda_2 / delta_lambda_1

    # Step 2: Use the physical relationship to solve for the period ratio (T₂ / T₁)
    # k_ratio_2_over_1 = (T₁ / T₂)^(1/3)
    # (k_ratio_2_over_1)³ = T₁ / T₂
    # T₂ / T₁ = 1 / (k_ratio_2_over_1)³
    
    calculated_period_ratio = 1 / (k_ratio_2_over_1 ** 3)

    # Step 3: Check if the calculated value matches the provided answer
    
    # Get the numerical value of the provided answer
    provided_answer_value = options.get(final_answer_letter)
    if provided_answer_value is None:
        return f"Error: The provided answer letter '{final_answer_letter}' is not a valid option."

    # Compare the calculated result with the value of the chosen option.
    # A tolerance is used because the options are given with "approximate" signs (~).
    # A relative tolerance of 5% is reasonable for this context.
    if math.isclose(calculated_period_ratio, provided_answer_value, rel_tol=0.05):
        return "Correct"
    else:
        # Find the best matching option for the calculated result
        best_match_letter = min(options, key=lambda k: abs(options[k] - calculated_period_ratio))
        
        reason = (
            f"Incorrect. The physical derivation shows that the ratio of the orbital periods (T₂ / T₁) should be "
            f"(Δλ₁ / Δλ₂)^3 = (5/7)^3 ≈ {calculated_period_ratio:.4f}. "
            f"This calculated value corresponds to option C (~0.36). "
            f"The provided answer was {final_answer_letter}, which is consistent with the calculation."
        )
        # This part of the logic is a bit tricky. The code is checking if the provided answer *is* correct.
        # If the calculation matches the provided answer, it should return "Correct".
        # If it doesn't, it should explain why. Let's refine the reason.
        
        reason_refined = (
            f"Incorrect. The calculated ratio of the orbital periods (T₂ / T₁) is approximately {calculated_period_ratio:.4f}. "
            f"The provided answer is {final_answer_letter}, which corresponds to a value of {provided_answer_value:.2f}. "
            f"The calculated value does not match the provided answer's value within a 5% tolerance."
        )
        # Since the provided answer C (0.36) IS the correct one, the code should return "Correct".
        # The logic above is correct. Let's re-verify the check.
        # calculated_period_ratio = 1 / (1.4**3) = 1 / 2.744 = 0.3644...
        # provided_answer_value = 0.36
        # math.isclose(0.3644, 0.36, rel_tol=0.05) -> True, because |0.3644 - 0.36| <= 0.05 * 0.36
        # 0.0044 <= 0.018, which is true.
        # So the code will correctly return "Correct".
        return reason_refined

# The code block to be returned
print(check_answer())