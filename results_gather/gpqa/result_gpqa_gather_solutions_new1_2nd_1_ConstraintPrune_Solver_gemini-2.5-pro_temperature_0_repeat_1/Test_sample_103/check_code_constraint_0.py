import math

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer by recalculating the result
    based on the physical principles outlined in the problem.
    """
    # 1. Define the given parameters from the question.
    delta_lambda_1 = 5.0  # Wavelength shift for system 1 in miliangstrom
    delta_lambda_2 = 7.0  # Wavelength shift for system 2 in miliangstrom

    # 2. Define the multiple-choice options and the provided answer.
    options = {
        'A': 1.40,
        'B': 1.96,
        'C': 0.36,
        'D': 0.85
    }
    provided_answer_key = 'C'
    
    # Check if the provided answer key is valid.
    if provided_answer_key not in options:
        return f"Invalid answer key '{provided_answer_key}'. Valid keys are {list(options.keys())}."
        
    provided_answer_value = options[provided_answer_key]

    # 3. Perform the step-by-step calculation based on physics.
    # Step A: Calculate the ratio of radial velocity semi-amplitudes (K₂ / K₁).
    # K is proportional to Δλ, so K₂ / K₁ = Δλ₂ / Δλ₁.
    try:
        K_ratio_2_over_1 = delta_lambda_2 / delta_lambda_1
    except ZeroDivisionError:
        return "Constraint violation: delta_lambda_1 cannot be zero."

    # Step B: Relate the velocity ratio to the period ratio (T₂ / T₁).
    # K ∝ T^(-1/3), so K₂ / K₁ = (T₁ / T₂)^(1/3).
    # We need to solve for T₂ / T₁.
    # (K₂ / K₁)^3 = T₁ / T₂
    T_ratio_1_over_2 = K_ratio_2_over_1 ** 3
    
    # T₂ / T₁ = 1 / (T₁ / T₂)
    calculated_T_ratio_2_over_1 = 1 / T_ratio_1_over_2

    # 4. Verify the final answer.
    # Check if the calculated result is close to the value of the chosen option.
    # A relative tolerance is used because the options are approximate (~).
    if math.isclose(calculated_T_ratio_2_over_1, provided_answer_value, rel_tol=0.05):
        # The logic in the provided answer is sound:
        # 1. It correctly calculates K₂/K₁ = 1.4.
        # 2. It correctly uses the relationship K ∝ T^(-1/3).
        # 3. It correctly calculates T₂/T₁ = 1 / (1.4)³ ≈ 0.3644.
        # 4. It correctly maps this result to option C.
        return "Correct"
    else:
        return (f"Incorrect. The derivation shows that the ratio T₂/T₁ should be "
                f"1 / (Δλ₂/Δλ₁)^3 = 1 / (7/5)^3 = 1 / (1.4)^3 ≈ {calculated_T_ratio_2_over_1:.4f}. "
                f"This value corresponds to option C ({options['C']}). "
                f"The provided answer was '{provided_answer_key}', which has a value of {provided_answer_value}, "
                f"but the calculation does not match this value within a reasonable tolerance.")

# Execute the check
result = check_answer_correctness()
print(result)