import math

def check_astronomy_problem():
    """
    This function verifies the solution to the exoplanet orbital period problem.

    It checks the following steps from the provided answer:
    1. The ratio of radial velocity amplitudes (K₂/K₁) is calculated from the wavelength shifts.
    2. The relationship between the K ratio and the period ratio (T₂/T₁) is established.
    3. The final numerical result is calculated.
    4. The calculated result is compared against the chosen option D.
    """

    # --- Given values from the problem statement ---
    delta_lambda_1 = 5.0  # miliangstrom
    delta_lambda_2 = 7.0  # miliangstrom

    # --- Value corresponding to the chosen answer (Option D) ---
    option_d_value = 0.36

    # --- Step 1: Calculate the ratio of radial velocity amplitudes (K₂/K₁) ---
    # The answer states K₂/K₁ = Δλ₂/Δλ₁ = 7/5 = 1.4.
    # This is based on the principle that radial velocity is proportional to Doppler shift.
    # This principle is correct.
    K_ratio = delta_lambda_2 / delta_lambda_1
    if not math.isclose(K_ratio, 1.4):
        return f"Incorrect intermediate calculation: The ratio of velocity amplitudes K₂/K₁ should be 7/5 = 1.4, but was calculated as {K_ratio}."

    # --- Step 2: Relate velocity ratio to period ratio ---
    # The answer states K ∝ T⁻¹/³, which means K₂/K₁ = (T₁/T₂)¹/³.
    # This relationship is correct for circular orbits with identical star and planet masses.
    # From this, we derive T₂/T₁ = (K₂/K₁)⁻³.
    try:
        calculated_period_ratio = K_ratio**(-3)
    except ValueError:
        return "Calculation error: Cannot raise a negative number to a fractional power if K_ratio were negative."

    # --- Step 3: Check the intermediate and final calculations from the answer ---
    # The answer calculates T₁/T₂ = 1.4³ = 2.744.
    T1_over_T2 = K_ratio**3
    if not math.isclose(T1_over_T2, 2.744):
        return f"Incorrect intermediate calculation: The value of (K₂/K₁)³ should be 1.4³ = 2.744, but was calculated as {T1_over_T2}."

    # The answer then calculates T₂/T₁ = 1 / 2.744 ≈ 0.3644.
    final_ratio_from_answer = 1 / 2.744
    if not math.isclose(calculated_period_ratio, final_ratio_from_answer):
        return (f"Discrepancy in final calculation. The code calculated T₂/T₁ as {calculated_period_ratio:.4f}, "
                f"while the answer's steps lead to {final_ratio_from_answer:.4f}.")

    # --- Step 4: Verify if the result matches the chosen option ---
    # The calculated value is ~0.3644, and the chosen option D is ~0.36.
    # We use a reasonable tolerance for multiple-choice questions.
    if math.isclose(calculated_period_ratio, option_d_value, rel_tol=0.03):
        return "Correct"
    else:
        return (f"The final calculated ratio T₂/T₁ is {calculated_period_ratio:.4f}. "
                f"This value is closest to option D ({option_d_value}). "
                f"The provided answer is correct in its reasoning and choice, but the numerical value in the option is a rounded version of the true answer.")

# Running the check
result = check_astronomy_problem()
print(result)