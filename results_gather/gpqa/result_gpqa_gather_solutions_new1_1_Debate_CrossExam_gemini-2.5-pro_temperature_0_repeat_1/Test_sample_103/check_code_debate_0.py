import math

def check_exoplanet_period_ratio():
    """
    This function verifies the calculation for the ratio of orbital periods of two exoplanets.

    The core physics principles are:
    1. The radial velocity semi-amplitude (K) is proportional to the observed wavelength shift (Δλ).
       K ∝ Δλ  =>  K₂ / K₁ = Δλ₂ / Δλ₁
    2. For a circular orbit, K is related to the orbital period (T) by:
       K ∝ T^(-1/3), assuming stellar and planetary masses are constant between the two systems.
    3. Combining these gives: Δλ₂ / Δλ₁ = (T₁ / T₂)^(1/3)
    4. The question asks for T₂ / T₁, which can be solved from the above equation.
    """

    # Given values from the question
    delta_lambda_1 = 5.0  # miliangstrom
    delta_lambda_2 = 7.0  # miliangstrom

    # The final answer provided by the LLM to be checked
    llm_final_answer_choice = 'A'
    llm_final_answer_value = 0.36

    # --- Step 1: Calculate the ratio of the radial velocity semi-amplitudes ---
    # K₂ / K₁ = Δλ₂ / Δλ₁
    K_ratio_2_over_1 = delta_lambda_2 / delta_lambda_1

    # --- Step 2: Use the physical relationship to find the period ratio ---
    # K₂ / K₁ = (T₁ / T₂)^(1/3)
    # We need to solve for T₂ / T₁
    
    # First, solve for T₁ / T₂ by cubing both sides
    T_ratio_1_over_2 = K_ratio_2_over_1 ** 3
    
    # Then, take the reciprocal to find T₂ / T₁
    calculated_T_ratio_2_over_1 = 1 / T_ratio_1_over_2

    # --- Step 3: Verify the LLM's answer ---
    # The LLM's reasoning is that T₂ / T₁ ≈ 0.3644, which corresponds to option A.
    
    # Check if the LLM's chosen option matches the calculation.
    options = {'A': 0.36, 'B': 0.85, 'C': 1.96, 'D': 1.40}
    
    # Find which option is numerically closest to our calculated result.
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_T_ratio_2_over_1))

    # Check if the LLM's final choice is the correct one.
    if llm_final_answer_choice != closest_option:
        return (f"Incorrect. The calculated ratio T₂/T₁ is approximately {calculated_T_ratio_2_over_1:.4f}. "
                f"This value is closest to option {closest_option} ({options[closest_option]}), "
                f"but the provided answer chose option {llm_final_answer_choice}.")

    # Check if the derivation logic in the provided answer is sound.
    # The provided answer correctly derives that T₂/T₁ = 1 / (1.4)³ ≈ 0.3644.
    # This matches our calculation.
    if not math.isclose(calculated_T_ratio_2_over_1, 1 / (1.4**3), rel_tol=1e-9):
        return (f"Incorrect. There is a flaw in the calculation logic. "
                f"The code calculated {calculated_T_ratio_2_over_1:.4f}, which does not match the expected 1 / (1.4**3).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_exoplanet_period_ratio()
print(result)