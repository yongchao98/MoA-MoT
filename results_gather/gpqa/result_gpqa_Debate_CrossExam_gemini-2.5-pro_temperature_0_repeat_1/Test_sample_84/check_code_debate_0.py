import math

def check_exoplanet_temperature_ratio():
    """
    Checks the correctness of the LLM's answer for the exoplanet temperature ratio problem.
    """
    # --- 1. Define given values from the question ---
    M_p1 = 7.0  # in Earth masses
    M_p2 = 5.0  # in Earth masses
    delta_lambda1 = 0.03  # in Angstroms
    delta_lambda2 = 0.04  # in Angstroms

    # --- 2. Replicate the derivation from the LLM's answer ---

    # The ratio of equilibrium temperatures T_eq1 / T_eq2 simplifies to sqrt(a2 / a1)
    # because T_star, R_star, and albedo are the same for both planets.
    # T_eq ∝ sqrt(1/a) => T_eq1 / T_eq2 = sqrt(a2/a1)

    # The radial velocity semi-amplitude K is related to the semi-major axis 'a' by:
    # K ∝ M_p / sqrt(a)
    # Rearranging for 'a': a ∝ (M_p / K)^2
    # Therefore, the ratio a2 / a1 is:
    # a2 / a1 = (M_p2 / K2)^2 / (M_p1 / K1)^2
    # a2 / a1 = (M_p2 / M_p1)^2 * (K1 / K2)^2

    # The ratio of RV semi-amplitudes K1 / K2 is equal to the ratio of Doppler shifts.
    # K = c * (Δλ / λ) => K1 / K2 = Δλ1 / Δλ2
    
    # --- 3. Perform the calculations ---

    # Calculate the ratio of planet masses (M_p2 / M_p1)
    mass_ratio_p2_over_p1 = M_p2 / M_p1
    
    # Calculate the ratio of RV semi-amplitudes (K1 / K2)
    K_ratio_1_over_2 = delta_lambda1 / delta_lambda2
    
    # Calculate the ratio of semi-major axes (a2 / a1)
    axis_ratio_a2_over_a1 = (mass_ratio_p2_over_p1**2) * (K_ratio_1_over_2**2)
    
    # Calculate the final ratio of equilibrium temperatures (T_eq1 / T_eq2)
    temp_ratio_t1_over_t2 = math.sqrt(axis_ratio_a2_over_a1)

    # --- 4. Verify the answer ---
    
    # The LLM's answer provides a step-by-step calculation. Let's check the final value.
    # LLM's calculation: sqrt((5/7)^2 * (3/4)^2) = sqrt(225/784) = 15/28
    expected_value = 15.0 / 28.0
    
    # The LLM's final numerical result is ~0.5357, which corresponds to option C (~0.53).
    llm_choice = 'C'
    llm_value = 0.5357

    # Check if our calculated value is close to the expected value.
    # A relative tolerance of 1e-5 is sufficient for this check.
    if not math.isclose(temp_ratio_t1_over_t2, expected_value, rel_tol=1e-5):
        return (f"Calculation mismatch. The derived logic leads to a value of {expected_value:.4f}, "
                f"but the code calculated {temp_ratio_t1_over_t2:.4f}.")

    # Check if the calculated value corresponds to the chosen option C.
    # Option C is ~0.53. Our value is ~0.5357. This is a good match.
    if llm_choice != 'C':
        return (f"Incorrect option chosen. The calculated ratio is {temp_ratio_t1_over_t2:.4f}, "
                f"which corresponds to option C, but the answer was {llm_choice}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_exoplanet_temperature_ratio()
print(result)