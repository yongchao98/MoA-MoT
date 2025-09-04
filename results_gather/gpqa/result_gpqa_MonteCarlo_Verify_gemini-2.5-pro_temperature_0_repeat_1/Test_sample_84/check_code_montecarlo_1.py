import math

def check_temperature_ratio():
    """
    This function checks the correctness of the LLM's answer by recalculating the ratio of equilibrium temperatures.
    """
    # --- Given values from the question ---
    # Planet 1
    Mp1 = 7.0  # in Earth masses
    delta_lambda1 = 0.03  # in Angstroms

    # Planet 2
    Mp2 = 5.0  # in Earth masses
    delta_lambda2 = 0.04  # in Angstroms

    # --- LLM's Answer ---
    llm_answer_option = 'C'
    llm_answer_value = 0.53

    # --- Derivation and Calculation ---
    # 1. The ratio of equilibrium temperatures T_eq1 / T_eq2 simplifies to sqrt(a2 / a1),
    #    where 'a' is the semi-major axis.
    # 2. The semi-major axis 'a' is proportional to (Mp / K)^2, where Mp is planet mass
    #    and K is the radial velocity semi-amplitude.
    # 3. The radial velocity K is proportional to the Doppler shift Δλ.
    # 4. Combining these, we get:
    #    T_eq1 / T_eq2 = sqrt( (Mp2/K2)^2 / (Mp1/K1)^2 )
    #    T_eq1 / T_eq2 = sqrt( (Mp2/Mp1)^2 * (K1/K2)^2 )
    #    T_eq1 / T_eq2 = (Mp2 / Mp1) * (K1 / K2)
    #    T_eq1 / T_eq2 = (Mp2 / Mp1) * (Δλ1 / Δλ2)

    # Calculate the ratios
    mass_ratio = Mp2 / Mp1
    doppler_shift_ratio = delta_lambda1 / delta_lambda2

    # Calculate the final temperature ratio
    calculated_temp_ratio = mass_ratio * doppler_shift_ratio

    # The exact fractional result is (5/7) * (0.03/0.04) = (5/7) * (3/4) = 15/28
    expected_exact_ratio = 15 / 28

    # --- Verification ---
    # Check if the calculated value matches the expected exact value
    if not math.isclose(calculated_temp_ratio, expected_exact_ratio):
        return f"Incorrect calculation logic. The calculated ratio is {calculated_temp_ratio:.4f}, but the physics implies a ratio of {expected_exact_ratio:.4f}."

    # Check if the calculated value is consistent with the chosen option C
    # We use a tolerance because the option is rounded. A tolerance of 0.01 is appropriate for an option with 2 decimal places.
    if math.isclose(calculated_temp_ratio, llm_answer_value, rel_tol=0.01):
        return "Correct"
    else:
        return f"The calculation is correct ({calculated_temp_ratio:.4f}), but the chosen option value ({llm_answer_value}) is not the closest or is poorly rounded. The calculated value rounds to {calculated_temp_ratio:.2f}, which matches option C."

# Run the check
result = check_temperature_ratio()
print(result)