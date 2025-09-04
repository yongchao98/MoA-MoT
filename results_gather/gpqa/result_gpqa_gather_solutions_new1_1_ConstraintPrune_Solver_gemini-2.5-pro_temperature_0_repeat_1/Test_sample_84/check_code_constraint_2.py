import math

def check_planetary_temperature_ratio():
    """
    This function checks the correctness of the provided answer for the planetary temperature ratio problem.
    It calculates the correct ratio based on the given physical parameters and compares it to the proposed answer.
    """
    # --- Given Parameters ---
    # Planet 1
    m1 = 7.0  # Mass in Earth masses
    d_lambda1 = 0.03  # Doppler shift in Angstroms

    # Planet 2
    m2 = 5.0  # Mass in Earth masses
    d_lambda2 = 0.04  # Doppler shift in Angstroms

    # --- Physics Principles & Derivation ---
    # 1. The equilibrium temperature (T_eq) of a planet is proportional to the inverse square root of its semi-major axis (a).
    #    T_eq ∝ 1 / sqrt(a)
    #    Therefore, the ratio of temperatures is: T_eq1 / T_eq2 = sqrt(a2 / a1)

    # 2. The radial velocity semi-amplitude (K) of the star is proportional to the planet's mass (M_p) and inversely
    #    proportional to the square root of the semi-major axis (a).
    #    K ∝ M_p / sqrt(a)
    #    Rearranging for 'a', we get: a ∝ (M_p / K)^2

    # 3. The ratio of the semi-major axes (a2 / a1) can be expressed as:
    #    a2 / a1 = (m2 / k2)^2 / (m1 / k1)^2 = (m2 / m1)^2 * (k1 / k2)^2

    # 4. The radial velocity semi-amplitude (K) is directly proportional to the observed Doppler shift (Δλ).
    #    Therefore, the ratio K1 / K2 is equal to the ratio Δλ1 / Δλ2.

    # 5. Combining these relationships, we get the final formula for the temperature ratio:
    #    T_eq1 / T_eq2 = sqrt(a2 / a1) = sqrt( (m2 / m1)^2 * (d_lambda1 / d_lambda2)^2 )
    #    T_eq1 / T_eq2 = (m2 / m1) * (d_lambda1 / d_lambda2)

    # --- Calculation ---
    calculated_ratio = (m2 / m1) * (d_lambda1 / d_lambda2)

    # --- Verification ---
    # The provided answer to check is A, which corresponds to a value of ~1.05.
    # The options are: A) ~1.05, B) ~0.98, C) ~1.30, D) ~0.53
    expected_value_for_A = 1.05
    correct_value_for_D = 0.53

    # Check if the calculated ratio matches the value for option A.
    if abs(calculated_ratio - expected_value_for_A) < 0.02:
        return "Correct"
    else:
        # Explain why the answer is incorrect and provide the correct calculation.
        error_reason = (
            f"The provided answer A (~1.05) is incorrect.\n"
            f"The correct formula for the temperature ratio is T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2).\n"
            f"Plugging in the values: T_eq1 / T_eq2 = (5 / 7) * (0.03 / 0.04) = 15 / 28 ≈ {calculated_ratio:.4f}.\n"
            f"This calculated value of ~0.536 corresponds to option D (~0.53).\n"
            f"The incorrect answer A (~1.05) is likely the result of inverting the mass ratio in the calculation: (7 / 5) * (0.03 / 0.04) = 21 / 20 = 1.05. This is based on a physically incorrect derivation."
        )
        return error_reason

# Run the check and print the result.
result = check_planetary_temperature_ratio()
print(result)