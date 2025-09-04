import math

def check_planet_temperature_ratio():
    """
    This function checks the correctness of the answer to the exoplanet problem.
    It calculates the ratio of equilibrium temperatures for two planets and compares
    it to the provided answer choice.
    """

    # --- Given values from the question ---
    # Planet 1 properties
    M_p1_earth_mass = 7.0
    delta_lambda1_A = 0.03  # Doppler shift in Angstroms

    # Planet 2 properties
    M_p2_earth_mass = 5.0
    delta_lambda2_A = 0.04  # Doppler shift in Angstroms

    # The expected value from answer B
    expected_value = 0.53

    # --- Physics Derivation Summary ---
    # 1. The equilibrium temperature (T_eq) of a planet is given by:
    #    T_eq = T_star * sqrt(R_star / (2 * a)) * (1 - A)^(1/4)
    #    where 'a' is the semi-major axis and 'A' is the albedo.
    #
    # 2. The ratio T_eq1 / T_eq2 simplifies because the star's properties (T_star, R_star)
    #    and the planets' albedo (A) are the same for both.
    #    T_eq1 / T_eq2 = sqrt(a2 / a1)
    #
    # 3. The maximum radial velocity (K) of the star induced by a planet in a circular orbit is
    #    approximately K ≈ M_p * sqrt(G / (M_star * a)) for M_p << M_star.
    #    Solving for 'a' gives: a ∝ (M_p / K)^2.
    #
    # 4. The ratio of semi-major axes (a2 / a1) is therefore:
    #    a2 / a1 = [(M_p2 / K2)^2] / [(M_p1 / K1)^2] = (M_p2 / M_p1)^2 * (K1 / K2)^2
    #
    # 5. Substituting this back into the temperature ratio formula:
    #    T_eq1 / T_eq2 = sqrt((M_p2 / M_p1)^2 * (K1 / K2)^2) = (M_p2 / M_p1) * (K1 / K2)
    #
    # 6. The radial velocity K is directly proportional to the Doppler shift Δλ (K = c * Δλ / λ₀).
    #    Therefore, the ratio K1 / K2 is equal to the ratio Δλ1 / Δλ2.
    #    Final formula: T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    # --- Calculation ---
    # Ratio of planet masses (M_p2 / M_p1)
    mass_ratio = M_p2_earth_mass / M_p1_earth_mass

    # Ratio of Doppler shifts (Δλ1 / Δλ2)
    shift_ratio = delta_lambda1_A / delta_lambda2_A

    # Calculate the final equilibrium temperature ratio
    calculated_temp_ratio = mass_ratio * shift_ratio

    # --- Verification ---
    # Compare the calculated result with the expected value from the answer.
    # A tolerance is used for floating-point comparison.
    tolerance = 0.01
    if abs(calculated_temp_ratio - expected_value) < tolerance:
        # The calculated value matches the expected answer.
        print("Correct")
    else:
        # The calculated value does not match the expected answer.
        error_message = (
            f"Incorrect.\n"
            f"The calculated ratio of equilibrium temperatures is {calculated_temp_ratio:.4f}.\n"
            f"The expected value from the answer is ~{expected_value}.\n"
            f"The calculation is based on the formula: T_eq1/T_eq2 = (M_p2/M_p1) * (Δλ1/Δλ2).\n"
            f"  - Mass Ratio (M_p2/M_p1) = {M_p2_earth_mass}/{M_p1_earth_mass} = {mass_ratio:.4f}\n"
            f"  - Doppler Shift Ratio (Δλ1/Δλ2) = {delta_lambda1_A}/{delta_lambda2_A} = {shift_ratio:.4f}\n"
            f"  - Result: {mass_ratio:.4f} * {shift_ratio:.4f} = {calculated_temp_ratio:.4f}\n"
            f"The provided answer value of {expected_value} does not match the calculated value of {calculated_temp_ratio:.4f} within the tolerance of {tolerance}."
        )
        print(error_message)

# Execute the check
check_planet_temperature_ratio()