import math

def check_exoplanet_period_ratio():
    """
    This function checks the correctness of the provided answer by recalculating the ratio of orbital periods.
    """
    # Given values from the question
    delta_lambda_1 = 5.0  # miliangstrom for planet #1
    delta_lambda_2 = 7.0  # miliangstrom for planet #2

    # The provided answer is 'D', which corresponds to a ratio of approximately 0.36.
    # The logic in the provided solution leads to the calculation (5/7)^3.
    # We will verify this logic and calculation.

    # Constraint Check: The problem states M_star1 = M_star2 and M_p1 = M_p2.
    # The provided solution correctly uses this by treating them as constants.

    # Step 1: The amplitude of the wavelength shift (delta_lambda) is proportional to the
    # amplitude of the star's radial velocity (K).
    # K2 / K1 = delta_lambda_2 / delta_lambda_1
    ratio_K = delta_lambda_2 / delta_lambda_1
    
    if not math.isclose(ratio_K, 7.0/5.0):
        return f"Incorrect calculation of velocity ratio. Expected {7.0/5.0}, but got {ratio_K}."

    # Step 2: The radial velocity amplitude K is proportional to T^(-1/3) when
    # star mass and planet mass are constant.
    # K ∝ T^(-1/3)
    # This implies: K2 / K1 = (T1 / T2)^(1/3)
    # The logic used in the provided solution is correct.

    # Step 3: Solve for the ratio of the periods, T2 / T1.
    # (K2 / K1)^3 = T1 / T2
    # T2 / T1 = 1 / (K2 / K1)^3 = (K1 / K2)^3
    # Since K2/K1 = 7/5, then K1/K2 = 5/7.
    calculated_ratio_T2_T1 = (5.0 / 7.0)**3

    # Step 4: Compare the result with the provided answer 'D' (~0.36).
    # The value for option D is approximately 0.36.
    # Our calculated value is ~0.3644.
    # The provided solution correctly identifies 'D' as the answer.
    
    # Let's check if the calculated value is close to the value for option D.
    # The options are A) 0.85, B) 1.40, C) 1.96, D) 0.36
    options = {'A': 0.85, 'B': 1.40, 'C': 1.96, 'D': 0.36}
    
    # Find the option closest to our calculated value
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_ratio_T2_T1))

    if closest_option == 'D':
        # The provided answer 'D' is indeed the closest and correct choice.
        # The reasoning and calculation in the solution are sound.
        return "Correct"
    else:
        return f"The calculated ratio is {calculated_ratio_T2_T1:.4f}, which is closest to option {closest_option} ({options[closest_option]}), not D (0.36). The provided answer is incorrect."

# Execute the check
result = check_exoplanet_period_ratio()
print(result)