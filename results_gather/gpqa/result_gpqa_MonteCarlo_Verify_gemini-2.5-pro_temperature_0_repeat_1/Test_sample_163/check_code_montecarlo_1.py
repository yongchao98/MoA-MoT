import math

def check_binary_star_mass_ratio():
    """
    Checks the correctness of the provided answer for the binary star mass ratio problem.
    """
    # Given parameters for system 1
    P1 = 2.0  # years
    K1_1 = 10.0  # km/s
    K1_2 = 5.0  # km/s

    # Given parameters for system 2
    P2 = 1.0  # year
    K2_1 = 15.0  # km/s
    K2_2 = 10.0  # km/s

    # The provided answer corresponds to option D
    provided_answer_value = 1.2

    # --- Step 1: Calculate the sum of radial velocity amplitudes for each system ---
    K1_sum = K1_1 + K1_2
    K2_sum = K2_1 + K2_2

    # --- Step 2: Calculate the actual mass ratio based on the correct formula ---
    # The question asks for the mass ratio (M1 / M2).
    # The formula is M_ratio = (P1 / P2) * (K1_sum / K2_sum)^3
    
    P_ratio = P1 / P2
    K_sum_ratio = K1_sum / K2_sum
    
    correct_mass_ratio = P_ratio * (K_sum_ratio)**3

    # --- Step 3: Check if the provided answer matches the calculated correct answer ---
    # We check if the provided answer is close to the calculated one.
    if not math.isclose(correct_mass_ratio, provided_answer_value, rel_tol=0.1):
        # The answer is incorrect. Explain why.
        
        # As a sanity check, let's see what physical quantity would give the provided answer.
        # The ratio of semi-major axes 'a' is proportional to P * K_sum.
        semi_major_axis_ratio = P_ratio * K_sum_ratio
        
        reason = (
            f"The provided answer is incorrect. The question asks for the mass ratio (M1/M2).\n"
            f"The correct formula for the mass ratio of two binary systems is M1/M2 = (P1/P2) * ((K1_1+K1_2)/(K2_1+K2_2))^3.\n"
            f"Calculation:\n"
            f"P1 = {P1}, P2 = {P2}  => P1/P2 = {P_ratio}\n"
            f"K1_sum = {K1_1} + {K1_2} = {K1_sum} km/s\n"
            f"K2_sum = {K2_1} + {K2_2} = {K2_sum} km/s\n"
            f"K1_sum/K2_sum = {K1_sum}/{K2_sum} = {K_sum_ratio}\n"
            f"Mass Ratio M1/M2 = {P_ratio} * ({K_sum_ratio})^3 = {P_ratio} * {K_sum_ratio**3:.4f} = {correct_mass_ratio:.4f}\n\n"
            f"The calculated correct mass ratio is approximately {correct_mass_ratio:.3f}, which corresponds to option B (~0.4).\n"
            f"The provided answer is D (~1.2). This value is not the mass ratio.\n"
            f"The value 1.2 is obtained by calculating the ratio of the semi-major axes (a1/a2), which follows the formula (P1/P2) * (K1_sum/K2_sum) = {semi_major_axis_ratio:.1f}. "
            f"However, the question explicitly asks for the mass ratio, not the semi-major axis ratio."
        )
        return reason
    else:
        return "Correct"

# Run the checker
result = check_binary_star_mass_ratio()
print(result)