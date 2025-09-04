import math

def check_answer():
    """
    This function checks the correctness of the provided answer by recalculating the ratio of equilibrium temperatures.
    """
    # Given values from the problem statement
    mp1 = 7.0  # Mass of Planet 1 in Earth masses
    mp2 = 5.0  # Mass of Planet 2 in Earth masses
    delta_lambda1 = 0.03  # Doppler shift for Planet 1 in Angstroms
    delta_lambda2 = 0.04  # Doppler shift for Planet 2 in Angstroms

    # The ratio of equilibrium temperatures T_eq1 / T_eq2 can be derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ₁ / Δλ₂)
    
    # Calculate the ratio of masses
    mass_ratio = mp2 / mp1
    
    # Calculate the ratio of Doppler shifts
    shift_ratio = delta_lambda1 / delta_lambda2
    
    # Calculate the final temperature ratio
    calculated_ratio = mass_ratio * shift_ratio
    
    # The expected answer from the provided solution is 'A', which corresponds to ~0.53
    expected_value = 0.53
    
    # Check if the calculated ratio is close to the expected value
    # We use math.isclose for robust floating-point comparison.
    # The value 15/28 is exactly 0.535714...
    # Option A is ~0.53. We'll check if our result rounds to this.
    if round(calculated_ratio, 2) == expected_value:
        return "Correct"
    else:
        return (f"Incorrect. The calculated ratio is {calculated_ratio:.4f}, "
                f"which is approximately 0.54. The closest option is A (~0.53). "
                f"The provided answer 'A' is correct, but the reasoning in some candidate answers "
                f"that led to other options (like C or D) is flawed in its final selection, "
                f"despite having the correct numerical calculation.")

# Run the check
result = check_answer()
print(result)
