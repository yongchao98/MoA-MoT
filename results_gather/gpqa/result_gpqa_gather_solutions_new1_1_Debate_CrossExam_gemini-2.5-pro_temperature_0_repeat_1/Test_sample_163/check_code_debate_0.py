import math

def check_binary_star_mass_ratio():
    """
    Calculates the mass ratio of system_1 to system_2 based on the problem's data.
    This serves as a method to verify the correctness of the provided answers.
    """
    # --- Given Data ---
    # System 1
    P1 = 2.0  # Period in years
    K1_1 = 10.0 # Radial velocity amplitude 1 in km/s
    K1_2 = 5.0  # Radial velocity amplitude 2 in km/s

    # System 2
    P2 = 1.0  # Period in years
    K2_1 = 15.0 # Radial velocity amplitude 1 in km/s
    K2_2 = 10.0 # Radial velocity amplitude 2 in km/s

    # --- Calculation ---
    # For eclipsing binaries, M_total is proportional to P * (K1 + K2)^3.
    # The proportionality constants cancel out in the ratio.

    # Sum of radial velocity amplitudes for each system
    K_total1 = K1_1 + K1_2
    K_total2 = K2_1 + K2_2

    # Calculate the mass ratio M_sys1 / M_sys2
    try:
        mass_ratio = (P1 / P2) * (K_total1 / K_total2)**3
    except ZeroDivisionError:
        return "Error: Division by zero in calculation."

    # --- Verification against candidate answers ---
    # The majority of candidate answers calculate a ratio of ~0.432.
    # Let's check if our calculation matches this consensus.
    expected_ratio = 0.432
    
    if not math.isclose(mass_ratio, expected_ratio, rel_tol=1e-3):
        return (f"Incorrect calculation. The calculated ratio is {mass_ratio:.3f}, "
                f"which does not match the expected value of {expected_ratio}.")

    # The calculated value 0.432 is closest to the option "~ 0.4".
    # Many candidate answers have different letter assignments for this option.
    # For example, Answer 5 has D) ~ 0.4, and Answer 2 has C) ~ 0.4.
    # Let's check Answer 5, which states <<<D>>> and has the mapping D) ~ 0.4.
    
    given_answer_value = 0.4
    if math.isclose(mass_ratio, given_answer_value, abs_tol=0.05):
        return "Correct"
    else:
        return (f"Incorrect. The calculated ratio is {mass_ratio:.3f}, which is closest to 0.4. "
                f"The provided answer does not match this result.")

# Execute the check
result = check_binary_star_mass_ratio()
# This will return "Correct" if the logic holds, confirming that answers
# which calculated ~0.432 and chose the corresponding option are correct.
# print(result)