import math

def check_binary_star_mass_ratio():
    """
    This function verifies the calculation for the mass ratio of two binary star systems.
    It uses the provided data and the physical relationship between mass, period, and radial velocity.
    """
    
    # --- Data from the question ---
    # System 1
    P1 = 2.0  # Period in years
    K1a = 10.0 # Radial velocity amplitude in km/s
    K1b = 5.0  # Radial velocity amplitude in km/s

    # System 2
    P2 = 1.0  # Period in years
    K2a = 15.0 # Radial velocity amplitude in km/s
    K2b = 10.0 # Radial velocity amplitude in km/s

    # --- Physics Principle ---
    # For a double-lined spectroscopic binary, the total mass (M_total) is given by:
    # M_total = (P * (Ka + Kb)^3) / (2 * pi * G * sin^3(i))
    # Since we are calculating a ratio of two masses (M1 / M2), the constant (2*pi*G) cancels out.
    # The problem states both systems are "eclipsing", which implies the inclination angle i is ~90 degrees.
    # Therefore, sin(i) is ~1, and sin^3(i) is also ~1.
    # The simplified relationship for the ratio is:
    # Ratio = (P1 * (K1a + K1b)^3) / (P2 * (K2a + K2b)^3)

    # --- Calculation as per the provided answer's logic ---
    
    # Sum of velocity amplitudes for each system
    K_total_1 = K1a + K1b
    K_total_2 = K2a + K2b
    
    # Check if the intermediate sums are correct
    if K_total_1 != 15.0:
        return f"Incorrect sum of velocities for system 1. Expected 15.0, but calculated {K_total_1}."
    if K_total_2 != 25.0:
        return f"Incorrect sum of velocities for system 2. Expected 25.0, but calculated {K_total_2}."

    # Calculate the mass ratio
    try:
        calculated_ratio = (P1 * (K_total_1)**3) / (P2 * (K_total_2)**3)
    except ZeroDivisionError:
        return "Error: Division by zero during calculation. Check input values for P2 and K_total_2."

    # The answer's derivation gives 54/125 = 0.432
    expected_ratio = 54 / 125

    # Verify the final numerical result
    if not math.isclose(calculated_ratio, expected_ratio, rel_tol=1e-9):
        return f"Calculation mismatch. The provided answer calculates the ratio to be {expected_ratio:.3f}, but the code calculates it as {calculated_ratio:.3f}."

    # --- Check against the multiple-choice options ---
    # The provided answer is 'B', which corresponds to a value of ~0.4.
    options = {'A': 1.2, 'B': 0.4, 'C': 0.6, 'D': 0.7}
    
    # Find the option that is numerically closest to our calculated ratio
    closest_option_label = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # The final answer selected was 'B'
    selected_answer_label = 'B'

    if closest_option_label == selected_answer_label:
        return "Correct"
    else:
        return (f"Incorrect choice selection. The calculated ratio is {calculated_ratio:.3f}. "
                f"The closest option is '{closest_option_label}' (value ~{options[closest_option_label]}), "
                f"but the provided answer selected '{selected_answer_label}' (value ~{options[selected_answer_label]}).")

# Execute the check and print the result
result = check_binary_star_mass_ratio()
print(result)