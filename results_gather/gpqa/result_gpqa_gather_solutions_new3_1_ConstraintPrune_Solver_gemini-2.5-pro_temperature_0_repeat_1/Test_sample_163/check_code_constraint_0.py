import math

def check_answer():
    """
    This function checks the correctness of the provided answer for the binary star system problem.
    """
    # --- 1. Define the given parameters from the question ---
    # System 1
    P1 = 2.0  # years
    K1_a = 10.0  # km/s
    K1_b = 5.0  # km/s

    # System 2
    P2 = 1.0  # year
    K2_a = 15.0  # km/s
    K2_b = 10.0  # km/s

    # --- 2. Perform the calculation based on physics principles ---
    # The total mass of an eclipsing binary system is proportional to P * (K_total)^3.
    # M_total ∝ P * (K_a + K_b)^3
    # We need to find the ratio M_sys1 / M_sys2.

    # Calculate the sum of radial velocity amplitudes for each system
    K_total1 = K1_a + K1_b
    K_total2 = K2_a + K2_b

    # Check for division by zero, although not possible with the given numbers
    if P2 == 0 or K_total2 == 0:
        return "Calculation error: Division by zero."

    # Calculate the mass ratio
    calculated_ratio = (P1 / P2) * (K_total1 / K_total2)**3

    # --- 3. Define the options and the provided answer ---
    # The final answer provided by the analyst uses the following options:
    # A) ~ 0.7, B) ~ 0.4, C) ~ 1.2, D) ~ 0.6
    options = {
        'A': 0.7,
        'B': 0.4,
        'C': 1.2,
        'D': 0.6
    }
    
    # The final answer given is <<<B>>>
    given_answer_letter = 'B'

    # --- 4. Verify the correctness ---
    # First, check if the calculation in the provided answer text (0.432) is correct.
    expected_calculation = 54.0 / 125.0
    if not math.isclose(calculated_ratio, expected_calculation, rel_tol=1e-6):
        return f"The calculation is incorrect. The ratio M_sys1 / M_sys2 should be (2/1) * (15/25)^3 = {expected_calculation:.3f}, but the code calculated {calculated_ratio:.3f}."

    # Second, find which option is numerically closest to the calculated ratio.
    closest_option_letter = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_letter = letter

    # Third, check if the chosen option in the answer matches the closest option.
    if closest_option_letter == given_answer_letter:
        return "Correct"
    else:
        return (f"The final answer is incorrect. The calculated mass ratio is {calculated_ratio:.3f}. "
                f"This value is closest to option {closest_option_letter} ({options[closest_option_letter]}), "
                f"but the provided answer was {given_answer_letter}.")

# Run the check
result = check_answer()
print(result)