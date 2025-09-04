import math

def check_answer():
    """
    This function checks the correctness of the provided answer for the binary star mass ratio problem.
    """
    # --- 1. Define the given parameters from the question ---
    # System 1
    P1 = 2.0  # Period in years
    K1_a = 10.0  # Radial velocity amplitude in km/s
    K1_b = 5.0   # Radial velocity amplitude in km/s

    # System 2
    P2 = 1.0  # Period in years
    K2_a = 15.0 # Radial velocity amplitude in km/s
    K2_b = 10.0  # Radial velocity amplitude in km/s

    # --- 2. Apply the relevant physical principle ---
    # For an eclipsing, double-lined spectroscopic binary, the total mass (M_total) is
    # proportional to the period (P) times the cube of the sum of the radial velocity amplitudes (K_sum).
    # M_total ∝ P * (K_sum)^3
    # The ratio of the masses M1/M2 is therefore:
    # (P1 * (K1_a + K1_b)^3) / (P2 * (K2_a + K2_b)^3)

    # Calculate the sum of radial velocities for each system
    K_sum1 = K1_a + K1_b
    K_sum2 = K2_a + K2_b

    # Calculate the mass ratio
    try:
        calculated_ratio = (P1 * (K_sum1**3)) / (P2 * (K_sum2**3))
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. P2 or K_sum2 cannot be zero."

    # --- 3. Define the options and the provided answer ---
    # The options as listed in the final analysis of the provided answer text:
    # A) ~ 0.6, B) ~ 0.4, C) ~ 0.7, D) ~ 1.2
    options = {
        'A': 0.6,
        'B': 0.4,
        'C': 0.7,
        'D': 1.2
    }
    
    # The final answer provided is <<<B>>>
    given_answer_letter = 'B'

    # --- 4. Verify the correctness ---
    # Find which option is numerically closest to the calculated ratio.
    closest_option = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option = letter
    
    # Check if the closest option matches the given answer.
    if closest_option == given_answer_letter:
        return "Correct"
    else:
        # If not, provide a reason for the discrepancy.
        reason = (
            f"The calculated mass ratio is {calculated_ratio:.3f}. "
            f"This value is closest to option {closest_option} ({options[closest_option]}), with a difference of {min_difference:.3f}. "
            f"The provided answer was {given_answer_letter} ({options[given_answer_letter]}), which has a difference of {abs(calculated_ratio - options[given_answer_letter]):.3f}. "
            f"Therefore, the provided answer is not the best fit for the calculated result."
        )
        return reason

# Execute the check and print the result
print(check_answer())