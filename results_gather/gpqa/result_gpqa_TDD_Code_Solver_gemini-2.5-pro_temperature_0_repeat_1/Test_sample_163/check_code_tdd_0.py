import math

def check_binary_star_mass_ratio():
    """
    This function checks the correctness of the provided answer to the binary star mass ratio problem.
    It recalculates the mass ratio based on the problem's data and compares it to the given options.
    """
    # --- 1. Define problem constraints and data from the question ---
    # System 1 parameters
    P1 = 2.0  # Period in years
    K1a = 10.0 # Radial velocity amplitude in km/s
    K1b = 5.0  # Radial velocity amplitude in km/s

    # System 2 parameters
    P2 = 1.0  # Period in years
    K2a = 15.0 # Radial velocity amplitude in km/s
    K2b = 10.0 # Radial velocity amplitude in km/s

    # The multiple-choice options provided in the question
    options = {
        'A': 0.6,
        'B': 0.7,
        'C': 1.2,
        'D': 0.4
    }

    # The answer provided by the LLM
    llm_answer_choice = 'D'

    # --- 2. Perform the calculation based on physics principles ---
    # The total mass (M) of a double-lined spectroscopic binary is given by:
    # M = (P * (K1 + K2)^3) / (2 * pi * G * sin^3(i))
    # Since both are eclipsing binaries, we assume inclination i is ~90 degrees, so sin(i) is ~1.
    # When taking the ratio M1/M2, the constant terms (2*pi*G*sin^3(i)) cancel out.
    # The formula for the ratio simplifies to:
    # M1 / M2 = (P1 * (K1a + K1b)^3) / (P2 * (K2a + K2b)^3)

    try:
        # Sum of radial velocity amplitudes for each system
        k1_sum = K1a + K1b
        k2_sum = K2a + K2b

        # Calculate the mass ratio
        calculated_ratio = (P1 * (k1_sum ** 3)) / (P2 * (k2_sum ** 3))
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- 3. Verify the correctness of the LLM's answer ---
    # Find which option is numerically closest to our calculated ratio
    closest_option = None
    min_difference = float('inf')

    for choice, value in options.items():
        difference = abs(calculated_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option = choice

    # Check if the LLM's chosen answer matches the closest option we found
    if llm_answer_choice == closest_option:
        # The answer is consistent with the calculation.
        # The calculation in the provided answer (0.432) is also correct.
        return "Correct"
    else:
        # The answer is incorrect. Provide a detailed reason.
        return (f"Incorrect. The calculated mass ratio is {calculated_ratio:.3f}. "
                f"This value is closest to option {closest_option} (which is ~{options[closest_option]}), "
                f"but the provided answer was option {llm_answer_choice} (which is ~{options[llm_answer_choice]}).")

# Execute the check and print the result
result = check_binary_star_mass_ratio()
print(result)