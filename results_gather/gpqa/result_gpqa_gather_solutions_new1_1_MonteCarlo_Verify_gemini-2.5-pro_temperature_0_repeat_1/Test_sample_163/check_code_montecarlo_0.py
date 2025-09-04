import math

def check_binary_star_mass_ratio():
    """
    This function checks the correctness of the provided answer for the binary star problem.
    It calculates the mass ratio based on the given physical parameters and compares it
    to the provided options to verify the selected answer.
    """

    # --- Problem Constraints and Data ---
    # System 1
    P1 = 2.0  # Period in years
    K1_1 = 10.0  # Radial velocity amplitude in km/s
    K1_2 = 5.0   # Radial velocity amplitude in km/s

    # System 2
    P2 = 1.0  # Period in years
    K2_1 = 15.0  # Radial velocity amplitude in km/s
    K2_2 = 10.0  # Radial velocity amplitude in km/s

    # The final answer provided by the LLM to be checked.
    # The last analysis block concludes with <<<B>>>.
    llm_answer_letter = "B"

    # The options as presented in the question.
    # Note: The final analysis block re-letters the options, but we will use the original lettering.
    # Original Question Options: A) ~ 0.6, B) ~ 0.4, C) ~ 0.7, D) ~ 1.2
    options = {
        "A": 0.6,
        "B": 0.4,
        "C": 0.7,
        "D": 1.2
    }

    # --- Physics Calculation ---
    # For a double-lined eclipsing binary, the inclination i ≈ 90°, so sin(i) ≈ 1.
    # The total mass M_total is proportional to P * (K1 + K2)³.
    # M_total ∝ P * (K_sum)³
    # The ratio of masses M_sys1 / M_sys2 is what we need to calculate.
    # The proportionality constant cancels out in the ratio.

    # Sum of radial velocity amplitudes for each system
    K1_sum = K1_1 + K1_2
    K2_sum = K2_1 + K2_2

    # Calculate the mass ratio
    try:
        calculated_ratio = (P1 * (K1_sum)**3) / (P2 * (K2_sum)**3)
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. P2 or K2_sum cannot be zero."

    # --- Verification ---
    # Find the option that is numerically closest to the calculated ratio.
    closest_option_letter = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_letter = letter

    # Check if the LLM's answer matches the calculated correct option.
    if llm_answer_letter == closest_option_letter:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect. "
            f"The calculation based on the provided data gives a mass ratio of M1/M2 = {calculated_ratio:.3f}. "
            f"The given options are A) ~0.6, B) ~0.4, C) ~0.7, D) ~1.2. "
            f"The calculated value {calculated_ratio:.3f} is closest to {options[closest_option_letter]} (Option {closest_option_letter}). "
            f"The provided answer was '{llm_answer_letter}', which does not match the correct option."
        )
        return reason

# Execute the check and print the result.
result = check_binary_star_mass_ratio()
print(result)