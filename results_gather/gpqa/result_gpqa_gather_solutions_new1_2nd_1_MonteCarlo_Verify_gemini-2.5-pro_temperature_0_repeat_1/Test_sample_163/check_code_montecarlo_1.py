import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the binary star system problem.
    """
    # --- Data from the question ---
    # System 1
    P1 = 2  # years
    K1a = 10  # km/s
    K1b = 5   # km/s

    # System 2
    P2 = 1  # year
    K2a = 15  # km/s
    K2b = 10  # km/s

    # --- The final answer to check ---
    # The provided final answer is <<<D>>>
    final_answer_letter = "D"

    # --- Options from the question ---
    # A) ~ 0.7, B) ~ 0.6, C) ~ 1.2, D) ~ 0.4
    options = {
        "A": 0.7,
        "B": 0.6,
        "C": 1.2,
        "D": 0.4
    }

    # --- Physics Calculation ---
    # For an eclipsing, double-lined spectroscopic binary, the total mass (M) is
    # proportional to the period (P) times the cube of the sum of the radial
    # velocity amplitudes (K_sum).
    # M ∝ P * (K_sum)^3
    # The question asks for the ratio M1 / M2.

    # Step 1: Calculate the sum of radial velocity amplitudes for each system.
    K_sum1 = K1a + K1b
    K_sum2 = K2a + K2b

    # Step 2: Calculate the mass ratio M1 / M2.
    # The constants of proportionality cancel out.
    try:
        calculated_ratio = (P1 * (K_sum1**3)) / (P2 * (K_sum2**3))
    except ZeroDivisionError:
        return "Calculation Error: Division by zero. P2 or K_sum2 cannot be zero."

    # --- Verification ---
    # Step 3: Find which option is numerically closest to the calculated ratio.
    closest_option_letter = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_letter = letter
            
    # Step 4: Check if the provided final answer matches the closest option.
    if final_answer_letter == closest_option_letter:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect. "
            f"The calculation for the mass ratio is M1/M2 = (P1 * (K1a+K1b)^3) / (P2 * (K2a+K2b)^3). "
            f"Plugging in the values: ({P1} * ({K_sum1})^3) / ({P2} * ({K_sum2})^3) = {calculated_ratio:.3f}. "
            f"The provided options are A) ~0.7, B) ~0.6, C) ~1.2, D) ~0.4. "
            f"The calculated value {calculated_ratio:.3f} is closest to {options[closest_option_letter]} (Option {closest_option_letter}). "
            f"The submitted answer was {final_answer_letter}, which corresponds to a value of ~{options[final_answer_letter]}, which is not the closest option."
        )
        return reason

# The final output of the check.
print(check_answer())