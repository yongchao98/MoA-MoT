import math

def check_binary_star_mass_ratio():
    """
    This function checks the correctness of the provided answer for the binary star mass ratio problem.
    
    The problem states:
    - For a double-lined, eclipsing spectroscopic binary, the total mass (M_total) is proportional to the orbital period (P) multiplied by the cube of the sum of the radial velocity amplitudes (K_sum).
    - M_total ∝ P * (K_sum)³
    - The ratio of masses is M1 / M2 = (P1 * (K_sum1)³) / (P2 * (K_sum2)³)

    Data:
    - System 1: P1 = 2 years, K_sum1 = 10 + 5 = 15 km/s
    - System 2: P2 = 1 year, K_sum2 = 15 + 10 = 25 km/s

    Options:
    A) ~ 0.7
    B) ~ 0.4
    C) ~ 0.6
    D) ~ 1.2

    The provided answer is <<<B>>>.
    """
    
    # --- Define problem parameters ---
    # System 1
    P1 = 2.0  # years
    K1_sum = 10.0 + 5.0  # km/s

    # System 2
    P2 = 1.0  # year
    K2_sum = 15.0 + 10.0  # km/s

    # --- Perform the calculation ---
    try:
        # Calculate the mass ratio M1 / M2
        mass_ratio = (P1 * (K1_sum**3)) / (P2 * (K2_sum**3))
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. P2 or K2_sum cannot be zero."
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verify the answer ---
    # The expected numerical result is 0.432
    expected_ratio = 0.432
    if not math.isclose(mass_ratio, expected_ratio, rel_tol=1e-5):
        return f"Incorrect calculation. The calculated mass ratio is {mass_ratio:.5f}, but it should be {expected_ratio}."

    # The options provided in the question
    options = {
        'A': 0.7,
        'B': 0.4,
        'C': 0.6,
        'D': 1.2
    }
    
    # The final answer given by the LLM
    llm_answer_letter = 'B'
    
    # Find which option is numerically closest to the calculated ratio
    closest_option_letter = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(mass_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_letter = letter
            
    # Check if the LLM's chosen letter matches the closest numerical option
    if llm_answer_letter == closest_option_letter:
        return "Correct"
    else:
        return (f"Incorrect final answer. The calculated ratio is {mass_ratio:.3f}. "
                f"The closest option is '{closest_option_letter}' ({options[closest_option_letter]}), "
                f"but the provided answer was '{llm_answer_letter}' ({options[llm_answer_letter]}).")

# Execute the check and print the result
result = check_binary_star_mass_ratio()
print(result)