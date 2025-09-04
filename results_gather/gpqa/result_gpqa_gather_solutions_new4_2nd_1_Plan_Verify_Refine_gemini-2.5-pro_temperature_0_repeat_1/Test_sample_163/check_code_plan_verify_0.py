import math

def check_answer():
    """
    Checks the correctness of the provided answer for the binary star mass ratio problem.
    """
    # Data from the question
    # System 1
    P1 = 2.0  # years
    K1a = 10.0  # km/s
    K1b = 5.0  # km/s

    # System 2
    P2 = 1.0  # year
    K2a = 15.0  # km/s
    K2b = 10.0  # km/s

    # The options as provided in the final answer's analysis
    options = {
        'A': 0.4,
        'B': 0.7,
        'C': 1.2,
        'D': 0.6
    }
    
    # The final answer provided by the LLM
    llm_answer_letter = "A"

    # The core physical principle is that for an eclipsing binary system,
    # the total mass M is proportional to P * (K1 + K2)^3.
    # The ratio M_sys1 / M_sys2 can be calculated directly.
    
    # Calculate the sum of radial velocity amplitudes for each system
    K_total1 = K1a + K1b
    K_total2 = K2a + K2b

    # Check if the sums are correct as per the reasoning
    if K_total1 != 15.0 or K_total2 != 25.0:
        return f"Incorrect sum of velocities. Calculated K_total1={K_total1} (expected 15) and K_total2={K_total2} (expected 25)."

    # Calculate the mass ratio
    try:
        calculated_ratio = (P1 * (K_total1)**3) / (P2 * (K_total2)**3)
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. P2 or K_total2 cannot be zero."

    # Expected ratio from the detailed analysis
    expected_ratio = 54.0 / 125.0  # 0.432

    # Check if the calculation is correct
    if not math.isclose(calculated_ratio, expected_ratio, rel_tol=1e-5):
        return f"Calculation is incorrect. The calculated ratio is {calculated_ratio:.4f}, but it should be {expected_ratio:.4f}."

    # Check if the selected option is the closest one to the calculated result
    if llm_answer_letter not in options:
        return f"The provided answer '{llm_answer_letter}' is not a valid option."

    # Find the option that is numerically closest to the calculated result
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_ratio))

    if llm_answer_letter == closest_option:
        return "Correct"
    else:
        return (f"The final answer is incorrect. "
                f"The calculated ratio is {calculated_ratio:.3f}. "
                f"This value is closest to option {closest_option} (~{options[closest_option]}). "
                f"However, the provided answer was {llm_answer_letter} (~{options[llm_answer_letter]}).")

# Run the check
result = check_answer()
print(result)