import math

def check_answer():
    """
    Checks the correctness of the final answer for the binary star system problem.
    """
    # Given data from the question
    # System 1
    P1 = 2  # years
    K1a = 10  # km/s
    K1b = 5   # km/s

    # System 2
    P2 = 1  # year
    K2a = 15  # km/s
    K2b = 10  # km/s

    # The final answer provided by the LLM to be checked
    final_answer_letter = "A"

    # The options provided in the question analysis
    options = {
        "A": 0.4,
        "B": 0.7,
        "C": 0.6,
        "D": 1.2
    }

    # The total mass of an eclipsing binary system is proportional to P * (K_total)^3
    # M_total ∝ P * (K_a + K_b)^3
    # We need to find the ratio M_sys1 / M_sys2

    # Calculate the sum of radial velocity amplitudes for each system
    K_total1 = K1a + K1b
    K_total2 = K2a + K2b

    # Check if the sums are calculated correctly in the provided analysis
    if K_total1 != 15:
        return f"Incorrect calculation for the sum of velocities in system 1. Expected 15 km/s, but calculation might be different."
    if K_total2 != 25:
        return f"Incorrect calculation for the sum of velocities in system 2. Expected 25 km/s, but calculation might be different."

    # Calculate the ratio of the masses
    try:
        mass_ratio = (P1 / P2) * (K_total1 / K_total2)**3
    except ZeroDivisionError:
        return "Error in calculation: Division by zero occurred. P2 or K_total2 cannot be zero."

    # The expected value is 2 * (15/25)^3 = 2 * (0.6)^3 = 2 * 0.216 = 0.432
    if not math.isclose(mass_ratio, 0.432, rel_tol=1e-5):
        return f"The calculated mass ratio is {mass_ratio:.5f}, which does not match the expected value of 0.432. The formula used might be incorrect."

    # Find which option is closest to the calculated value
    closest_option_letter = min(options.keys(), key=lambda k: abs(options[k] - mass_ratio))

    # Check if the selected answer letter matches the closest option
    if final_answer_letter == closest_option_letter:
        return "Correct"
    else:
        return (f"The final answer is incorrect. "
                f"The calculated mass ratio is {mass_ratio:.3f}. "
                f"The closest option value is {options[closest_option_letter]} which corresponds to option '{closest_option_letter}'. "
                f"The provided answer was '{final_answer_letter}'.")

# Run the check
result = check_answer()
print(result)