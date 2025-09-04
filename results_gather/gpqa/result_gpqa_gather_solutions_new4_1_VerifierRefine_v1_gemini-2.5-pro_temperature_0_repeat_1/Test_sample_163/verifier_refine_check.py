import math

def check_correctness():
    """
    Checks the correctness of the provided answer for the binary star mass ratio problem.

    The total mass of a double-lined spectroscopic, eclipsing binary system is proportional to:
    M ∝ P * (K1 + K2)³
    where P is the orbital period and K1, K2 are the radial velocity amplitudes.

    The ratio of the mass of system 1 to system 2 is therefore:
    M_sys1 / M_sys2 = [P1 * (K1a + K1b)³] / [P2 * (K2a + K2b)³]
    """

    # --- Define problem constraints and given values ---
    # System 1
    P1 = 2  # years
    K1a = 10 # km/s
    K1b = 5  # km/s

    # System 2
    P2 = 1  # year
    K2a = 15 # km/s
    K2b = 10 # km/s

    # The final answer from the LLM to be checked
    llm_answer_option = 'D'

    # --- Perform the calculation ---
    try:
        # Sum of velocity amplitudes for each system
        K_total1 = K1a + K1b
        K_total2 = K2a + K2b

        # Calculate the mass ratio
        calculated_ratio = (P1 * (K_total1)**3) / (P2 * (K_total2)**3)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verify the answer ---
    # The expected numerical result is 2 * (15/25)^3 = 2 * (3/5)^3 = 2 * 27/125 = 54/125 = 0.432
    expected_ratio = 0.432

    # Check if the calculation is correct
    if not math.isclose(calculated_ratio, expected_ratio, rel_tol=1e-6):
        return f"The calculation is incorrect. The calculated ratio is {calculated_ratio}, but the expected value is {expected_ratio}."

    # Define the options from the question
    options = {
        'A': 0.6,
        'B': 1.2,
        'C': 0.7,
        'D': 0.4
    }

    # Find which option is numerically closest to the calculated ratio
    closest_option = min(options.keys(), key=lambda option: abs(options[option] - calculated_ratio))

    # Check if the LLM's chosen option matches the closest option
    if llm_answer_option == closest_option:
        return "Correct"
    else:
        return (f"The final answer is incorrect. The calculated mass ratio is {calculated_ratio:.3f}. "
                f"This value is closest to option {closest_option} ({options[closest_option]}), "
                f"but the provided answer was option {llm_answer_option} ({options[llm_answer_option]}).")

# The final output of the checker code will be the return value of this function.
# print(check_correctness())