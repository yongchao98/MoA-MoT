import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the binary star system problem.

    The total mass of a double-lined spectroscopic, eclipsing binary system is proportional to:
    M_total ∝ P * (K1 + K2)³
    where P is the period and K1, K2 are the radial velocity amplitudes.

    The factor by which system_1 is more massive than system_2 is the ratio of their masses:
    Ratio = M_sys1 / M_sys2 = [P1 * (K1_total)³] / [P2 * (K2_total)³]
    """

    # --- Define the given parameters from the question ---
    # System 1
    P1 = 2.0  # years
    K1_a = 10.0  # km/s
    K1_b = 5.0  # km/s

    # System 2
    P2 = 1.0  # year
    K2_a = 15.0  # km/s
    K2_b = 10.0  # km/s

    # --- Perform the calculation based on the physical formula ---
    # Sum of radial velocity amplitudes for each system
    K_total1 = K1_a + K1_b
    K_total2 = K2_a + K2_b

    # Calculate the mass ratio
    try:
        mass_ratio = (P1 / P2) * (K_total1 / K_total2)**3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Check the result against the provided options ---
    # The options from the question are:
    # A) ~ 0.6, B) ~ 0.7, C) ~ 1.2, D) ~ 0.4
    options = {
        "A": 0.6,
        "B": 0.7,
        "C": 1.2,
        "D": 0.4
    }
    
    # The final answer provided by the LLM is 'D'
    llm_answer_letter = "D"

    # Find the option closest to the calculated result
    closest_option_letter = min(options.keys(), key=lambda k: abs(options[k] - mass_ratio))
    
    # --- Verify the correctness ---
    # 1. Check if the calculation is correct
    expected_value = 0.432
    if not math.isclose(mass_ratio, expected_value, rel_tol=1e-3):
        return f"Incorrect calculation. The calculated mass ratio is {mass_ratio:.3f}, but it should be approximately {expected_value}."

    # 2. Check if the chosen option corresponds to the calculation
    if closest_option_letter != "D":
        return f"The calculated value is {mass_ratio:.3f}, which is closest to option {closest_option_letter} ({options[closest_option_letter]}), not option D (0.4)."

    # 3. Check if the final answer matches the correct option
    if llm_answer_letter == closest_option_letter:
        return "Correct"
    else:
        return f"The final answer is incorrect. The calculated ratio is {mass_ratio:.3f}, which is closest to option {closest_option_letter} ({options[closest_option_letter]}). The provided answer was {llm_answer_letter}."

# Run the check
result = check_correctness()
print(result)