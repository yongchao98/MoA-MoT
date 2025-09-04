import math

def check_binary_star_mass_ratio():
    """
    This function checks the correctness of the answer to the binary star mass ratio problem.
    It recalculates the ratio based on the provided data and compares it to the given options.
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

    # --- The final answer provided by the LLM to be checked ---
    llm_final_answer_letter = "D"

    # --- Step-by-step calculation based on physics principles ---
    # For eclipsing, double-lined spectroscopic binaries, the total mass M is proportional to P * (K_total)^3.
    # M_total ∝ P * (K₁ + K₂)^3
    # We need to find the ratio M_sys1 / M_sys2.

    # 1. Calculate the sum of radial velocity amplitudes for each system.
    K_total1 = K1a + K1b
    K_total2 = K2a + K2b

    # 2. Check if the intermediate sums are correct.
    if K_total1 != 15:
        return f"Constraint not satisfied: The sum of radial velocities for system 1 should be 10 + 5 = 15 km/s, but was calculated as {K_total1}."
    if K_total2 != 25:
        return f"Constraint not satisfied: The sum of radial velocities for system 2 should be 15 + 10 = 25 km/s, but was calculated as {K_total2}."

    # 3. Calculate the mass ratio.
    # The ratio M_sys1 / M_sys2 = [P₁ * (K_total1)³] / [P₂ * (K_total2)³]
    try:
        calculated_ratio = (P1 / P2) * (K_total1 / K_total2)**3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # 4. Compare the calculated ratio with the options to find the closest one.
    options = {
        "A": 0.7,
        "B": 0.6,
        "C": 1.2,
        "D": 0.4
    }

    # Find the option that is numerically closest to our calculated result.
    closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # 5. Verify the final answer.
    # The calculated value should be approximately 0.432.
    if not math.isclose(calculated_ratio, 0.432, rel_tol=1e-5):
        return f"Incorrect calculation. The expected ratio is ~0.432, but the calculated value is {calculated_ratio:.5f}."

    # The closest option to 0.432 is D (0.4).
    if closest_option_letter != "D":
        return f"Incorrect logic in selecting the final option. The calculated ratio is {calculated_ratio:.3f}, which is closest to option D (0.4), not {closest_option_letter}."

    # Check if the LLM's final answer matches the correct option.
    if llm_final_answer_letter == closest_option_letter:
        return "Correct"
    else:
        return f"Incorrect. The calculated ratio is {calculated_ratio:.3f}, which corresponds to option {closest_option_letter} (~{options[closest_option_letter]}). The provided answer was {llm_final_answer_letter}."

# Execute the check and print the result.
result = check_binary_star_mass_ratio()
print(result)