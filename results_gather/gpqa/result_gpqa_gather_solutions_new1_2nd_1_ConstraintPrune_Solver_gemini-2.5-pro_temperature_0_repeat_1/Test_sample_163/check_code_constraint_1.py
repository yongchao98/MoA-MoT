import re
import numpy as np

def check_binary_star_mass_ratio():
    """
    Checks the correctness of the LLM's answer for the binary star mass ratio problem.

    The problem requires calculating the ratio of the total mass of system_1 to system_2.
    For an eclipsing, double-lined spectroscopic binary, the total mass (M_total) is
    proportional to the orbital period (P) multiplied by the cube of the sum of the
    radial velocity amplitudes (K_sum).
    
    M_total ∝ P * (K_sum)^3

    Therefore, the ratio is:
    Ratio = (P1 * (K_sum1)^3) / (P2 * (K_sum2)^3)
    """
    
    # --- Data from the question ---
    # System 1
    P1 = 2.0  # years
    K1_sys1 = 10.0  # km/s
    K2_sys1 = 5.0   # km/s
    K_sum1 = K1_sys1 + K2_sys1

    # System 2
    P2 = 1.0  # year
    K1_sys2 = 15.0  # km/s
    K2_sys2 = 10.0  # km/s
    K_sum2 = K1_sys2 + K2_sys2

    # --- Theoretical Calculation ---
    # Calculate the ratio of the masses M1 / M2
    try:
        calculated_ratio = (P1 * (K_sum1)**3) / (P2 * (K_sum2)**3)
    except ZeroDivisionError:
        return "Error in calculation: Division by zero."

    # --- Options from the question prompt ---
    # A) ~ 0.4, B) ~ 1.2, C) ~ 0.6, D) ~ 0.7
    options = {
        'A': 0.4,
        'B': 1.2,
        'C': 0.6,
        'D': 0.7
    }

    # --- Extract the LLM's final answer ---
    llm_answer_text = """
The analysis and code verification have confirmed the correct answer.

**Step-by-Step Breakdown:**

1.  **Physics Principle:** The total mass of an eclipsing, double-lined spectroscopic binary is proportional to its orbital period (`P`) multiplied by the cube of the sum of its radial velocity amplitudes (`K_sum`). The formula simplifies to `M_total ∝ P * (K_sum)³`.
2.  **Data Extraction:**
    *   System 1: `P1 = 2` years, `K_sum1 = 10 + 5 = 15` km/s.
    *   System 2: `P2 = 1` year, `K_sum2 = 15 + 10 = 25` km/s.
3.  **Calculation:** The ratio of the masses is:
    *   Ratio = `M1 / M2 = [P1 * (K_sum1)³] / [P2 * (K_sum2)³]`
    *   Ratio = `[2 * (15)³] / [1 * (25)³]`
    *   Ratio = `2 * (15/25)³`
    *   Ratio = `2 * (3/5)³`
    *   Ratio = `2 * (27/125)`
    *   Ratio = `54 / 125`
    *   Ratio = `0.432`
4.  **Conclusion:** The calculated ratio is `0.432`. This value is closest to the option `~ 0.4`. The code output confirms that `0.4` is the best fit among the given choices.

<<<A>>>
"""
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect. The final answer is not in the required format '<<<X>>>'."
    
    llm_choice_letter = match.group(1)

    # --- Verification ---
    # Find which option is numerically closest to the calculated ratio
    option_values = np.array(list(options.values()))
    differences = np.abs(option_values - calculated_ratio)
    closest_index = np.argmin(differences)
    correct_option_letter = list(options.keys())[closest_index]
    
    # Check if the LLM's reasoning and calculation are correct
    # The LLM calculated 0.432, which is correct.
    # 2 * (15/25)**3 = 2 * (0.6)**3 = 2 * 0.216 = 0.432
    if not np.isclose(calculated_ratio, 0.432):
        return f"Incorrect. The calculation is wrong. The expected ratio is 0.432, but the code calculated {calculated_ratio:.3f}."

    # Check if the LLM's final choice matches the correct choice
    if llm_choice_letter == correct_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The calculated mass ratio is {calculated_ratio:.3f}. "
                f"This value is closest to option {correct_option_letter} (value ~{options[correct_option_letter]}). "
                f"The provided answer was option {llm_choice_letter} (value ~{options[llm_choice_letter]}).")

# Execute the check and print the result
result = check_binary_star_mass_ratio()
print(result)