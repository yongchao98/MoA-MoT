import math

def check_binary_star_mass_ratio():
    """
    This function checks the correctness of the provided answer for the binary star mass ratio problem.
    """
    # --- Given Data ---
    # System 1
    P1 = 2.0  # Period in years
    K1a = 10.0 # Radial velocity amplitude in km/s
    K1b = 5.0  # Radial velocity amplitude in km/s

    # System 2
    P2 = 1.0  # Period in years
    K2a = 15.0 # Radial velocity amplitude in km/s
    K2b = 10.0 # Radial velocity amplitude in km/s

    # --- LLM's Answer and Options ---
    llm_answer_choice = 'B'
    options = {
        'A': 0.6,
        'B': 0.4,
        'C': 1.2,
        'D': 0.7
    }

    # --- Physics Calculation ---
    # The total mass M of an eclipsing, double-lined spectroscopic binary is proportional to P * (Ka + Kb)^3.
    # The constants of proportionality (like G and sin(i)) cancel out in the ratio M1/M2.
    # The units also cancel out, so no conversion is needed.

    # Calculate the sum of velocity amplitudes for each system
    K1_total = K1a + K1b
    K2_total = K2a + K2b

    # Calculate the mass ratio M1 / M2
    try:
        mass_ratio = (P1 * (K1_total**3)) / (P2 * (K2_total**3))
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. P2 or K2_total cannot be zero."

    # --- Verification ---
    # Find which option is numerically closest to the calculated ratio.
    # We use a lambda function to calculate the absolute difference for each option.
    closest_option = min(options.keys(), key=lambda option: abs(options[option] - mass_ratio))

    # Check if the LLM's answer matches the closest calculated option.
    if llm_answer_choice == closest_option:
        return "Correct"
    else:
        # If not correct, provide a detailed explanation.
        reason = (
            f"The provided answer is '{llm_answer_choice}', which corresponds to a value of ~{options[llm_answer_choice]}.\n"
            f"However, the calculated mass ratio is {mass_ratio:.3f}.\n"
            f"The calculation is based on the formula M ∝ P * (Ka + Kb)^3 for eclipsing binaries.\n"
            f"Calculation steps:\n"
            f"1. For System 1: P1 = {P1}, K1_total = {K1a} + {K1b} = {K1_total} km/s.\n"
            f"2. For System 2: P2 = {P2}, K2_total = {K2a} + {K2b} = {K2_total} km/s.\n"
            f"3. Mass Ratio M1/M2 = (P1 * K1_total^3) / (P2 * K2_total^3)\n"
            f"4. M1/M2 = ({P1} * {K1_total}^3) / ({P2} * {K2_total}^3) = {mass_ratio:.3f}.\n"
            f"The calculated value {mass_ratio:.3f} is closest to option '{closest_option}' (~{options[closest_option]}), not '{llm_answer_choice}'."
        )
        return reason

# Execute the check and print the result
result = check_binary_star_mass_ratio()
print(result)