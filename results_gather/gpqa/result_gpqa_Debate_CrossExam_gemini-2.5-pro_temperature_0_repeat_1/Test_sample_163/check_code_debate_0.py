import math

def check_binary_star_mass_ratio():
    """
    This function checks the correctness of the LLM's answer for the binary star mass ratio problem.
    It recalculates the ratio based on the provided formula and data and compares it to the LLM's result and chosen option.
    """
    # Given values from the question
    P1 = 2.0  # Period of system_1 in years
    K1a = 10.0  # RV amplitude of star 1a in km/s
    K1b = 5.0   # RV amplitude of star 1b in km/s
    
    P2 = 1.0  # Period of system_2 in years
    K2a = 15.0  # RV amplitude of star 2a in km/s
    K2b = 10.0  # RV amplitude of star 2b in km/s

    # The LLM's final answer and the corresponding option value
    llm_answer_option = 'D'
    llm_calculated_value = 0.432
    
    # The options provided in the question
    options = {
        'A': 0.6,
        'B': 0.7,
        'C': 1.2,
        'D': 0.4
    }

    # The formula for the ratio of total masses M1/M2 is derived from the mass function for a double-lined spectroscopic binary:
    # M_total * sin^3(i) = P * (K_a + K_b)^3 / (2 * pi * G)
    # Since both systems are eclipsing, we assume inclination i ≈ 90°, so sin(i) ≈ 1.
    # The term sin^3(i) and the constant (2 * pi * G) cancel out in the ratio.
    # Therefore, M1 / M2 = (P1 * (K1a + K1b)^3) / (P2 * (K2a + K2b)^3)

    # Step 1: Calculate the sum of the radial velocity amplitudes for each system.
    K_total1 = K1a + K1b
    K_total2 = K2a + K2b

    # Step 2: Substitute the values into the ratio formula.
    # The units (years and km/s) are consistent for the ratio and cancel out.
    try:
        calculated_ratio = (P1 * (K_total1**3)) / (P2 * (K_total2**3))
    except ZeroDivisionError:
        return "Error in calculation: Division by zero. P2 or K_total2 cannot be zero."

    # Step 3: Verify the LLM's calculation.
    # We check if our calculated ratio is close to the LLM's calculated value.
    if not math.isclose(calculated_ratio, llm_calculated_value, rel_tol=1e-3):
        return (f"The LLM's intermediate calculation is incorrect. "
                f"It calculated the ratio to be {llm_calculated_value}, but the correct value is {calculated_ratio:.4f}.")

    # Step 4: Verify if the chosen option is the most appropriate one.
    # Find which option is closest to our calculated value.
    closest_option = min(options.keys(), key=lambda option: abs(options[option] - calculated_ratio))

    if closest_option != llm_answer_option:
        return (f"The LLM chose the wrong option. "
                f"The calculated ratio is {calculated_ratio:.4f}, which is closest to option {closest_option} (~{options[closest_option]}), "
                f"but the LLM chose option {llm_answer_option} (~{options[llm_answer_option]}).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
print(check_binary_star_mass_ratio())