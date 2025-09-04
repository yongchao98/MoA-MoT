import math

def check_answer():
    """
    This function verifies the correctness of the provided answer for the binary star system problem.
    """
    # --- Define problem constraints and given values ---
    # System 1
    P1 = 2.0  # Period in years
    K1_1 = 10.0 # Radial velocity amplitude in km/s
    K1_2 = 5.0  # Radial velocity amplitude in km/s

    # System 2
    P2 = 1.0  # Period in years
    K2_1 = 15.0 # Radial velocity amplitude in km/s
    K2_2 = 10.0 # Radial velocity amplitude in km/s

    # The provided answer from the LLM
    llm_answer_option = 'B'
    options = {'A': 0.6, 'B': 0.4, 'C': 1.2, 'D': 0.7}
    
    # --- Physical Principle ---
    # The total mass (M = m1 + m2) of a double-lined spectroscopic binary is given by:
    # M = [P * (K1 + K2)^3] / [2 * pi * G * sin^3(i)]
    # The problem states both systems are "eclipsing", which implies the orbital inclination i is approximately 90 degrees.
    # Therefore, sin(i) ≈ 1, and sin^3(i) ≈ 1.
    # This simplifies the mass calculation and allows for a direct ratio comparison.

    # --- Calculation ---
    # We need to find the ratio M_system1 / M_system2.
    # The constant term (2 * pi * G) cancels out.
    # M1 / M2 = (P1 / P2) * [ (K1_1 + K1_2) / (K2_1 + K2_2) ]^3

    # Sum of radial velocity amplitudes for each system
    K_sum_1 = K1_1 + K1_2
    K_sum_2 = K2_1 + K2_2

    # Calculate the ratio
    try:
        calculated_ratio = (P1 / P2) * (K_sum_1 / K_sum_2)**3
    except ZeroDivisionError:
        return "Incorrect. The calculation involves a division by zero based on the input values."

    # --- Verification ---
    # Find the option that is numerically closest to our calculated result.
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_ratio))

    # Check if the LLM's chosen option is the closest one.
    if llm_answer_option == closest_option:
        # Additionally, check if the calculated value is reasonably close to the option's value.
        # A 10% relative tolerance is reasonable for a question with "~" values.
        if math.isclose(calculated_ratio, options[llm_answer_option], rel_tol=0.1):
            return "Correct"
        else:
            # This case is unlikely but possible if options are very far apart.
            return f"Incorrect. The LLM chose the closest option ({llm_answer_option}), but the calculated value ({calculated_ratio:.4f}) is not numerically close to the option's value ({options[llm_answer_option]})."
    else:
        return f"Incorrect. The calculated mass ratio is {calculated_ratio:.4f}. This is closest to option {closest_option} ({options[closest_option]}), not the provided answer of option {llm_answer_option} ({options[llm_answer_option]})."

# Run the check and print the result
print(check_answer())