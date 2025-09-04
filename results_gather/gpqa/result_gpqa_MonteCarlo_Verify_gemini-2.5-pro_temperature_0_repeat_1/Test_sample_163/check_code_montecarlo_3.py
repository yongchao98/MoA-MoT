import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the binary star mass ratio problem.
    It calculates the theoretical mass ratio based on the given data and compares it to the answer.
    """
    try:
        # --- Given Data from the Question ---
        # System 1
        P1 = 2.0  # Period in years
        K1_1 = 10.0  # Radial velocity amplitude of star 1 in km/s
        K1_2 = 5.0   # Radial velocity amplitude of star 2 in km/s

        # System 2
        P2 = 1.0  # Period in years
        K2_1 = 15.0  # Radial velocity amplitude of star 1 in km/s
        K2_2 = 10.0  # Radial velocity amplitude of star 2 in km/s

        # The provided answer is 'B', which corresponds to a value of approximately 0.4
        llm_answer_option = 'B'
        options = {'A': 0.7, 'B': 0.4, 'C': 0.6, 'D': 1.2}
        
        # --- Physics Principles and Calculation ---
        # For a double-lined spectroscopic binary, the total mass (M = m1 + m2) is related to the
        # orbital period (P) and the radial velocity semi-amplitudes (K1, K2) by the formula:
        # M = (P / (2 * pi * G)) * ((K1 + K2) / sin(i))^3
        # where G is the gravitational constant and 'i' is the orbital inclination.

        # A critical constraint is that both systems exhibit eclipses. For eclipses to occur,
        # the inclination angle 'i' must be close to 90 degrees.
        # Therefore, sin(i) is approximately 1, and sin^3(i) is also approximately 1.
        # This simplifies the mass relationship to: M is proportional to P * (K1 + K2)^3.

        # We need to find the ratio of the mass of system_1 (M1) to the mass of system_2 (M2).
        # The ratio is M1 / M2 = [P1 * (K1_1 + K1_2)^3] / [P2 * (K2_1 + K2_2)^3].
        # The constants (like 1/(2*pi*G)) and the units cancel out in the ratio.

        # Sum of velocity amplitudes for system 1
        K_sum_1 = K1_1 + K1_2  # 10 + 5 = 15 km/s

        # Sum of velocity amplitudes for system 2
        K_sum_2 = K2_1 + K2_2  # 15 + 10 = 25 km/s

        # Calculate the mass ratio
        calculated_ratio = (P1 * (K_sum_1**3)) / (P2 * (K_sum_2**3))
        # calculated_ratio = (2.0 * (15**3)) / (1.0 * (25**3))
        # calculated_ratio = (2 * 3375) / 15625
        # calculated_ratio = 6750 / 15625
        # calculated_ratio = 0.432

        # --- Verification ---
        # Find which of the given options is closest to our calculated value.
        closest_option = min(options.keys(), key=lambda option: abs(options[option] - calculated_ratio))

        if llm_answer_option == closest_option:
            # The LLM's answer corresponds to the option that is mathematically closest to the calculated result.
            return "Correct"
        else:
            return (f"Incorrect. The provided answer is {llm_answer_option}, but the calculation shows the correct option is {closest_option}.\n"
                    f"Reasoning:\n"
                    f"1. The total mass M of an eclipsing binary system is proportional to P * (K1 + K2)^3, where P is the period and K1, K2 are the radial velocity amplitudes.\n"
                    f"2. For System 1: P1 = 2 years, K1_sum = 10 + 5 = 15 km/s.\n"
                    f"3. For System 2: P2 = 1 year, K2_sum = 15 + 10 = 25 km/s.\n"
                    f"4. The mass ratio M1/M2 = [P1 * (K1_sum)^3] / [P2 * (K2_sum)^3].\n"
                    f"5. Plugging in the values: M1/M2 = [2 * (15)^3] / [1 * (25)^3] = [2 * 3375] / [15625] = 6750 / 15625 = 0.432.\n"
                    f"6. The calculated value 0.432 is closest to option {closest_option} (~{options[closest_option]}), not option {llm_answer_option} (~{options[llm_answer_option]}).")

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# The final output of the code block will be the return value of this function.
# print(check_correctness())