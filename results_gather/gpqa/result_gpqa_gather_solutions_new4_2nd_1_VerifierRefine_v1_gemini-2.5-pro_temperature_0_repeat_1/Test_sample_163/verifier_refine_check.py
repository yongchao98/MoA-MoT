import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the binary star mass ratio problem.
    """
    # --- Data from the question ---
    # System 1
    P1 = 2.0  # Period in years
    K1a = 10.0 # RV amplitude in km/s
    K1b = 5.0  # RV amplitude in km/s

    # System 2
    P2 = 1.0  # Period in years
    K2a = 15.0 # RV amplitude in km/s
    K2b = 10.0 # RV amplitude in km/s

    # --- The final answer provided by the LLM ---
    # The LLM's final conclusion is 'C' which corresponds to ~0.4
    llm_answer_letter = 'C'
    llm_reasoning_value = 0.432
    
    # The options as presented in the final LLM response being checked
    options = {
        'A': 1.2,
        'B': 0.6,
        'C': 0.4,
        'D': 0.7
    }

    # --- Verification ---
    # The total mass of an eclipsing binary system is proportional to P * (K_total)^3.
    # The ratio of masses is therefore (P1 * (K1_total)^3) / (P2 * (K2_total)^3).

    # 1. Calculate the sum of radial velocity amplitudes for each system.
    K_total1 = K1a + K1b
    K_total2 = K2a + K2b

    if K_total1 != 15.0 or K_total2 != 25.0:
        return f"Incorrect. The sum of radial velocities is miscalculated. Expected K_total1=15 and K_total2=25, but got {K_total1} and {K_total2}."

    # 2. Calculate the mass ratio.
    try:
        calculated_ratio = (P1 * (K_total1**3)) / (P2 * (K_total2**3))
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # 3. Check if the LLM's calculated value in its reasoning is correct.
    if not math.isclose(calculated_ratio, llm_reasoning_value, rel_tol=1e-3):
        return (f"Incorrect. The numerical calculation in the provided answer's reasoning is wrong. "
                f"The correct ratio is {calculated_ratio:.3f}, but the answer's reasoning states it is {llm_reasoning_value}.")

    # 4. Find the option that is numerically closest to the calculated ratio.
    closest_option_letter = None
    min_difference = float('inf')

    for letter, value in options.items():
        difference = abs(calculated_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option_letter = letter
            
    # 5. Check if the LLM's final answer letter matches the closest option.
    if llm_answer_letter == closest_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The final answer choice is wrong. "
                f"The calculated ratio is {calculated_ratio:.3f}, which is closest to the value of option '{closest_option_letter}' ({options[closest_option_letter]}). "
                f"The provided answer selected option '{llm_answer_letter}' ({options[llm_answer_letter]}).")

# The final output of the code will be the return value of the check_answer() function.
# To display the result, you would typically call: print(check_answer())
result = check_answer()
print(result)