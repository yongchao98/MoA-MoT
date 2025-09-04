import math

def check_binary_star_mass_ratio():
    """
    Checks the correctness of the calculated mass ratio for the two binary star systems.
    """
    # --- Define problem parameters ---
    # System 1
    P1 = 2.0  # years
    K1_1 = 10.0 # km/s
    K1_2 = 5.0  # km/s

    # System 2
    P2 = 1.0  # year
    K2_1 = 15.0 # km/s
    K2_2 = 10.0 # km/s

    # --- Define the LLM's answer and the available options ---
    llm_answer_option = 'B'
    options = {'A': 0.7, 'B': 0.4, 'C': 0.6, 'D': 1.2}
    
    # --- Perform the calculation based on established physics ---
    # For an eclipsing binary system (inclination i ≈ 90°, so sin(i) ≈ 1),
    # the total mass M is proportional to P * (K1 + K2)³.
    # M_total ∝ P * (K_sum)³
    # Therefore, the ratio of the masses is:
    # M1 / M2 = (P1 / P2) * ((K1_1 + K1_2) / (K2_1 + K2_2))³

    # Calculate the sum of radial velocity amplitudes for each system
    K_sum_1 = K1_1 + K1_2  # 10 + 5 = 15 km/s
    K_sum_2 = K2_1 + K2_2  # 15 + 10 = 25 km/s

    # Calculate the exact mass ratio
    try:
        calculated_ratio = (P1 / P2) * (K_sum_1 / K_sum_2)**3
    except ZeroDivisionError:
        return "Error: Division by zero in the calculation. P2 or K_sum_2 cannot be zero."
    except Exception as e:
        return f"An unexpected error occurred during calculation: {e}"

    # --- Verify the LLM's answer ---
    # Find the option that is numerically closest to our calculated result.
    # This is necessary because the options are approximate (indicated by ~).
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # Check if the LLM's chosen option is indeed the closest one.
    if llm_answer_option == closest_option_key:
        return "Correct"
    else:
        return (f"Incorrect. The calculated mass ratio is {calculated_ratio:.4f}. "
                f"The closest option to this value is '{closest_option_key}' ({options[closest_option_key]}), "
                f"but the provided answer was '{llm_answer_option}' ({options[llm_answer_option]}).")

# Run the checker and print the result
result = check_binary_star_mass_ratio()
print(result)