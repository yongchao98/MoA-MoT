import numpy as np

def check_planet_temperature_ratio():
    """
    This function checks the correctness of the LLM's answer by independently
    calculating the ratio of equilibrium temperatures based on the problem's parameters.
    """
    # --- Given parameters from the question ---
    # Planet 1 mass in Earth masses
    mp1_earth = 7.0
    # Planet 2 mass in Earth masses
    mp2_earth = 5.0
    # Doppler shift for Planet 1 in Angstroms
    d_lambda1 = 0.03
    # Doppler shift for Planet 2 in Angstroms
    d_lambda2 = 0.04

    # The multiple-choice options
    options = {'A': 0.98, 'B': 1.05, 'C': 0.53, 'D': 1.30}
    # The final answer provided by the LLM
    llm_answer_key = 'C'

    # --- Physics Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 can be derived as:
    # T_eq1 / T_eq2 = (Mp2 / Mp1) * (d_lambda1 / d_lambda2)
    # This derivation correctly ignores irrelevant information like star properties
    # and planet radii, as these factors cancel out in the ratio.

    # Calculate the ratio of planet masses (Mp2 / Mp1)
    mass_ratio = mp2_earth / mp1_earth

    # Calculate the ratio of Doppler shifts (d_lambda1 / d_lambda2)
    shift_ratio = d_lambda1 / d_lambda2

    # Calculate the final theoretical ratio of equilibrium temperatures
    calculated_ratio = mass_ratio * shift_ratio

    # --- Verification ---
    # Check if the provided answer key is valid
    if llm_answer_key not in options:
        return f"Incorrect. The provided answer key '{llm_answer_key}' is not one of the valid options {list(options.keys())}."

    # Find the option that is numerically closest to our calculated result.
    # This is the most robust way to select a multiple-choice answer from a calculation.
    differences = {key: abs(value - calculated_ratio) for key, value in options.items()}
    closest_option_key = min(differences, key=differences.get)

    # Compare the LLM's answer with the derived closest option
    if closest_option_key == llm_answer_key:
        # The LLM correctly identified the closest numerical option.
        # Let's check if the rounding is reasonable.
        error = differences[closest_option_key]
        # The relative error between the calculated value and the option value.
        relative_error = error / options[closest_option_key]
        
        # A small relative error (e.g., < 2%) indicates a good match due to rounding.
        # 15/28 is ~0.5357. |0.5357 - 0.53| / 0.53 is ~1.07%
        if relative_error < 0.02:
            return "Correct"
        else:
            # This case would catch if the "closest" option was still very far off.
            return (f"Correct. The LLM chose the mathematically closest option '{llm_answer_key}'. "
                    f"The calculated value is {calculated_ratio:.4f} and the option value is {options[llm_answer_key]}.")
    else:
        # The LLM chose the wrong option.
        return (f"Incorrect. The calculated ratio is {calculated_ratio:.4f}. "
                f"This value is closest to option {closest_option_key} ({options[closest_option_key]}), "
                f"but the provided answer was option {llm_answer_key} ({options[llm_answer_key]}).")

# Execute the check
result = check_planet_temperature_ratio()
print(result)