import numpy as np
import re

def check_correctness():
    """
    This function checks the correctness of the final answer to the physics problem.
    It recalculates the angle based on the principles of special relativity and compares
    it to the provided answer choice.
    """
    
    # Given parameters from the question
    gamma_e = 4  # Lorentz factor for the electron
    gamma_p = 2  # Lorentz factor for the positron
    
    # The options provided in the question
    options = {
        'A': 172,
        'B': 96,
        'C': 138,
        'D': 74
    }
    
    # The final answer provided by the LLM
    llm_answer_string = "<<<C>>>"

    # --- Step 1: Calculate the correct angle ---
    # The formula is derived from the conservation of energy and momentum.
    # cos(theta) = (sqrt(gamma_e^2 - 1) - sqrt(gamma_p^2 - 1)) / (gamma_e + gamma_p)
    # The angle between the photons is 2 * theta.
    try:
        numerator = np.sqrt(gamma_e**2 - 1) - np.sqrt(gamma_p**2 - 1)
        denominator = gamma_e + gamma_p
        cos_theta = numerator / denominator
        
        # Calculate theta in radians, then convert to degrees
        theta_rad = np.arccos(cos_theta)
        theta_deg = np.degrees(theta_rad)
        
        # The total angle between the two photons is 2 * theta
        calculated_angle = 2 * theta_deg
    except Exception as e:
        return f"An error occurred during the calculation: {e}"

    # --- Step 2: Parse and evaluate the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_string)
    if not match:
        return f"The answer format is invalid: {llm_answer_string}"
    
    chosen_option_letter = match.group(1)
    chosen_option_value = options[chosen_option_letter]

    # --- Step 3: Compare the calculated result with the chosen option ---
    # We check if the chosen option is the closest integer value to our precise calculation.
    
    # Find which option is numerically closest to our calculated angle
    closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_angle))

    if chosen_option_letter == closest_option_letter:
        return "Correct"
    else:
        closest_option_value = options[closest_option_letter]
        return (f"Incorrect. The calculation shows the angle between the photons is approximately {calculated_angle:.2f} degrees. "
                f"This value is closest to option {closest_option_letter}) {closest_option_value}. "
                f"The provided answer was {chosen_option_letter}) {chosen_option_value}.")

# Run the check and print the result
result = check_correctness()
print(result)