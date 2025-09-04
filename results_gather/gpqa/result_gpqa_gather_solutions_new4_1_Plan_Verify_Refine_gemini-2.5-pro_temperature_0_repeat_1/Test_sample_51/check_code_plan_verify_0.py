import numpy as np
import re

def check_correctness():
    """
    This function checks the correctness of the final answer by recalculating the physical quantity from the problem statement.
    """
    # --- Problem Parameters ---
    # Temperatures
    T_nospots = 6000.0  # K
    T_spots = 5500.0    # K
    # Wavelength
    wavelength_A = 1448.0  # Angstroms
    wavelength_m = wavelength_A * 1e-10  # meters

    # --- Physical Constants (CODATA 2018) ---
    h = 6.62607015e-34   # Planck's constant in J·s
    c = 2.99792458e8     # Speed of light in m/s
    k = 1.380649e-23      # Boltzmann constant in J/K

    # --- Calculation ---
    # The factor is F = exp[ (ΔE / k) * (1/T_spots - 1/T_nospots) ]
    # where ΔE = hc/λ
    try:
        exponent = ((h * c) / (wavelength_m * k)) * (1/T_spots - 1/T_nospots)
        calculated_factor = np.exp(exponent)
    except Exception as e:
        return f"An error occurred during the calculation: {e}"

    # --- Answer Verification ---
    # The final answer provided by the LLM
    llm_answer_str = "<<<A>>>"
    
    # Options from the question
    options = {
        "A": 4.5,
        "B": 7.8,
        "C": 2.9,
        "D": 1.1
    }

    # Extract the chosen option letter
    match = re.search(r'<<<([A-D])>>>', llm_answer_str)
    if not match:
        return f"The answer format is incorrect. Expected '<<<X>>>' but got '{llm_answer_str}'."
    
    chosen_option = match.group(1)
    
    # Get the value corresponding to the chosen option
    chosen_value = options.get(chosen_option)
    if chosen_value is None:
        return f"The chosen option '{chosen_option}' is not one of the valid options (A, B, C, D)."

    # Check if the calculated factor matches the chosen option's value within a tolerance
    tolerance = 0.1
    if abs(calculated_factor - chosen_value) < tolerance:
        return "Correct"
    else:
        # Find which option the calculation actually matches
        correct_option = 'None'
        for opt, val in options.items():
            if abs(calculated_factor - val) < tolerance:
                correct_option = opt
                break
        
        return (f"The answer is incorrect. "
                f"The calculated factor is approximately {calculated_factor:.3f}. "
                f"This value corresponds to option {correct_option} (~{options.get(correct_option, 'N/A')}). "
                f"The provided answer was option {chosen_option} (~{chosen_value}).")

# Run the check
result = check_correctness()
print(result)