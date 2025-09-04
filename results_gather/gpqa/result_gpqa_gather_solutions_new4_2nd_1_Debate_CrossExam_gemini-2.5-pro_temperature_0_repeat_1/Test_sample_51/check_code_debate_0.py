import numpy as np
from scipy import constants

def check_correctness():
    """
    Checks the correctness of the final answer by recalculating the factor.
    """
    # --- Define physical constants from scipy for accuracy ---
    h = constants.h  # Planck's constant in J·s
    c = constants.c  # Speed of light in m/s
    k = constants.k  # Boltzmann constant in J/K

    # --- Define the given parameters from the problem ---
    T_nospots = 6000.0  # Temperature without spots (T1) in Kelvin
    T_spots = 5500.0    # Temperature with spots (T2) in Kelvin
    wavelength_A = 1448  # Wavelength in Angstroms
    wavelength_m = wavelength_A * 1e-10 # Convert wavelength to meters

    # --- Define the options as presented in the final answer ---
    # A) ~1.1, B) ~7.8, C) ~4.5, D) ~2.9
    options = {'A': 1.1, 'B': 7.8, 'C': 4.5, 'D': 2.9}
    chosen_answer_option = 'C'
    
    # --- Perform the calculation ---
    # The formula for the factor is:
    # Factor = exp[ (ΔE / k) * (1/T_spots - 1/T_nospots) ]
    # where ΔE = hc/λ
    try:
        # Calculate the energy term (ΔE/k)
        energy_term = (h * c) / (wavelength_m * k)

        # Calculate the temperature term
        temp_term = (1 / T_spots) - (1 / T_nospots)

        # Calculate the exponent
        exponent = energy_term * temp_term

        # Calculate the final factor
        calculated_factor = np.exp(exponent)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verify the result ---
    # 1. Check if the calculated factor matches the value of the chosen option.
    # We use np.isclose to account for minor rounding differences.
    # A relative tolerance of 5% is reasonable for "approximately equal to".
    expected_value = options[chosen_answer_option]
    
    if np.isclose(calculated_factor, expected_value, rtol=0.05):
        return "Correct"
    else:
        # If the calculation does not match the chosen option, find the correct option.
        correct_option = None
        for option, value in options.items():
            if np.isclose(calculated_factor, value, rtol=0.05):
                correct_option = option
                break
        
        reason = (f"Incorrect. The calculation yields a factor of approximately {calculated_factor:.2f}. "
                  f"This value corresponds to option {correct_option} (~{options.get(correct_option, 'N/A')}). "
                  f"The provided answer selected option {chosen_answer_option} (~{expected_value}).")
        return reason

# Run the check
result = check_correctness()
print(result)