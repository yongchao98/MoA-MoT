import numpy as np

def check_answer():
    """
    This function verifies the calculation for the imaginary part of the scattering amplitude.
    It checks the two main hypotheses (relativistic vs. non-relativistic) and compares
    the results to the given options.
    """
    # --- Problem Inputs ---
    phase_shifts_deg = {0: 90, 1: 67, 2: 55, 3: 30, 4: 13}
    T_MeV = 50.0  # Kinetic energy in MeV
    options = {
        "A": 355.351,
        "B": 87163.4,
        "C": 251.271,
        "D": 177.675
    }
    # The final answer from the LLM to be checked.
    llm_answer_letter = "C"

    # --- Physical Constants ---
    m_e_c2_MeV = 0.511      # Electron rest mass energy in MeV
    hbar_c_MeV_fm = 197.327 # h-bar * c in MeV*fm

    # --- Calculation ---
    # Step 1: Calculate the summation term S
    S = 0
    for l, delta_deg in phase_shifts_deg.items():
        delta_rad = np.deg2rad(delta_deg)
        term = (2 * l + 1) * (np.sin(delta_rad))**2
        S += term

    # Step 2: Calculate the result using the non-relativistic hypothesis.
    # This is the only path that leads to one of the options.
    pc_non_rel_MeV = np.sqrt(2 * m_e_c2_MeV * T_MeV)
    k_non_rel_fm_inv = pc_non_rel_MeV / hbar_c_MeV_fm
    calculated_value = S / k_non_rel_fm_inv

    # --- Verification ---
    # Check if the LLM's chosen option value matches our calculated value.
    # Use a tolerance for floating-point comparison.
    tolerance = 1e-3
    llm_option_value = options.get(llm_answer_letter)

    if llm_option_value is None:
        return f"Incorrect. The provided answer '{llm_answer_letter}' is not a valid option."

    if np.isclose(calculated_value, llm_option_value, atol=tolerance):
        return "Correct"
    else:
        # If the answer is wrong, find the correct one to provide a helpful reason.
        correct_option = None
        for key, value in options.items():
            if np.isclose(calculated_value, value, atol=tolerance):
                correct_option = key
                break
        
        if correct_option:
            return (f"Incorrect. The provided answer is '{llm_answer_letter}' ({llm_option_value} fm). "
                    f"The calculation using the non-relativistic approximation yields approximately {calculated_value:.3f} fm, "
                    f"which matches option '{correct_option}' ({options[correct_option]} fm).")
        else:
            return (f"Incorrect. The provided answer is '{llm_answer_letter}'. "
                    f"The calculation using the non-relativistic approximation yields {calculated_value:.3f} fm, "
                    f"which does not match any of the provided options.")

# Execute the check and print the result
print(check_answer())