import numpy as np

def check_correctness():
    """
    This function verifies the answer to the particle physics problem.
    It calculates the theoretical resolution required and checks if the provided
    answer satisfies the problem's constraints.
    """
    # --- Given constants and problem values ---
    E_GeV = 27.0  # Energy in GeV
    m_GeV = 3.41  # Mass in GeV/c^2
    tau_0 = 8e-16  # Proper lifetime in seconds
    c = 299792458.0  # Speed of light in m/s
    target_fraction = 0.30

    # --- Options provided in the question ---
    options = {
        "A": 2.08e-9,
        "B": 2.08e-3,
        "C": 2.08e-1,
        "D": 2.08e-6,
    }
    
    # The answer from the LLM to be checked
    llm_answer_key = "D"
    
    if llm_answer_key not in options:
        return f"Incorrect. The provided answer key '{llm_answer_key}' is not one of the options A, B, C, or D."
        
    llm_answer_value = options[llm_answer_key]

    # --- Physics Calculations ---
    # 1. Lorentz factor (gamma)
    gamma = E_GeV / m_GeV
    
    # 2. Mean decay length (L) in the lab frame.
    # L = v * tau_lab = (beta*c) * (gamma*tau_0)
    # Using the identity beta*gamma = sqrt(gamma^2 - 1) for better precision.
    mean_decay_length = c * tau_0 * np.sqrt(gamma**2 - 1)

    # --- Verification ---
    # The question asks for the minimum resolution 'R' to observe AT LEAST 30% of decays.
    # This means P(decay length > R) >= 0.30, where P(l > R) = exp(-R / L).
    
    # First, check if the LLM's chosen answer satisfies the core condition.
    fraction_observed_for_llm_answer = np.exp(-llm_answer_value / mean_decay_length)

    if fraction_observed_for_llm_answer < target_fraction:
        return (f"Incorrect. The resolution from answer {llm_answer_key} ({llm_answer_value:.2e} m) "
                f"would observe only {fraction_observed_for_llm_answer:.2%} of decays, "
                f"which is less than the required {target_fraction:.0%}. "
                f"The calculated mean decay length is {mean_decay_length:.3e} m.")

    # The answer satisfies the "at least 30%" condition.
    # Now, we must check if it's the best choice among the valid options.
    # The "minimum resolution required" corresponds to the largest possible resolution value R
    # from the options that still satisfies the condition.
    
    valid_options = []
    for key, r_val in options.items():
        frac = np.exp(-r_val / mean_decay_length)
        if frac >= target_fraction:
            valid_options.append((key, r_val))
            
    if not valid_options:
        # This case should not be reached if the LLM answer was already validated above.
        return "Error: No options satisfy the condition."

    # Find the valid option with the largest resolution value (i.e., the tightest constraint).
    best_option = max(valid_options, key=lambda item: item[1])
    
    if llm_answer_key == best_option[0]:
        return "Correct"
    else:
        return (f"Incorrect. While answer {llm_answer_key} observes more than 30% of decays, "
                f"it is not the best answer. The question asks for the minimum resolution required, "
                f"which corresponds to the largest possible resolution value from the options that works. "
                f"Option {best_option[0]} ({best_option[1]:.2e} m) also satisfies the condition and is a larger value, "
                f"making it a better fit for the 'minimum required' threshold.")

# print(check_correctness())