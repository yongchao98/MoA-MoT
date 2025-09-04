import math

def check_answer():
    """
    Checks the correctness of the answer for the exoplanet temperature ratio problem.
    """
    # --- Given values from the problem ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- Derivation and Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 is derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)
    # where K is the radial velocity semi-amplitude.
    # Since K is proportional to the Doppler shift Δλ, we have:
    # K1 / K2 = Δλ₁ / Δλ₂
    # So, the final formula is:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ₁ / Δλ₂)

    mass_ratio = m_p2 / m_p1
    shift_ratio = delta_lambda1 / delta_lambda2
    
    calculated_ratio = mass_ratio * shift_ratio
    
    # --- Options from the question ---
    options = {
        "A": 0.53,
        "B": 0.98,
        "C": 1.30,
        "D": 1.05
    }
    
    # The provided answer is 'A'
    provided_answer_label = 'A'
    provided_answer_value = options[provided_answer_label]

    # --- Verification ---
    # Check if the calculated ratio is close to the value of the chosen option 'A'
    if math.isclose(calculated_ratio, provided_answer_value, rel_tol=1e-2):
        return "Correct"
    else:
        # Find the correct option based on the calculation
        correct_option = ''
        min_diff = float('inf')
        for label, value in options.items():
            diff = abs(calculated_ratio - value)
            if diff < min_diff:
                min_diff = diff
                correct_option = label
        
        return (f"Incorrect. The calculated ratio is {calculated_ratio:.4f}. "
                f"This value is approximately {options[correct_option]:.2f}, which corresponds to option {correct_option}. "
                f"The provided answer was {provided_answer_label} ({provided_answer_value:.2f}).")

# Run the check
result = check_answer()
print(result)