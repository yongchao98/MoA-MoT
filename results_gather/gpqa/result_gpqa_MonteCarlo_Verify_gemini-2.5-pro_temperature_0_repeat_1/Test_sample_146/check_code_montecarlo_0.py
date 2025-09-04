import numpy as np

def check_annihilation_velocity():
    """
    Checks the correctness of the answer for the particle annihilation problem.
    The problem is: p + p_bar -> 2A+ + 2A-
    Given: m_A*c^2 = 300 MeV, antiproton is slowly moving.
    The LLM's answer is A) 0.77c.
    """
    # --- Define constants and the given answer ---
    # Rest mass energy of a proton in MeV (CODATA 2018)
    m_p_c2 = 938.272
    # Rest mass energy of particle A in MeV
    m_A_c2 = 300.0
    # The velocity from the selected answer (A)
    beta_answer = 0.77

    # --- Physics Calculation ---
    # From conservation of energy: 2 * m_p * c^2 = 4 * gamma * m_A * c^2
    # We can solve for the Lorentz factor, gamma.
    # gamma = (2 * m_p * c^2) / (4 * m_A * c^2)
    try:
        gamma_calculated = (2 * m_p_c2) / (4 * m_A_c2)
    except ZeroDivisionError:
        return "Incorrect. The mass of particle A (m_A) cannot be zero."

    # The Lorentz factor must be >= 1 for a real velocity.
    if gamma_calculated < 1:
        return (f"Incorrect. The calculated Lorentz factor is {gamma_calculated:.4f}, which is less than 1. "
                "This implies an imaginary velocity, meaning the process is not physically possible with the given masses.")

    # Now, solve for beta (v/c) using the definition of gamma:
    # gamma = 1 / sqrt(1 - beta^2)  =>  beta = sqrt(1 - 1/gamma^2)
    beta_squared = 1 - (1 / gamma_calculated**2)
    beta_calculated = np.sqrt(beta_squared)

    # --- Verification ---
    # Check if the calculated velocity is close to the answer's velocity.
    # A tolerance is used because the options in the question are rounded.
    # A tolerance of 0.01 is appropriate for options given to two decimal places.
    if np.isclose(beta_calculated, beta_answer, atol=0.01):
        return "Correct"
    else:
        # Provide a detailed reason for the failure.
        # Find which option is actually the closest.
        options = {'A': 0.77, 'B': 0.96, 'C': 0.91, 'D': 0.86}
        errors = {label: abs(beta_calculated - val) for label, val in options.items()}
        closest_option = min(errors, key=errors.get)
        
        return (f"Incorrect. The provided answer corresponds to a velocity of {beta_answer}c. "
                f"However, the velocity calculated from the conservation of energy is {beta_calculated:.4f}c. "
                f"This calculated value is not within a reasonable tolerance of the answer. "
                f"The correct answer should be option {closest_option} ({options[closest_option]}c).")

# Execute the check and print the result.
result = check_annihilation_velocity()
print(result)