import math

def check_correctness():
    """
    Checks the correctness of the provided answer by recalculating the velocity of particle A.
    """
    # --- Given Constants ---
    # Rest energy of a proton (m_p * c^2) in MeV. Using the CODATA 2018 value for precision.
    m_p_c2 = 938.27208816
    # Rest energy of particle A (m_A * c^2) in MeV, from the question.
    m_A_c2 = 300.0

    # --- Physics Calculation ---
    # 1. Calculate the Lorentz factor (gamma) from the conservation of energy.
    # 2 * m_p_c2 = 4 * gamma * m_A_c2
    # gamma = (2 * m_p_c2) / (4 * m_A_c2)
    try:
        gamma = (2 * m_p_c2) / (4 * m_A_c2)
    except ZeroDivisionError:
        return "Incorrect. The mass of particle A cannot be zero."

    # 2. Check if the process is physically possible (gamma must be >= 1).
    if gamma < 1:
        return f"Incorrect. The calculated Lorentz factor (gamma) is {gamma:.4f}, which is less than 1. This implies the initial energy is not sufficient to create the four A particles at rest, making the process physically impossible."

    # 3. Calculate the velocity (v) as a fraction of the speed of light (c), also known as beta.
    # beta = v/c = sqrt(1 - 1/gamma^2)
    beta = math.sqrt(1 - (1 / gamma**2))

    # --- Verification ---
    # The provided answer is <<<D>>>, which corresponds to 0.77c from the question's options.
    # Options: A) 0.91c, B) 0.96c, C) 0.86c, D) 0.77c
    expected_beta = 0.77

    # Check if the calculated beta is close to the expected beta.
    # A tolerance is used to account for rounding in the option value.
    tolerance = 0.01
    if abs(beta - expected_beta) < tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is D, which corresponds to a velocity of {expected_beta}c. "
                f"However, the calculation based on conservation of energy yields a velocity of {beta:.4f}c. "
                f"The calculation steps are: Initial Energy = 2 * {m_p_c2:.2f} MeV = {2*m_p_c2:.2f} MeV. "
                f"This requires a Lorentz factor (gamma) of {gamma:.4f}, which corresponds to a velocity of {beta:.4f}c.")

# Run the check
result = check_correctness()
print(result)