import math

def check_annihilation_answer():
    """
    Checks the correctness of the given answer for the particle annihilation problem.
    """
    # --- Constants and Given Values ---
    # Rest energy of a proton in MeV. An antiproton has the same rest energy.
    # Using a standard value from physics data.
    m_p_c2 = 938.272  # MeV
    
    # Rest energy of particle A in MeV, as given in the question.
    m_A_c2 = 300.0    # MeV
    
    # The velocity from the proposed answer (C: 0.77c), as a fraction of c.
    beta_answer = 0.77
    
    # --- Step 1: Calculate the initial energy of the system ---
    # Assuming the proton and antiproton are effectively at rest.
    E_initial = 2 * m_p_c2
    
    # --- Step 2: Calculate the final energy based on the proposed answer ---
    # First, calculate the Lorentz factor (gamma) for the given velocity.
    try:
        gamma = 1 / math.sqrt(1 - beta_answer**2)
    except ValueError:
        return "The velocity from the answer is greater than or equal to the speed of light, which is physically impossible."

    # The final energy is the sum of the total energies of the four A particles.
    E_final_from_answer = 4 * gamma * m_A_c2
    
    # --- Step 3: Check for energy conservation ---
    # We use a relative tolerance to account for potential rounding in the problem's
    # given values or the multiple-choice options. A 1.5% tolerance is reasonable.
    if math.isclose(E_initial, E_final_from_answer, rel_tol=0.015):
        return "Correct"
    else:
        # If the answer is incorrect, provide the reason.
        # Let's calculate the correct velocity for comparison.
        # From 2*m_p*c^2 = 4*gamma*m_A*c^2, we get gamma = (2*m_p*c^2) / (4*m_A*c^2)
        gamma_correct = (2 * m_p_c2) / (4 * m_A_c2)
        # From gamma = 1/sqrt(1-beta^2), we get beta = sqrt(1 - 1/gamma^2)
        beta_correct = math.sqrt(1 - 1/gamma_correct**2)
        
        reason = (
            f"The answer is incorrect because it violates the principle of conservation of energy.\n"
            f"Initial Energy (from proton/antiproton rest mass): {E_initial:.2f} MeV.\n"
            f"Final Energy (calculated from the answer's velocity v={beta_answer}c): {E_final_from_answer:.2f} MeV.\n"
            f"These two energy values are not equal within a reasonable tolerance.\n"
            f"The correct velocity should be approximately {beta_correct:.3f}c, which corresponds to the calculated final energy of {E_initial:.2f} MeV."
        )
        return reason

# Run the check
result = check_annihilation_answer()
print(result)