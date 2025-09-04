import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the gamma-ray annihilation problem.
    
    The problem is to find the threshold energy of a high-energy gamma-ray (E_gamma)
    colliding with a CMB photon (E_CMB) to produce an electron-positron pair.
    The reaction is: gamma + gamma_CMB -> e+ + e-

    The threshold condition for a head-on collision (which gives the minimum E_gamma) is:
    E_gamma = (m_e * c^2)^2 / E_CMB
    """
    
    # --- Step 1: Define constants and the given values ---
    
    # Electron rest mass energy (m_e * c^2) in eV.
    # The LLM used 0.511 MeV, which is a standard rounded value. We will use a more precise one for verification.
    # PDG 2022 value: 0.51099895000 MeV
    m_e_c2_eV = 0.51099895 * 1e6
    
    # Average CMB photon energy in eV, as given in the question.
    E_CMB_eV = 1e-3
    
    # The LLM's chosen answer is A, which corresponds to 2.6 * 1e5 GeV.
    llm_answer_value_GeV = 2.6 * 1e5
    
    # --- Step 2: Perform the theoretical calculation ---
    
    # Calculate the threshold energy in eV using the derived formula.
    try:
        calculated_E_gamma_eV = (m_e_c2_eV**2) / E_CMB_eV
    except ZeroDivisionError:
        return "Error: Division by zero. The CMB photon energy cannot be zero."
    except Exception as e:
        return f"An error occurred during calculation: {e}"
        
    # Convert the calculated energy from eV to GeV for comparison.
    # 1 GeV = 1e9 eV
    calculated_E_gamma_GeV = calculated_E_gamma_eV / 1e9
    
    # --- Step 3: Compare the calculated result with the LLM's answer ---
    
    # The LLM's reasoning and calculation are based on the same formula and similar constants.
    # The LLM calculated E_gamma ≈ 2.61e5 GeV and correctly matched it to option A (2.6e5 GeV).
    # We will verify this by checking if our more precise calculation is close to the value of option A.
    # A relative tolerance of 2% is reasonable to account for rounding in the option's value and the constants used.
    
    is_correct = math.isclose(calculated_E_gamma_GeV, llm_answer_value_GeV, rel_tol=0.02)
    
    if is_correct:
        return "Correct"
    else:
        # If the answer is not correct, provide a detailed explanation.
        reason = (
            f"The answer is incorrect. The chosen option A does not sufficiently match the theoretical calculation.\n"
            f"Constraint: The threshold energy for pair production is given by E_gamma = (m_e*c^2)^2 / E_CMB.\n"
            f"Calculation:\n"
            f"  - Electron rest mass energy (m_e*c^2) = {m_e_c2_eV:.4e} eV\n"
            f"  - CMB photon energy (E_CMB) = {E_CMB_eV:.1e} eV\n"
            f"  - Calculated E_gamma = ({m_e_c2_eV:.4e})^2 / {E_CMB_eV:.1e} = {calculated_E_gamma_eV:.4e} eV\n"
            f"  - Calculated E_gamma in GeV = {calculated_E_gamma_GeV:.4e} GeV\n"
            f"Comparison:\n"
            f"  - The calculated value is {calculated_E_gamma_GeV:.4e} GeV.\n"
            f"  - The value from the chosen option A is {llm_answer_value_GeV:.4e} GeV.\n"
            f"  - These values are not close enough (relative difference is {abs(calculated_E_gamma_GeV - llm_answer_value_GeV)/llm_answer_value_GeV:.2%}, which is outside the tolerance)."
        )
        return reason
