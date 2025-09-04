import numpy as np

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the quantum mechanics problem.
    It calculates the uncertainty Delta_Pz from first principles and compares it to the value
    of the option selected by the LLM.
    """
    # Let hbar = 1 for simplicity in calculation. The final answer will be in units of hbar.
    hbar = 1.0

    # Define the state vector |psi> as a column vector
    # psi = [-1/2, 1/sqrt(2), -1/2]^T
    psi = np.array([[-1/2], [1/np.sqrt(2)], [-1/2]])

    # Define the operator Pz as a matrix
    # Pz = hbar * [[1, 0, 0], [0, 0, 0], [0, 0, -1]]
    Pz = hbar * np.array([[1, 0, 0],
                          [0, 0, 0],
                          [0, 0, -1]])

    # --- Step 1: Verify the state is normalized ---
    # The bra vector <psi| is the conjugate transpose of |psi>
    psi_bra = psi.conj().T
    norm_squared = (psi_bra @ psi)[0, 0]
    if not np.isclose(norm_squared, 1.0):
        return f"Incorrect. The provided state vector is not normalized. <psi|psi> = {norm_squared:.4f}, but it should be 1."

    # --- Step 2: Calculate the expectation value <Pz> ---
    # <Pz> = <psi|Pz|psi>
    exp_Pz = (psi_bra @ Pz @ psi)[0, 0]

    # --- Step 3: Calculate the expectation value <Pz^2> ---
    # First, find the matrix for Pz^2
    Pz_squared = Pz @ Pz
    # Then, calculate <Pz^2> = <psi|Pz^2|psi>
    exp_Pz_squared = (psi_bra @ Pz_squared @ psi)[0, 0]

    # --- Step 4: Calculate the uncertainty Delta_Pz ---
    # (Delta_Pz)^2 = <Pz^2> - <Pz>^2
    variance = exp_Pz_squared - exp_Pz**2
    # The uncertainty is the square root of the variance
    uncertainty = np.sqrt(variance)

    # --- Step 5: Compare the calculated result with the LLM's answer ---
    # The options as listed in the final LLM's analysis
    options = {
        'A': np.sqrt(2) * hbar,
        'B': hbar / 2,
        'C': hbar / np.sqrt(2),
        'D': hbar
    }

    # The final answer provided by the LLM
    llm_answer_choice = 'C'
    
    # Check if the calculated uncertainty matches the value of the chosen option
    if llm_answer_choice not in options:
        return f"Incorrect. The provided answer '{llm_answer_choice}' is not a valid option."

    llm_answer_value = options[llm_answer_choice]

    if np.isclose(uncertainty, llm_answer_value):
        return "Correct"
    else:
        # If the answer is wrong, find the correct option
        correct_choice = None
        for choice, value in options.items():
            if np.isclose(uncertainty, value):
                correct_choice = choice
                break
        
        reason = f"The calculation from first principles gives:\n"
        reason += f"<Pz> = {exp_Pz:.4f} hbar\n"
        reason += f"<Pz^2> = {exp_Pz_squared:.4f} hbar^2\n"
        reason += f"Uncertainty Delta_Pz = sqrt(<Pz^2> - <Pz>^2) = sqrt({exp_Pz_squared:.4f} - {exp_Pz**2:.4f}) = {uncertainty:.4f} hbar.\n"
        reason += f"This calculated value ({uncertainty:.4f} hbar) corresponds to hbar/sqrt(2).\n"
        
        if correct_choice:
            reason += f"The correct option is {correct_choice}, which has a value of {options[correct_choice]:.4f} hbar.\n"
            reason += f"The provided answer was {llm_answer_choice}, which corresponds to a value of {llm_answer_value:.4f} hbar."
        else:
            reason += "The calculated value does not match any of the given options."
            
        return f"Incorrect. {reason}"

# Run the check and print the result
print(check_correctness())