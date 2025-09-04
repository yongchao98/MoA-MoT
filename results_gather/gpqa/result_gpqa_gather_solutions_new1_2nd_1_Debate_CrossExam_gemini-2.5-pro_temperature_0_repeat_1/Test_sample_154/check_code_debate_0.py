import numpy as np

def check_answer():
    """
    Checks the correctness of the final answer for the quantum mechanics problem.
    """
    # Set hbar to 1 for simplicity, as it's a common factor.
    # The final answer will be in units of hbar.
    hbar = 1.0

    # Define the operators and the state vector from the question
    Pz = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ], dtype=float)

    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ], dtype=float)

    # --- Constraint Checks ---

    # 1. Check if the state is normalized (|psi|^2 = 1)
    norm_sq = np.vdot(psi, psi)
    if not np.isclose(norm_sq, 1.0):
        return f"Incorrect: The state vector is not normalized. <psi|psi> = {norm_sq}, but it should be 1."

    # 2. Check the consistency of the problem statement (optional but good practice)
    # The problem states psi is an eigenstate of Px with eigenvalue -hbar.
    Px = (hbar / np.sqrt(2)) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ], dtype=float)
    
    Px_psi = Px @ psi
    expected_Px_psi = -hbar * psi
    if not np.allclose(Px_psi, expected_Px_psi):
        return f"Incorrect: The problem statement is inconsistent. The given state is not an eigenstate of Px with eigenvalue -hbar."

    # --- Calculation of Uncertainty ---

    # Calculate the expectation value of Pz, <Pz>
    exp_Pz = np.vdot(psi, Pz @ psi)

    # Calculate the matrix for Pz^2
    Pz2 = Pz @ Pz

    # Calculate the expectation value of Pz^2, <Pz^2>
    exp_Pz2 = np.vdot(psi, Pz2 @ psi)

    # Calculate the variance and uncertainty (Delta Pz)
    variance = exp_Pz2 - exp_Pz**2
    uncertainty = np.sqrt(variance)

    # --- Verification of the Final Answer ---

    # The calculated numerical result should be hbar / sqrt(2)
    expected_uncertainty = hbar / np.sqrt(2)
    if not np.isclose(uncertainty, expected_uncertainty):
        return f"Incorrect: The calculated uncertainty is {uncertainty:.4f} hbar, but it should be {expected_uncertainty:.4f} hbar (hbar/sqrt(2))."

    # The final answer provided is 'D'.
    # The options provided in the final answer's analysis are:
    # A) hbar, B) sqrt(2)*hbar, C) hbar/2, D) hbar/sqrt(2)
    final_answer_choice = 'D'
    options = {
        'A': hbar,
        'B': np.sqrt(2) * hbar,
        'C': hbar / 2,
        'D': hbar / np.sqrt(2)
    }

    # Check if the value corresponding to the chosen answer 'D' matches the calculation
    if not np.isclose(options[final_answer_choice], uncertainty):
        return f"Incorrect: The final answer is '{final_answer_choice}', which corresponds to a value of {options[final_answer_choice]:.4f} hbar. However, the correctly calculated uncertainty is {uncertainty:.4f} hbar."

    # Check if the final answer's reasoning is sound.
    # The final answer correctly calculates the uncertainty as hbar/sqrt(2) and correctly maps it to option D from its list.
    # Therefore, the reasoning within the final answer is correct.
    
    # However, the question asks to decide the final answer based on the options provided in the *question prompt*.
    # The prompt's options are: A) hbar, B) sqrt(2)*hbar, C) hbar/2, D) hbar/sqrt(2)
    # The final answer's analysis uses the same options. So the mapping is correct.
    
    return "Correct"

# Run the check
result = check_answer()
print(result)