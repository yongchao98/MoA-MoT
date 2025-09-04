import numpy as np

def check_answer():
    """
    Checks the correctness of the given answer for the quantum mechanics problem.
    """
    # For numerical calculations, we can set hbar = 1.0 and add it back symbolically at the end.
    hbar = 1.0

    # Define the operators P_x and P_z as numpy arrays
    P_x = (hbar / np.sqrt(2)) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    P_z = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])

    # Define the state vector |psi> as a column vector
    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ])

    # --- Constraint Check 1: Verify the state is normalized ---
    # The norm squared should be <psi|psi> = 1
    norm_sq = np.vdot(psi, psi)
    if not np.isclose(norm_sq, 1.0):
        return f"Incorrect. The given state vector psi is not normalized. Its norm squared is {norm_sq:.4f}, but it should be 1."

    # --- Constraint Check 2: Verify the state is an eigenstate of P_x ---
    # The problem states psi is an eigenstate of P_x with eigenvalue -hbar.
    # Let's check if P_x|psi> = -hbar|psi>
    eigenvalue = -hbar
    Px_psi = P_x @ psi
    expected_result = eigenvalue * psi
    if not np.allclose(Px_psi, expected_result):
        return f"Incorrect. The problem statement is inconsistent. The given state is not an eigenstate of P_x with eigenvalue -hbar."

    # --- Calculation of Uncertainty ---

    # 1. Calculate the expectation value of P_z: <P_z> = <psi|P_z|psi>
    # <psi| is the conjugate transpose of |psi>
    psi_dagger = psi.T.conj()
    exp_Pz = (psi_dagger @ P_z @ psi)[0, 0]

    # 2. Calculate the expectation value of P_z^2: <P_z^2> = <psi|P_z^2|psi>
    Pz_squared = P_z @ P_z
    exp_Pz_squared = (psi_dagger @ Pz_squared @ psi)[0, 0]

    # 3. Calculate the variance and uncertainty (standard deviation)
    # Variance = <P_z^2> - <P_z>^2
    variance_Pz = exp_Pz_squared - (exp_Pz**2)
    
    # Uncertainty = sqrt(Variance)
    uncertainty_Pz = np.sqrt(variance_Pz)

    # --- Final Check ---
    # The given answer is B, which corresponds to hbar/sqrt(2)
    expected_value = hbar / np.sqrt(2)

    if np.isclose(uncertainty_Pz, expected_value):
        return "Correct"
    else:
        return (f"Incorrect. The calculated uncertainty is {uncertainty_Pz:.4f} * hbar. "
                f"The value from option B is hbar/sqrt(2) which is approximately {expected_value:.4f} * hbar. "
                f"The calculated expectation value <P_z> was {exp_Pz:.4f} * hbar. "
                f"The calculated expectation value <P_z^2> was {exp_Pz_squared:.4f} * hbar^2.")

# Run the check
result = check_answer()
print(result)