import numpy as np

def check_correctness_of_answer():
    """
    This function checks the correctness of the answer to the quantum mechanics problem.
    It calculates the uncertainty Delta P_z and compares it to the expected answer.
    The expected correct answer is D) hbar/sqrt(2).
    """
    # Let hbar = 1 for simplicity in calculation. The final answer will be in units of hbar.
    hbar = 1.0

    # Define the operator P_x as a numpy array
    P_x = (hbar / np.sqrt(2)) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    # Define the operator P_z as a numpy array
    P_z = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])

    # Define the state of the system as a column vector
    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ])

    # --- Step 1: Verify the problem's constraints ---

    # Check if the state is normalized: <psi|psi> = 1
    # psi.conj().T is the Hermitian conjugate (bra vector)
    norm_squared = (psi.conj().T @ psi)[0, 0]
    if not np.isclose(norm_squared, 1.0):
        return f"Constraint failed: The given state vector is not normalized. <psi|psi> = {norm_squared:.4f}"

    # Check if the state is an eigenstate of P_x with eigenvalue -hbar
    # Calculate P_x |psi>
    Px_psi = P_x @ psi
    # The expected result is -hbar * |psi>
    expected_eigenvector_eq = -hbar * psi
    if not np.allclose(Px_psi, expected_eigenvector_eq):
        return (f"Constraint failed: The given state is not an eigenstate of P_x with eigenvalue -hbar.\n"
                f"P_x|psi> = {Px_psi.flatten()}\n"
                f"but -hbar|psi> = {expected_eigenvector_eq.flatten()}")

    # --- Step 2: Calculate the uncertainty Delta P_z ---

    # Calculate the expectation value of P_z: <P_z> = <psi|P_z|psi>
    exp_Pz = (psi.conj().T @ P_z @ psi)[0, 0]

    # Calculate the matrix for P_z^2
    Pz_squared = P_z @ P_z

    # Calculate the expectation value of P_z^2: <P_z^2> = <psi|P_z^2|psi>
    exp_Pz_squared = (psi.conj().T @ Pz_squared @ psi)[0, 0]

    # Calculate the variance: (Delta P_z)^2 = <P_z^2> - <P_z>^2
    variance_Pz = exp_Pz_squared - exp_Pz**2
    
    # The uncertainty is the square root of the variance
    delta_Pz = np.sqrt(variance_Pz)

    # --- Step 3: Compare with the expected answer ---
    # The correct answer is D) hbar/sqrt(2)
    expected_answer_value = hbar / np.sqrt(2)

    if np.isclose(delta_Pz, expected_answer_value):
        return "Correct"
    else:
        return (f"Incorrect: The calculated uncertainty Delta P_z is {delta_Pz:.4f} hbar, "
                f"which does not match the expected answer D) hbar/sqrt(2) ({expected_answer_value:.4f} hbar).\n"
                f"Details: <Pz> = {exp_Pz:.4f} hbar, <Pz^2> = {exp_Pz_squared:.4f} hbar^2.")

# Execute the checker function and print the result
result = check_correctness_of_answer()
print(result)