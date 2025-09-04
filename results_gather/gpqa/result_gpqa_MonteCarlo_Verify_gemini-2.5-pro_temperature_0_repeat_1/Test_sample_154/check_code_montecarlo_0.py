import numpy as np

def check_answer():
    """
    Checks the correctness of the given answer by performing the quantum mechanical calculations.
    """
    # Set hbar = 1 for simplicity. The final answer will be in units of hbar.
    hbar = 1.0

    # --- 1. Define operators and state vector from the problem ---
    # Operator P_x
    Px = (hbar / np.sqrt(2)) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    # Operator P_z
    Pz = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])

    # State vector |psi>
    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ])

    # --- 2. Verify problem constraints ---
    # a) Check if the state is normalized
    norm_squared = np.vdot(psi, psi)
    if not np.isclose(norm_squared, 1.0):
        return f"Incorrect: The provided state vector is not normalized. The squared norm is {norm_squared}, but it should be 1."

    # b) Check if the state is an eigenstate of Px with eigenvalue -hbar
    # Calculate Px|psi>
    Px_psi = Px @ psi
    # The expected result is -hbar * |psi>
    expected_eigenstate_result = -hbar * psi
    if not np.allclose(Px_psi, expected_eigenstate_result):
        return (f"Incorrect: The provided state is not an eigenstate of Px with eigenvalue -hbar as claimed.\n"
                f"Px|psi> = {Px_psi.flatten()}\n"
                f"but -hbar*|psi> = {expected_eigenstate_result.flatten()}")

    # --- 3. Calculate expectation values ---
    # a) Calculate <Pz> = <psi|Pz|psi>
    # psi.conj().T is the Hermitian conjugate (bra vector <psi|)
    exp_Pz = (psi.conj().T @ Pz @ psi)[0, 0]

    # b) Calculate <Pz^2> = <psi|Pz^2|psi>
    Pz2 = Pz @ Pz
    exp_Pz2 = (psi.conj().T @ Pz2 @ psi)[0, 0]

    # --- 4. Calculate the uncertainty ---
    # The variance is <Pz^2> - <Pz>^2
    variance_Pz = exp_Pz2 - exp_Pz**2
    # The uncertainty (standard deviation) is the square root of the variance
    uncertainty_Pz = np.sqrt(variance_Pz)

    # --- 5. Compare with the given answer ---
    # The answer 'D' corresponds to hbar/sqrt(2)
    answer_value = hbar / np.sqrt(2)

    if np.isclose(uncertainty_Pz, answer_value):
        return "Correct"
    else:
        return (f"Incorrect: The calculated uncertainty does not match the answer.\n"
                f"Calculated <Pz> = {exp_Pz}\n"
                f"Calculated <Pz^2> = {exp_Pz2}\n"
                f"Calculated uncertainty Delta Pz = sqrt({exp_Pz2} - {exp_Pz}^2) = {uncertainty_Pz}\n"
                f"The value from answer D is {answer_value}.\n"
                f"These values do not match.")

# Run the check
result = check_answer()
print(result)