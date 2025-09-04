import numpy as np

def check_answer():
    """
    Checks the correctness of the calculated uncertainty for the Pz operator.
    """
    # Set hbar = 1 for numerical calculations. It can be added back at the end.
    hbar = 1.0

    # Define the operators and the state vector from the problem description.
    # Px operator matrix
    Px = (hbar / np.sqrt(2)) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    # Pz operator matrix
    Pz = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])

    # The state of the system |psi> as a column vector
    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ])

    # The conjugate transpose of psi, <psi|
    psi_dagger = psi.conj().T

    # --- Constraint Checks ---

    # 1. Check if the state vector is normalized: <psi|psi> = 1
    norm = (psi_dagger @ psi).item()
    if not np.isclose(norm, 1.0):
        return f"The state vector is not normalized. <psi|psi> = {norm}, but it should be 1."

    # 2. Check if the state is an eigenstate of Px with eigenvalue -hbar
    # Calculate Px|psi>
    px_psi = Px @ psi
    # The expected result is -hbar * |psi>
    expected_eigenvalue_psi = -hbar * psi
    if not np.allclose(px_psi, expected_eigenvalue_psi):
        return f"The given state is not an eigenstate of Px with eigenvalue -hbar as claimed. Px|psi> = {px_psi.flatten()}, but expected {-hbar*psi.flatten()}."

    # --- Calculation of Uncertainty ---

    # 3. Calculate the expectation value of Pz: <Pz> = <psi|Pz|psi>
    exp_Pz = (psi_dagger @ Pz @ psi).item()

    # 4. Calculate the Pz^2 operator
    Pz_sq = Pz @ Pz

    # 5. Calculate the expectation value of Pz^2: <Pz^2> = <psi|Pz^2|psi>
    exp_Pz_sq = (psi_dagger @ Pz_sq @ psi).item()

    # 6. Calculate the variance: Var(Pz) = <Pz^2> - <Pz>^2
    # Check for negative variance which can happen due to floating point errors
    variance_Pz = exp_Pz_sq - exp_Pz**2
    if variance_Pz < 0:
        variance_Pz = 0

    # 7. The uncertainty is the square root of the variance
    uncertainty_Pz = np.sqrt(variance_Pz)

    # --- Final Answer Check ---

    # 8. The answer from option A is hbar/sqrt(2)
    expected_answer = hbar / np.sqrt(2)

    if np.isclose(uncertainty_Pz, expected_answer):
        return "Correct"
    else:
        return f"The calculated uncertainty is {uncertainty_Pz} * hbar, but the answer A suggests {expected_answer} * hbar. The calculation is as follows: <Pz> = {exp_Pz}, <Pz^2> = {exp_Pz_sq}, Variance = {variance_Pz}."

# Run the check
result = check_answer()
print(result)