import numpy as np

def check_quantum_uncertainty_answer():
    """
    This function verifies the correctness of the given answer by:
    1. Defining the quantum operators (Px, Pz) and the state vector (psi) from the problem statement.
    2. Verifying that the given state vector is normalized and is an eigenstate of Px with eigenvalue -hbar, as stated in the problem.
    3. Calculating the expectation values <Pz> and <Pz^2>.
    4. Calculating the uncertainty Delta Pz = sqrt(<Pz^2> - <Pz>^2).
    5. Comparing the calculated results with the values provided in the answer.
    """
    # Use hbar = 1 for simplicity in calculation, and add it back symbolically at the end.
    # All results will be in units of hbar.
    hbar = 1.0
    sqrt2 = np.sqrt(2)

    # Define the matrix form of the operators
    Px = (hbar / sqrt2) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    Pz = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])

    # Define the state of the system as a column vector
    psi = np.array([
        [-1/2],
        [1/sqrt2],
        [-1/2]
    ])

    # --- Step 1: Verify the problem's constraints ---

    # Check if the state vector is normalized (<psi|psi> = 1)
    # <psi| is the conjugate transpose of |psi>
    psi_bra = psi.T.conj()
    norm_squared = (psi_bra @ psi).item()
    if not np.isclose(norm_squared, 1.0):
        return f"Constraint check failed: The state vector is not normalized. <psi|psi> = {norm_squared:.4f}, but it should be 1."

    # Check if the state is an eigenstate of Px with eigenvalue -hbar
    # We need to check if Px|psi> = -hbar|psi>
    px_on_psi = Px @ psi
    expected_result = -hbar * psi
    if not np.allclose(px_on_psi, expected_result):
        return f"Constraint check failed: The given state is not an eigenstate of Px with eigenvalue -hbar as claimed."

    # --- Step 2: Verify the calculations in the answer ---

    # The answer provides the following intermediate and final results:
    # <Pz> = 0
    # <Pz^2> = hbar^2 / 2
    # Delta Pz = hbar / sqrt(2)

    # Calculate the expectation value of Pz: <Pz> = <psi|Pz|psi>
    exp_Pz = (psi_bra @ Pz @ psi).item()
    if not np.isclose(exp_Pz, 0.0):
        return f"Calculation check failed: The expectation value of Pz is incorrect. Calculated <Pz> = {exp_Pz:.4f}, but the answer states it is 0."

    # Calculate the expectation value of Pz^2: <Pz^2> = <psi|Pz^2|psi>
    Pz_squared = Pz @ Pz
    exp_Pz_squared = (psi_bra @ Pz_squared @ psi).item()
    expected_exp_Pz_squared = hbar**2 / 2.0
    if not np.isclose(exp_Pz_squared, expected_exp_Pz_squared):
        return f"Calculation check failed: The expectation value of Pz^2 is incorrect. Calculated <Pz^2> = {exp_Pz_squared:.4f}, but the answer states it is hbar^2/2 = {expected_exp_Pz_squared:.4f}."

    # Calculate the uncertainty Delta Pz
    variance_Pz = exp_Pz_squared - exp_Pz**2
    delta_Pz = np.sqrt(variance_Pz)
    expected_delta_Pz = hbar / sqrt2

    if not np.isclose(delta_Pz, expected_delta_Pz):
        return f"Final result check failed: The uncertainty Delta Pz is incorrect. Calculated Delta Pz = {delta_Pz:.4f}, but the answer states it is hbar/sqrt(2) ≈ {expected_delta_Pz:.4f}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
# This will return "Correct" if the LLM's answer is valid, or a reason why it is not.
print(check_quantum_uncertainty_answer())