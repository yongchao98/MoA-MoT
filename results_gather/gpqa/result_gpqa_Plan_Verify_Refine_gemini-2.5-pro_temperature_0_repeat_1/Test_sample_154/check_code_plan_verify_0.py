import numpy as np

def check_correctness():
    """
    Checks the correctness of the LLM's answer by recalculating the uncertainty of Pz.
    """
    # Set hbar = 1 for numerical calculations. It will be reintroduced at the end.
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

    # The state vector of the system, given as a column vector.
    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ])

    # --- Constraint Check 1: Verify the state is an eigenstate of Px ---
    # The problem states the system is in an eigenstate of Px with eigenvalue -ħ.
    # Let's check if Px|ψ> = -ħ|ψ>.
    eigenvalue = -hbar
    result_vector = Px @ psi
    expected_vector = eigenvalue * psi
    if not np.allclose(result_vector, expected_vector):
        return (f"Incorrect: The provided state vector is not an eigenstate of Px with eigenvalue -ħ as claimed.\n"
                f"Px|ψ> = {result_vector.flatten()}\n"
                f"but -ħ|ψ> = {expected_vector.flatten()}")

    # --- Constraint Check 2: Verify the state is normalized ---
    # The norm of a state vector must be 1. <ψ|ψ> = 1.
    norm_squared = (psi.conj().T @ psi).item()
    if not np.isclose(norm_squared, 1.0):
        return f"Incorrect: The provided state vector is not normalized. <ψ|ψ> = {norm_squared:.4f}, but it should be 1."

    # --- Recalculate the uncertainty ΔPz ---
    # The formula is ΔPz = sqrt(<Pz^2> - <Pz>^2)

    # 1. Calculate the expectation value <Pz>
    psi_dagger = psi.conj().T
    exp_Pz = (psi_dagger @ Pz @ psi).item()

    # 2. Calculate the operator Pz^2
    Pz_squared = Pz @ Pz

    # 3. Calculate the expectation value <Pz^2>
    exp_Pz_squared = (psi_dagger @ Pz_squared @ psi).item()

    # 4. Calculate the variance and uncertainty
    variance_Pz = exp_Pz_squared - exp_Pz**2
    # Ensure variance is not negative due to floating point errors
    if variance_Pz < 0 and np.isclose(variance_Pz, 0):
        variance_Pz = 0
    uncertainty_Pz = np.sqrt(variance_Pz)

    # --- Compare with the given answer ---
    # The answer A corresponds to ħ/sqrt(2).
    # With hbar=1, this is 1/sqrt(2).
    expected_uncertainty = hbar / np.sqrt(2)

    if np.isclose(uncertainty_Pz, expected_uncertainty):
        return "Correct"
    else:
        # Provide detailed calculation results if the final answer is wrong.
        return (f"Incorrect: The final calculated uncertainty does not match the answer.\n"
                f"Calculated <Pz> = {exp_Pz:.4f} ħ\n"
                f"Calculated <Pz^2> = {exp_Pz_squared:.4f} ħ^2\n"
                f"Calculated Variance (<Pz^2> - <Pz>^2) = {variance_Pz:.4f} ħ^2\n"
                f"Calculated Uncertainty (sqrt(Variance)) = {uncertainty_Pz:.4f} ħ\n"
                f"Expected Uncertainty from answer A = {expected_uncertainty:.4f} ħ")

# Run the check and print the result.
result = check_correctness()
print(result)