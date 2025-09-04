import numpy as np

def check_answer():
    """
    This function checks the correctness of the provided LLM answer.
    It defines the quantum operators and the state vector, then calculates
    the uncertainty in P_z and compares it to the given answer.
    """
    # Use hbar = 1 for numerical calculations, and add it back symbolically at the end.
    hbar = 1.0

    # Define the matrix form of the operators based on the problem description.
    # P_x matrix
    Px = (hbar / np.sqrt(2)) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    # P_y matrix
    Py = (hbar / np.sqrt(2)) * np.array([
        [0, -1j, 0],
        [1j, 0, -1j],
        [0, 1j, 0]
    ], dtype=complex)

    # P_z matrix
    Pz = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])

    # Define the state of the system as a column vector |psi>
    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ])

    # --- Constraint Check ---
    # The question states that the system is in an eigenstate of Px with eigenvalue -hbar.
    # Let's verify this constraint.
    # Calculate Px|psi>
    result_vector = Px @ psi
    # Calculate -hbar|psi>
    expected_vector = -hbar * psi
    if not np.allclose(result_vector, expected_vector):
        return (f"Constraint check failed: The given state is not an eigenstate of Px with eigenvalue -hbar.\n"
                f"Px|psi> = {result_vector.flatten()}\n"
                f"but -hbar|psi> = {expected_vector.flatten()}")

    # --- Calculation of Uncertainty ---

    # 1. Calculate the expectation value of Pz: <Pz> = <psi|Pz|psi>
    # <psi| is the conjugate transpose of |psi>
    psi_dagger = psi.conj().T
    
    # The result of the matrix multiplication is a 1x1 matrix, so we extract the element.
    exp_Pz = (psi_dagger @ Pz @ psi)[0, 0]

    # 2. Calculate the expectation value of Pz^2: <Pz^2> = <psi|Pz^2|psi>
    # First, calculate the matrix for Pz^2
    Pz_squared = Pz @ Pz
    exp_Pz_squared = (psi_dagger @ Pz_squared @ psi)[0, 0]

    # 3. Calculate the uncertainty Delta Pz = sqrt(<Pz^2> - <Pz>^2)
    # We use np.real to handle any tiny imaginary parts from floating point errors.
    radicand = np.real(exp_Pz_squared - exp_Pz**2)
    uncertainty_Pz = np.sqrt(radicand)

    # --- Comparison with the LLM's Answer ---
    
    # Values from the LLM's answer (with hbar=1)
    llm_exp_Pz = 0.0
    llm_exp_Pz_squared = hbar**2 / 2.0
    llm_uncertainty = hbar / np.sqrt(2)

    # Check if the calculated values match the LLM's intermediate and final results.
    if not np.isclose(exp_Pz, llm_exp_Pz):
        return (f"Incorrect: The expectation value of P_z, <P_z>, is incorrect.\n"
                f"Calculated value: {exp_Pz:.4f}\n"
                f"LLM's value: {llm_exp_Pz:.4f}")

    if not np.isclose(exp_Pz_squared, llm_exp_Pz_squared):
        return (f"Incorrect: The expectation value of P_z^2, <P_z^2>, is incorrect.\n"
                f"Calculated value (in units of hbar^2): {exp_Pz_squared:.4f}\n"
                f"LLM's value (in units of hbar^2): {llm_exp_Pz_squared:.4f}")

    if not np.isclose(uncertainty_Pz, llm_uncertainty):
        return (f"Incorrect: The final uncertainty, Delta P_z, is incorrect.\n"
                f"Calculated value (in units of hbar): {uncertainty_Pz:.4f}\n"
                f"LLM's value (in units of hbar): {llm_uncertainty:.4f}")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
print(check_answer())