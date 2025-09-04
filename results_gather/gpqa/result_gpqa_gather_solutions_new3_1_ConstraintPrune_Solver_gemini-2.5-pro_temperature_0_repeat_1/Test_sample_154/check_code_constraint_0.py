import numpy as np

def check_correctness_of_answer():
    """
    This function checks the correctness of the final answer by performing the quantum mechanical calculation.
    """
    # For simplicity in calculation, let hbar = 1. The final result will be in units of hbar.
    hbar = 1.0

    # Define the operator P_z from the problem description.
    Pz = hbar * np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ], dtype=float)

    # Define the state vector |psi> from the problem description.
    psi = np.array([
        [-1/2],
        [1/np.sqrt(2)],
        [-1/2]
    ], dtype=float)

    # The corresponding bra vector <psi| is the conjugate transpose of |psi>.
    psi_bra = psi.conj().T

    # Constraint 1: The state vector must be normalized, i.e., <psi|psi> = 1.
    norm_squared = (psi_bra @ psi)[0, 0]
    if not np.isclose(norm_squared, 1.0):
        return f"Incorrect: The provided state vector is not normalized. The inner product <psi|psi> is {norm_squared}, but it should be 1."

    # Constraint 2: Calculate the expectation value of P_z, <P_z>.
    # <P_z> = <psi|P_z|psi>
    exp_Pz = (psi_bra @ Pz @ psi)[0, 0]
    
    # The analysis correctly states that <P_z> should be 0.
    if not np.isclose(exp_Pz, 0.0):
        return f"Incorrect: The calculated expectation value <P_z> is {exp_Pz}, which contradicts the analysis result of 0."

    # Constraint 3: Calculate the expectation value of P_z^2, <P_z^2>.
    # First, find the matrix for the operator P_z^2.
    Pz_squared_matrix = Pz @ Pz
    # <P_z^2> = <psi|P_z^2|psi>
    exp_Pz_squared = (psi_bra @ Pz_squared_matrix @ psi)[0, 0]
    
    # The analysis correctly states that <P_z^2> should be hbar^2 / 2.
    expected_exp_Pz_squared = hbar**2 / 2
    if not np.isclose(exp_Pz_squared, expected_exp_Pz_squared):
        return f"Incorrect: The calculated expectation value <P_z^2> is {exp_Pz_squared}, which contradicts the analysis result of {expected_exp_Pz_squared} (hbar^2/2)."

    # Constraint 4: Calculate the uncertainty Delta P_z.
    # Delta P_z = sqrt(<P_z^2> - <P_z>^2)
    variance = exp_Pz_squared - exp_Pz**2
    delta_Pz = np.sqrt(variance)

    # The final answer is C, which corresponds to hbar / sqrt(2).
    # We compare our calculated uncertainty with this value.
    expected_delta_Pz = hbar / np.sqrt(2)
    
    if not np.isclose(delta_Pz, expected_delta_Pz):
        return f"Incorrect: The final calculated uncertainty Delta P_z is {delta_Pz}*hbar. The correct answer C corresponds to {expected_delta_Pz}*hbar (hbar/sqrt(2)). The final answer is therefore incorrect."

    # If all checks pass, the answer is correct.
    return "Correct"

# print(check_correctness_of_answer())