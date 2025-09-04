import numpy as np

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the quantum mechanics problem.

    The problem asks for the uncertainty Delta P_z for a system in a given state.
    The uncertainty is calculated using the formula: Delta P_z = sqrt(<P_z^2> - <P_z>^2).

    The function will:
    1. Define the operators and the state vector from the problem description.
    2. Verify that the given state satisfies the problem's constraints (i.e., it's a normalized
       eigenstate of P_x with eigenvalue -hbar).
    3. Re-calculate the expectation values <P_z> and <P_z^2>.
    4. Compute the uncertainty Delta P_z.
    5. Compare the computed result with the LLM's answer (hbar/sqrt(2)).
    """
    # For numerical calculations, we can set hbar = 1 and check the coefficient.
    # The final answer will be in units of hbar.
    hbar = 1.0

    # Define the operator P_x
    P_x = (hbar / np.sqrt(2)) * np.array([
        [0, 1, 0],
        [1, 0, 1],
        [0, 1, 0]
    ])

    # Define the operator P_z
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

    # --- Constraint Verification ---

    # 1. Verify that the state is normalized, i.e., <psi|psi> = 1
    # <psi| is the conjugate transpose of |psi>
    psi_bra = psi.conj().T
    norm_squared = (psi_bra @ psi)[0, 0]
    if not np.isclose(norm_squared, 1.0):
        return f"Constraint not satisfied: The state vector is not normalized. <psi|psi> = {norm_squared:.4f}, but it should be 1."

    # 2. Verify that the state is an eigenstate of P_x with eigenvalue -hbar
    eigenvalue = -hbar
    # Calculate P_x |psi>
    Px_psi = P_x @ psi
    # Calculate (-hbar) * |psi>
    eigenvalue_psi = eigenvalue * psi
    if not np.allclose(Px_psi, eigenvalue_psi):
        return (f"Constraint not satisfied: The given state is not an eigenstate of P_x with eigenvalue -hbar.\n"
                f"P_x|psi> results in:\n{Px_psi}\n"
                f"But (-hbar)|psi> is:\n{eigenvalue_psi}")

    # --- Calculation of Uncertainty ---

    # Calculate the expectation value of P_z: <P_z> = <psi|P_z|psi>
    exp_Pz = (psi_bra @ P_z @ psi)[0, 0]

    # Calculate the operator P_z^2
    Pz_squared = P_z @ P_z

    # Calculate the expectation value of P_z^2: <P_z^2> = <psi|P_z^2|psi>
    exp_Pz_squared = (psi_bra @ Pz_squared @ psi)[0, 0]

    # Calculate the variance (Delta P_z)^2
    variance_Pz = exp_Pz_squared - exp_Pz**2
    
    # The variance must be non-negative
    if variance_Pz < 0:
        return f"Calculation error: The variance is negative ({variance_Pz:.4f}), which is physically impossible."

    # Calculate the uncertainty (standard deviation)
    delta_Pz = np.sqrt(variance_Pz)

    # The expected answer from the LLM is hbar/sqrt(2)
    expected_answer = hbar / np.sqrt(2)

    # Compare the calculated result with the expected answer
    if np.isclose(delta_Pz, expected_answer):
        return "Correct"
    else:
        # The LLM's reasoning and steps were correct, so we check the final value.
        llm_answer_val = hbar / np.sqrt(2)
        return (f"Incorrect final answer. The calculated uncertainty is Delta P_z = {delta_Pz:.4f} * hbar. "
                f"The LLM's answer is {llm_answer_val:.4f} * hbar. The values do not match.")

# Execute the check and print the result
result = check_correctness()
# The LLM's answer is correct, so the code should return "Correct".
# If it were incorrect, it would return a reason.
# For the final output, we just need the code block.
# print(result)