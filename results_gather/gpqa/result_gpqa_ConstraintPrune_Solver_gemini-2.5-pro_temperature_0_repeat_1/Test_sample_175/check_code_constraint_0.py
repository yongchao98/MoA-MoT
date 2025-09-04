import numpy as np

def check_correctness():
    """
    This function verifies the solution to the quantum mechanics problem.
    It calculates the joint probability of two sequential measurements and
    compares it to the expected answer.
    """
    # --- 1. Define the initial state and operators ---
    # The initial state of the system (unnormalized)
    psi_initial_unnormalized = np.array([-1, 2, 1], dtype=complex).reshape(3, 1)

    # The operator P
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=complex)

    # The operator Q
    Q = np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ], dtype=complex)

    # The expected answer from option C
    expected_probability = 1/6

    # --- 2. Normalize the initial state ---
    norm = np.linalg.norm(psi_initial_unnormalized)
    if np.isclose(norm, 0):
        return "Error: The initial state vector is a zero vector and cannot be normalized."
    psi_normalized = psi_initial_unnormalized / norm

    # --- 3. First Measurement: P=0 ---
    # Find eigenvalues and eigenvectors of P
    eigenvalues_P, eigenvectors_P = np.linalg.eig(P)

    # Find the eigenvector for the eigenvalue 0
    target_eigenvalue_P = 0
    try:
        idx_P = np.where(np.isclose(eigenvalues_P, target_eigenvalue_P))[0][0]
        eigenvector_P0 = eigenvectors_P[:, idx_P].reshape(3, 1)
    except IndexError:
        return f"Constraint failed: Operator P does not have an eigenvalue of {target_eigenvalue_P}."

    # Calculate the probability of measuring P=0
    # Prob(P=0) = |<eigenvector_P0 | psi_normalized>|^2
    prob_P0 = np.abs(np.vdot(eigenvector_P0, psi_normalized))**2

    # The state collapses to the eigenvector of P for eigenvalue 0
    state_after_P = eigenvector_P0

    # --- 4. Second Measurement: Q=-1 ---
    # Find eigenvalues and eigenvectors of Q
    eigenvalues_Q, eigenvectors_Q = np.linalg.eig(Q)

    # Find the eigenvector for the eigenvalue -1
    target_eigenvalue_Q = -1
    try:
        idx_Q = np.where(np.isclose(eigenvalues_Q, target_eigenvalue_Q))[0][0]
        eigenvector_Q_neg1 = eigenvectors_Q[:, idx_Q].reshape(3, 1)
    except IndexError:
        return f"Constraint failed: Operator Q does not have an eigenvalue of {target_eigenvalue_Q}."

    # Calculate the probability of measuring Q=-1 on the new state
    # Prob(Q=-1 after P=0) = |<eigenvector_Q_neg1 | state_after_P>|^2
    prob_Q_neg1_after_P0 = np.abs(np.vdot(eigenvector_Q_neg1, state_after_P))**2

    # --- 5. Calculate the final joint probability ---
    final_probability = prob_P0 * prob_Q_neg1_after_P0

    # --- 6. Verify the result ---
    # Check if the calculated probability matches the expected answer (1/6)
    if not np.isclose(final_probability, expected_probability):
        # Check intermediate steps to provide a more detailed reason
        if not np.isclose(prob_P0, 1/3):
            return (f"Incorrect. The probability of the first measurement (P=0) is calculated as {prob_P0:.4f}, "
                    f"but the correct value is 1/3. This leads to an incorrect final answer.")
        if not np.isclose(prob_Q_neg1_after_P0, 1/2):
            return (f"Incorrect. The probability of the second measurement (Q=-1) is calculated as {prob_Q_neg1_after_P0:.4f}, "
                    f"but the correct value is 1/2. This leads to an incorrect final answer.")
        return (f"Incorrect. The final calculated probability is {final_probability:.4f}, "
                f"which does not match the expected answer of 1/6 ({expected_probability:.4f}).")

    return "Correct"

# Run the check
result = check_correctness()
print(result)