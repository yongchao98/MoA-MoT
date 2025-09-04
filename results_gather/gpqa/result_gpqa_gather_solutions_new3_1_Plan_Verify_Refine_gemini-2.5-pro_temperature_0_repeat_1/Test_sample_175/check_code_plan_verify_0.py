import numpy as np

def check_quantum_measurement_probability():
    """
    This function verifies the calculation for a sequential quantum measurement problem.
    It calculates the probability of measuring P=0 and then Q=-1.
    """
    try:
        # Step 1: Define the initial state and operators
        # The initial state vector (unnormalized)
        psi = np.array([[-1], [2], [1]], dtype=complex)

        # The matrix for operator P
        P = np.array([
            [0, 1/np.sqrt(2), 0],
            [1/np.sqrt(2), 0, 1/np.sqrt(2)],
            [0, 1/np.sqrt(2), 0]
        ], dtype=complex)

        # The matrix for operator Q
        Q = np.array([
            [1, 0, 0],
            [0, 0, 0],
            [0, 0, -1]
        ], dtype=complex)

        # Step 2: Normalize the initial state
        norm_psi = np.linalg.norm(psi)
        if np.isclose(norm_psi, 0):
            return "Incorrect: The initial state vector has zero norm."
        psi_norm = psi / norm_psi

        # Step 3: Find the eigenvector of P for eigenvalue 0
        eigenvalues_P, eigenvectors_P = np.linalg.eig(P)
        # Find the index corresponding to the eigenvalue 0
        p0_indices = np.where(np.isclose(eigenvalues_P, 0))[0]
        if len(p0_indices) == 0:
            return "Incorrect: Operator P does not have an eigenvalue of 0."
        # Get the normalized eigenvector for P=0
        p0_vec = eigenvectors_P[:, p0_indices[0]].reshape(3, 1)

        # Step 4: Calculate the probability of measuring P=0
        # Prob(P=0) = |<p0|psi_norm>|^2
        # np.vdot is the complex-conjugate dot product
        amplitude1 = np.vdot(p0_vec, psi_norm)
        prob_P0 = np.abs(amplitude1)**2

        # Step 5: The state collapses to the eigenvector p0_vec
        psi_prime = p0_vec

        # Step 6: Find the eigenvector of Q for eigenvalue -1
        eigenvalues_Q, eigenvectors_Q = np.linalg.eig(Q)
        # Find the index corresponding to the eigenvalue -1
        q_neg1_indices = np.where(np.isclose(eigenvalues_Q, -1))[0]
        if len(q_neg1_indices) == 0:
            return "Incorrect: Operator Q does not have an eigenvalue of -1."
        # Get the normalized eigenvector for Q=-1
        q_neg1_vec = eigenvectors_Q[:, q_neg1_indices[0]].reshape(3, 1)

        # Step 7: Calculate the probability of measuring Q=-1 from the collapsed state
        # Prob(Q=-1 | P=0) = |<q_neg1|psi_prime>|^2
        amplitude2 = np.vdot(q_neg1_vec, psi_prime)
        prob_Q_neg1_given_P0 = np.abs(amplitude2)**2

        # Step 8: Calculate the total joint probability
        total_prob = prob_P0 * prob_Q_neg1_given_P0

        # Step 9: Compare with the expected answer (1/6)
        expected_prob = 1/6
        if np.isclose(total_prob, expected_prob):
            return "Correct"
        else:
            # Provide a detailed reason for the failure
            reason = (
                f"The final calculated probability is {total_prob:.6f}, which is not the expected {expected_prob:.6f}.\n"
                f"Details:\n"
                f"  - Probability of P=0: {prob_P0:.6f} (Expected: {1/3:.6f})\n"
                f"  - Probability of Q=-1 after P=0: {prob_Q_neg1_given_P0:.6f} (Expected: {1/2:.6f})"
            )
            return f"Incorrect: {reason}"

    except Exception as e:
        return f"An error occurred during the calculation: {e}"

# Execute the check
result = check_quantum_measurement_probability()
print(result)