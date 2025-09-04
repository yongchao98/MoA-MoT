import numpy as np

def check_correctness():
    """
    This function checks the correctness of the provided answer to the quantum mechanics problem.
    It recalculates the required probability from the problem statement's initial conditions.
    """
    try:
        # --- Problem Definition ---
        # Initial state vector (unnormalized)
        psi_initial = np.array([-1, 2, 1], dtype=float)

        # Operator P
        P = np.array([
            [0, 1/np.sqrt(2), 0],
            [1/np.sqrt(2), 0, 1/np.sqrt(2)],
            [0, 1/np.sqrt(2), 0]
        ], dtype=float)

        # Operator Q
        Q = np.array([
            [1, 0, 0],
            [0, 0, 0],
            [0, 0, -1]
        ], dtype=float)

        # --- Step 1: First Measurement (P=0) ---

        # Normalize the initial state
        psi_norm = psi_initial / np.linalg.norm(psi_initial)

        # Find eigenvalues and eigenvectors of P. P is Hermitian (real and symmetric).
        eigenvalues_P, eigenvectors_P = np.linalg.eigh(P)

        # Find the eigenvector corresponding to the eigenvalue 0
        target_eigenvalue_P = 0
        # Use np.isclose for safe floating-point comparison
        p0_indices = np.where(np.isclose(eigenvalues_P, target_eigenvalue_P))[0]
        if len(p0_indices) == 0:
            return f"Constraint not satisfied: Operator P does not have an eigenvalue of {target_eigenvalue_P}."
        # The eigenvector is the column corresponding to the eigenvalue's index
        p0_eigenvector = eigenvectors_P[:, p0_indices[0]]

        # Calculate the probability of measuring P=0
        # Prob(P=0) = |<p0|ψ_norm>|^2
        prob_P0 = np.abs(np.dot(p0_eigenvector.conj(), psi_norm))**2

        # --- Step 2: State Collapse ---

        # After measuring P=0, the state collapses to the corresponding eigenvector
        psi_after_P = p0_eigenvector

        # --- Step 3: Second Measurement (Q=-1) ---

        # Find eigenvalues and eigenvectors of Q. Q is Hermitian.
        eigenvalues_Q, eigenvectors_Q = np.linalg.eigh(Q)

        # Find the eigenvector corresponding to the eigenvalue -1
        target_eigenvalue_Q = -1
        q_neg1_indices = np.where(np.isclose(eigenvalues_Q, target_eigenvalue_Q))[0]
        if len(q_neg1_indices) == 0:
            return f"Constraint not satisfied: Operator Q does not have an eigenvalue of {target_eigenvalue_Q}."
        q_neg1_eigenvector = eigenvectors_Q[:, q_neg1_indices[0]]

        # Calculate the probability of measuring Q=-1 on the new state |ψ_after_P>
        # Prob(Q=-1 | P=0) = |<q_-1|ψ_after_P>|^2
        prob_Q_neg1_after_P = np.abs(np.dot(q_neg1_eigenvector.conj(), psi_after_P))**2

        # --- Step 4: Total Probability ---

        # The total probability of the sequence is the product of the individual probabilities
        total_probability = prob_P0 * prob_Q_neg1_after_P

        # The expected answer is 1/6 (Option C)
        expected_value = 1/6

        # Check if the calculated probability matches the expected answer
        if np.isclose(total_probability, expected_value):
            return "Correct"
        else:
            # Provide a detailed reason for the incorrectness
            return (f"Incorrect. The calculated probability is {total_probability:.6f}, which is not close to the expected 1/6 ≈ {expected_value:.6f}. "
                    f"The probability of measuring P=0 was calculated as {prob_P0:.4f} (expected 1/3). "
                    f"The subsequent probability of measuring Q=-1 was calculated as {prob_Q_neg1_after_P:.4f} (expected 1/2).")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness()
print(result)