import numpy as np

def check_correctness():
    """
    Checks the correctness of the answer to the quantum mechanics problem.

    The procedure is as follows:
    1. Define the initial state and operators from the problem description.
    2. Normalize the initial state vector.
    3. Calculate the probability of the first measurement (P=0) by projecting the
       normalized state onto the corresponding eigenvector of P.
    4. Determine the new state of the system after it "collapses" to this eigenvector.
    5. Calculate the probability of the second measurement (Q=-1) by projecting the
       new, collapsed state onto the corresponding eigenvector of Q.
    6. Multiply the two probabilities to get the final joint probability.
    7. Compare the result with the expected answer (1/6).
    """
    # --- Define problem constants ---
    # Initial state vector (unnormalized)
    psi_initial = np.array([-1, 2, 1], dtype=complex)

    # Operator P
    P = np.array([[0, 1/np.sqrt(2), 0],
                  [1/np.sqrt(2), 0, 1/np.sqrt(2)],
                  [0, 1/np.sqrt(2), 0]], dtype=complex)

    # Operator Q
    Q = np.array([[1, 0, 0],
                  [0, 0, 0],
                  [0, 0, -1]], dtype=complex)

    # The expected numerical result from the provided answer
    expected_prob_val = 1/6
    
    # --- Step 1: Normalize the initial state ---
    norm_psi = np.linalg.norm(psi_initial)
    # The squared norm should be (-1)^2 + 2^2 + 1^2 = 6
    if not np.isclose(norm_psi**2, 6):
        return f"Constraint check failed: The squared norm of the initial state should be 6, but was calculated as {norm_psi**2}."
    psi_norm = psi_initial / norm_psi

    # --- Step 2: First Measurement (P=0) ---
    # Find eigenvalues and eigenvectors of P. eigh is used for Hermitian matrices.
    eigenvalues_P, eigenvectors_P = np.linalg.eigh(P)
    
    # Find the eigenvector for eigenvalue 0
    try:
        idx_P0 = np.where(np.isclose(eigenvalues_P, 0))[0][0]
        eigenvector_P0 = eigenvectors_P[:, idx_P0]
    except IndexError:
        return f"Constraint check failed: Could not find an eigenvector for P with eigenvalue 0. Found eigenvalues: {eigenvalues_P}"

    # Calculate the probability of measuring P=0
    # Prob(P=0) = |<eigenvector_P0|psi_norm>|^2
    prob_P0 = np.abs(np.dot(eigenvector_P0.conj(), psi_norm))**2
    
    if not np.isclose(prob_P0, 1/3):
        return f"Intermediate check failed: The probability of measuring P=0 should be 1/3, but was calculated as {prob_P0:.4f}."

    # --- Step 3: State Collapse ---
    # The new state is the eigenvector corresponding to the measurement outcome
    psi_collapsed = eigenvector_P0

    # --- Step 4: Second Measurement (Q=-1) ---
    # Find eigenvalues and eigenvectors of Q
    eigenvalues_Q, eigenvectors_Q = np.linalg.eigh(Q)

    # Find the eigenvector for eigenvalue -1
    try:
        idx_Q_neg1 = np.where(np.isclose(eigenvalues_Q, -1))[0][0]
        eigenvector_Q_neg1 = eigenvectors_Q[:, idx_Q_neg1]
    except IndexError:
        return f"Constraint check failed: Could not find an eigenvector for Q with eigenvalue -1. Found eigenvalues: {eigenvalues_Q}"

    # Calculate the conditional probability of measuring Q=-1
    # Prob(Q=-1|P=0) = |<eigenvector_Q_neg1|psi_collapsed>|^2
    prob_Q_neg1_given_P0 = np.abs(np.dot(eigenvector_Q_neg1.conj(), psi_collapsed))**2

    if not np.isclose(prob_Q_neg1_given_P0, 1/2):
        return f"Intermediate check failed: The probability of measuring Q=-1 after P=0 should be 1/2, but was calculated as {prob_Q_neg1_given_P0:.4f}."

    # --- Step 5: Combined Probability ---
    total_prob = prob_P0 * prob_Q_neg1_given_P0

    # --- Step 6: Final Check ---
    # Check if the calculated total probability matches the expected value
    if np.isclose(total_prob, expected_prob_val):
        # The numerical calculation is correct.
        # The question lists the options as: A) 1/2, B) 1/3, C) 1/6, D) 2/3.
        # The provided answer is <<<C>>>, which corresponds to 1/6.
        # Since the calculated value is 1/6, the answer is correct.
        return "Correct"
    else:
        return (f"Incorrect: The final calculated probability is {total_prob:.4f}, "
                f"which does not match the expected value of {expected_prob_val:.4f} (1/6).")

# Run the check
result = check_correctness()
print(result)