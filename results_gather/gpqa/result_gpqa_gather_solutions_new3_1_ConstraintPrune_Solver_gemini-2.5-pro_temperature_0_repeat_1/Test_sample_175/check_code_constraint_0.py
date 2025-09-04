import numpy as np

def check_quantum_measurement_probability():
    """
    This function verifies the calculation for a sequential quantum measurement problem.
    It calculates the probability of measuring P=0 and then Q=-1.
    """
    
    # --- Step 1: Define initial state and operators ---
    # Initial state vector |ψ⟩
    psi = np.array([[-1], [2], [1]], dtype=float)
    
    # Operator P
    P = np.array([[0, 1/np.sqrt(2), 0],
                  [1/np.sqrt(2), 0, 1/np.sqrt(2)],
                  [0, 1/np.sqrt(2), 0]], dtype=float)
                  
    # Operator Q
    Q = np.array([[1, 0, 0],
                  [0, 0, 0],
                  [0, 0, -1]], dtype=float)

    # --- Step 2: Normalize the initial state ---
    norm_psi = np.linalg.norm(psi)
    psi_norm = psi / norm_psi
    
    # --- Step 3: First Measurement (P=0) ---
    
    # Find eigenvalues and eigenvectors of P. np.linalg.eigh is used for Hermitian matrices.
    eigenvalues_P, eigenvectors_P = np.linalg.eigh(P)
    
    # Find the index of the eigenvalue that is close to 0.
    try:
        p0_index = np.where(np.isclose(eigenvalues_P, 0))[0][0]
    except IndexError:
        return "Constraint not satisfied: Could not find eigenvalue 0 for operator P."
        
    # Get the corresponding normalized eigenvector |p₀⟩.
    p0_vec = eigenvectors_P[:, p0_index].reshape(3, 1)

    # Calculate the probability of measuring P=0: Prob(P=0) = |⟨p₀|ψ_norm⟩|²
    # np.vdot calculates the inner product.
    inner_product_P = np.vdot(p0_vec, psi_norm)
    prob_P0 = np.abs(inner_product_P)**2

    # Verify this intermediate probability. It should be 1/3.
    if not np.isclose(prob_P0, 1/3):
        return f"Constraint not satisfied: The probability of measuring P=0 was calculated as {prob_P0:.4f}, but it should be 1/3."

    # --- Step 4: State Collapse ---
    # The new state |ψ'⟩ is the eigenvector corresponding to the measurement outcome.
    psi_prime = p0_vec

    # --- Step 5: Second Measurement (Q=-1) ---
    
    # Find eigenvalues and eigenvectors of Q.
    eigenvalues_Q, eigenvectors_Q = np.linalg.eigh(Q)
    
    # Find the index of the eigenvalue -1.
    try:
        q_neg1_index = np.where(np.isclose(eigenvalues_Q, -1))[0][0]
    except IndexError:
        return "Constraint not satisfied: Could not find eigenvalue -1 for operator Q."
        
    # Get the corresponding normalized eigenvector |q₋₁⟩.
    q_neg1_vec = eigenvectors_Q[:, q_neg1_index].reshape(3, 1)

    # Calculate the probability of measuring Q=-1 in the new state: Prob(Q=-1|P=0) = |⟨q₋₁|ψ'⟩|²
    inner_product_Q = np.vdot(q_neg1_vec, psi_prime)
    prob_Q_neg1_given_P0 = np.abs(inner_product_Q)**2

    # Verify this intermediate probability. It should be 1/2.
    if not np.isclose(prob_Q_neg1_given_P0, 1/2):
        return f"Constraint not satisfied: The conditional probability of measuring Q=-1 was calculated as {prob_Q_neg1_given_P0:.4f}, but it should be 1/2."

    # --- Step 6: Calculate the Total Probability ---
    total_prob = prob_P0 * prob_Q_neg1_given_P0

    # --- Step 7: Check against the final answer ---
    # The provided answer is C, which corresponds to 1/6.
    expected_answer_value = 1/6
    
    if np.isclose(total_prob, expected_answer_value):
        return "Correct"
    else:
        return f"Incorrect: The final calculated probability is {total_prob:.4f}, which does not match the expected answer of 1/6 ({expected_answer_value:.4f})."

# Run the check
result = check_quantum_measurement_probability()
print(result)