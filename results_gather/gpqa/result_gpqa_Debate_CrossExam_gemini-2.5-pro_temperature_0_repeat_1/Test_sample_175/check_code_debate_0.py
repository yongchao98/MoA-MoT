import numpy as np

def check_correctness():
    """
    This function verifies the step-by-step calculation for the quantum measurement problem.
    """
    # --- 1. Define initial state and operators ---
    # Initial state vector |ψ⟩
    psi = np.array([-1, 2, 1], dtype=float)

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

    # --- 2. Normalize the initial state ---
    norm_psi = np.linalg.norm(psi)
    psi_norm = psi / norm_psi
    
    # Verify the squared norm is 6, as calculated in the solution
    if not np.isclose(norm_psi**2, 6):
        return f"Constraint check failed: The squared norm of the initial state is {norm_psi**2}, not 6."

    # --- 3. Find the eigenvector of P for eigenvalue 0 ---
    # Since P is a real symmetric matrix, we can use np.linalg.eigh
    eigvals_P, eigvecs_P = np.linalg.eigh(P)
    
    # Find the index of the eigenvalue that is close to 0
    try:
        p0_index = np.where(np.isclose(eigvals_P, 0))[0][0]
    except IndexError:
        return "Constraint check failed: Operator P does not have an eigenvalue of 0."
        
    # Get the corresponding normalized eigenvector |p_0⟩
    p0_vec = eigvecs_P[:, p0_index]

    # --- 4. Calculate the probability of measuring P=0 ---
    # Prob(P=0) = |⟨p_0|ψ_norm⟩|^2
    prob_P0 = np.abs(np.dot(p0_vec.conj(), psi_norm))**2
    
    # Verify this intermediate probability
    if not np.isclose(prob_P0, 1/3):
        return f"Calculation error: The probability of measuring P=0 is {prob_P0:.4f}, not 1/3."

    # --- 5. State collapse ---
    # After the measurement, the new state |ψ'⟩ is the eigenvector |p_0⟩
    psi_collapsed = p0_vec

    # --- 6. Find the eigenvector of Q for eigenvalue -1 ---
    eigvals_Q, eigvecs_Q = np.linalg.eigh(Q)
    
    # Find the index of the eigenvalue that is close to -1
    try:
        q_minus1_index = np.where(np.isclose(eigvals_Q, -1))[0][0]
    except IndexError:
        return "Constraint check failed: Operator Q does not have an eigenvalue of -1."
        
    # Get the corresponding normalized eigenvector |q_-1⟩
    q_minus1_vec = eigvecs_Q[:, q_minus1_index]

    # --- 7. Calculate the probability of measuring Q=-1 in the collapsed state ---
    # Prob(Q=-1 | after P=0) = |⟨q_-1|ψ'⟩|^2
    prob_Q_minus1 = np.abs(np.dot(q_minus1_vec.conj(), psi_collapsed))**2

    # Verify this conditional probability
    if not np.isclose(prob_Q_minus1, 1/2):
        return f"Calculation error: The conditional probability of measuring Q=-1 is {prob_Q_minus1:.4f}, not 1/2."

    # --- 8. Calculate the final total probability ---
    total_prob = prob_P0 * prob_Q_minus1
    
    # --- 9. Compare with the expected answer (1/6) ---
    expected_answer = 1/6
    if np.isclose(total_prob, expected_answer):
        return "Correct"
    else:
        return f"Incorrect final answer: The calculated total probability is {total_prob:.4f}, which is not {expected_answer:.4f} (1/6)."

# Run the check
result = check_correctness()
print(result)