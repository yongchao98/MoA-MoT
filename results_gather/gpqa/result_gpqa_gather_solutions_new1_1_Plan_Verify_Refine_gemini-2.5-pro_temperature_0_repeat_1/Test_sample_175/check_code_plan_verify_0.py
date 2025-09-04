import numpy as np
import math

def check_quantum_measurement_probability():
    """
    This function checks the correctness of the calculated probability for a sequence of quantum measurements.
    
    It follows these steps:
    1. Defines the initial state and operators from the problem description.
    2. Normalizes the initial state vector.
    3. Calculates the probability of measuring P=0 by finding the corresponding eigenvector and projecting the state.
    4. Determines the collapsed state after the first measurement.
    5. Calculates the probability of measuring Q=-1 in the collapsed state.
    6. Computes the final joint probability by multiplying the two probabilities.
    7. Compares the result with the expected answer (1/6).
    """
    
    # 1. Define initial state and operators
    psi = np.array([-1, 2, 1], dtype=float)
    P = np.array([
        [0, 1/math.sqrt(2), 0],
        [1/math.sqrt(2), 0, 1/math.sqrt(2)],
        [0, 1/math.sqrt(2), 0]
    ], dtype=float)
    Q = np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ], dtype=float)
    
    # Expected final answer from the provided solution
    expected_prob = 1/6

    # 2. Normalize the initial state
    norm_psi = np.linalg.norm(psi)
    if not np.isclose(norm_psi, math.sqrt(6)):
        return f"Incorrect normalization constant. Expected sqrt(6), got {norm_psi}."
    psi_norm = psi / norm_psi

    # 3. First Measurement: Probability of P=0
    # Find eigenvalues and eigenvectors of P
    p_eigenvalues, p_eigenvectors = np.linalg.eig(P)
    
    # Find the eigenvector corresponding to eigenvalue 0
    try:
        zero_idx = np.where(np.isclose(p_eigenvalues, 0))[0][0]
        v_p0 = p_eigenvectors[:, zero_idx]
    except IndexError:
        return "Could not find an eigenvector for P with eigenvalue 0."

    # Calculate the probability of measuring P=0
    # Prob(P=0) = |<v_p0|psi_norm>|^2
    inner_product_p = np.dot(v_p0.conj(), psi_norm)
    prob_p0 = np.abs(inner_product_p)**2
    
    if not np.isclose(prob_p0, 1/3):
        return f"Incorrect probability for P=0. Expected 1/3, but calculated {prob_p0:.4f}."

    # 4. State Collapse
    # The state collapses to the eigenvector v_p0
    psi_prime = v_p0

    # 5. Second Measurement: Probability of Q=-1
    # Find eigenvalues and eigenvectors of Q
    q_eigenvalues, q_eigenvectors = np.linalg.eig(Q)

    # Find the eigenvector corresponding to eigenvalue -1
    try:
        neg_one_idx = np.where(np.isclose(q_eigenvalues, -1))[0][0]
        v_q_neg1 = q_eigenvectors[:, neg_one_idx]
    except IndexError:
        return "Could not find an eigenvector for Q with eigenvalue -1."

    # Calculate the probability of measuring Q=-1 in the collapsed state psi_prime
    # Prob(Q=-1 | P=0) = |<v_q_neg1|psi_prime>|^2
    inner_product_q = np.dot(v_q_neg1.conj(), psi_prime)
    prob_q_neg1 = np.abs(inner_product_q)**2

    if not np.isclose(prob_q_neg1, 1/2):
        return f"Incorrect probability for Q=-1 after P=0. Expected 1/2, but calculated {prob_q_neg1:.4f}."

    # 6. Calculate the final joint probability
    total_prob = prob_p0 * prob_q_neg1

    # 7. Check against the expected answer
    if np.isclose(total_prob, expected_prob):
        return "Correct"
    else:
        return f"Incorrect final probability. The calculated probability is {total_prob:.6f}, which is {total_prob}. The expected probability is 1/6 ≈ {1/6:.6f}."

# Run the check
result = check_quantum_measurement_probability()
print(result)