import numpy as np

def check_quantum_measurement_probability():
    """
    This function checks the correctness of the given quantum mechanics problem.
    It calculates the joint probability of measuring P=0 and then Q=-1.
    """
    # Given values from the problem
    # Initial state vector (unnormalized)
    psi = np.array([-1, 2, 1], dtype=complex)

    # Operator P
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=complex)

    # Operator Q
    Q = np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ], dtype=complex)

    # --- Step 1: Normalize the initial state ---
    norm_psi = np.linalg.norm(psi)
    if norm_psi == 0:
        return "Error: The initial state vector has zero norm."
    psi_norm = psi / norm_psi

    # --- Step 2: First measurement (P=0) ---
    # Find eigenvalues and eigenvectors of P
    eigenvalues_P, eigenvectors_P = np.linalg.eig(P)

    # Find the eigenvector corresponding to eigenvalue 0
    target_eigenvalue_P = 0
    # Use a tolerance for floating point comparison
    indices_p0 = np.where(np.isclose(eigenvalues_P, target_eigenvalue_P))[0]
    if len(indices_p0) == 0:
        return f"Constraint not satisfied: Operator P does not have an eigenvalue of {target_eigenvalue_P}."
    
    # The state collapses to the eigenvector corresponding to the measurement
    eigenvector_p0 = eigenvectors_P[:, indices_p0[0]]
    
    # Calculate the probability of measuring P=0
    # Prob(P=0) = |<eigenvector_p0 | psi_norm>|^2
    prob_p0 = np.abs(np.vdot(eigenvector_p0, psi_norm))**2
    
    # The state after measurement collapses to the eigenvector
    state_after_P = eigenvector_p0

    # --- Step 3: Second measurement (Q=-1) ---
    # Find eigenvalues and eigenvectors of Q
    eigenvalues_Q, eigenvectors_Q = np.linalg.eig(Q)

    # Find the eigenvector corresponding to eigenvalue -1
    target_eigenvalue_Q = -1
    indices_q_neg1 = np.where(np.isclose(eigenvalues_Q, target_eigenvalue_Q))[0]
    if len(indices_q_neg1) == 0:
        return f"Constraint not satisfied: Operator Q does not have an eigenvalue of {target_eigenvalue_Q}."
    
    eigenvector_q_neg1 = eigenvectors_Q[:, indices_q_neg1[0]]

    # Calculate the probability of measuring Q=-1 on the collapsed state
    # Prob(Q=-1 | P=0) = |<eigenvector_q_neg1 | state_after_P>|^2
    prob_q_neg1 = np.abs(np.vdot(eigenvector_q_neg1, state_after_P))**2

    # --- Step 4: Calculate the joint probability ---
    total_probability = prob_p0 * prob_q_neg1

    # --- Step 5: Check against the expected answer ---
    expected_answer = 1/6
    
    # Check if the calculated probability is close to the expected answer
    if np.isclose(total_probability, expected_answer):
        return "Correct"
    else:
        reason = (f"The answer is incorrect.\n"
                  f"The calculated total probability is {total_probability:.5f}, but the expected answer is 1/6 ≈ {expected_answer:.5f}.\n"
                  f"Breakdown of calculation:\n"
                  f"1. Probability of measuring P=0: {prob_p0:.5f} (Expected: 1/3 ≈ 0.33333)\n"
                  f"2. Probability of measuring Q=-1 after P=0: {prob_q_neg1:.5f} (Expected: 1/2 = 0.50000)\n"
                  f"3. Joint Probability = {prob_p0:.5f} * {prob_q_neg1:.5f} = {total_probability:.5f}")
        return reason

# Run the check
result = check_quantum_measurement_probability()
print(result)