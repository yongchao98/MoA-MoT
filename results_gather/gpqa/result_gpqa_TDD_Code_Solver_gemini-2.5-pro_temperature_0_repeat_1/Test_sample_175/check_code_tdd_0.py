import numpy as np

def check_quantum_measurement_probability():
    """
    This function verifies the calculation for a sequential quantum measurement problem.
    It calculates the probability of measuring P=0 and then Q=-1 and compares it
    to the expected answer of 1/6.
    """
    # --- 1. Define the problem's constants and matrices ---
    sqrt2 = np.sqrt(2)
    
    # The initial state of the system, |ψ⟩
    initial_state = np.array([-1, 2, 1], dtype=complex)
    
    # The operator matrix for observable P
    P_matrix = np.array([
        [0, 1/sqrt2, 0],
        [1/sqrt2, 0, 1/sqrt2],
        [0, 1/sqrt2, 0]
    ], dtype=float)
    
    # The operator matrix for observable Q
    Q_matrix = np.array([
        [1, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ], dtype=float)
    
    # The target measurement outcomes (eigenvalues)
    p_target_value = 0
    q_target_value = -1
    
    # The expected final probability from the provided answer (Option D: 1/6)
    expected_probability = 1/6

    # --- 2. Perform the calculation step-by-step ---

    # Step A: Normalize the initial state vector
    # The probability amplitudes must be calculated with a state of norm 1.
    # ⟨ψ|ψ⟩ = (-1)*(-1) + 2*2 + 1*1 = 1 + 4 + 1 = 6
    norm_squared = np.vdot(initial_state, initial_state).real
    if norm_squared == 0:
        return "Error: The initial state vector cannot be a zero vector."
    normalized_state = initial_state / np.sqrt(norm_squared)

    # Step B: Find the eigenvector of P for the eigenvalue p=0
    p_eigenvalues, p_eigenvectors = np.linalg.eig(P_matrix)
    try:
        p_index = np.where(np.isclose(p_eigenvalues, p_target_value))[0][0]
        p0_eigenvector = p_eigenvectors[:, p_index]
    except IndexError:
        return f"Constraint not satisfied: The operator P does not have an eigenvalue of {p_target_value}."

    # Step C: Calculate the probability of measuring P=0
    # Prob(p=0) = |⟨p=0|ψ_normalized⟩|^2
    prob_p0 = np.abs(np.vdot(p0_eigenvector, normalized_state))**2

    # Step D: The state collapses to the eigenvector of P. This is the new state.
    state_after_p_measurement = p0_eigenvector

    # Step E: Find the eigenvector of Q for the eigenvalue q=-1
    q_eigenvalues, q_eigenvectors = np.linalg.eig(Q_matrix)
    try:
        q_index = np.where(np.isclose(q_eigenvalues, q_target_value))[0][0]
        q_neg1_eigenvector = q_eigenvectors[:, q_index]
    except IndexError:
        return f"Constraint not satisfied: The operator Q does not have an eigenvalue of {q_target_value}."

    # Step F: Calculate the probability of measuring Q=-1 from the collapsed state
    # Prob(q=-1 | p=0) = |⟨q=-1|state_after_p⟩|^2
    prob_q_given_p = np.abs(np.vdot(q_neg1_eigenvector, state_after_p_measurement))**2

    # Step G: The total probability is the product of the individual probabilities
    total_probability = prob_p0 * prob_q_given_p

    # --- 3. Compare the result with the expected answer ---
    if np.isclose(total_probability, expected_probability):
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability is {total_probability:.4f}, "
                f"while the expected answer is {expected_probability:.4f} (1/6). "
                f"The breakdown is: Prob(P=0) ≈ {prob_p0:.4f} (1/3) and "
                f"Prob(Q=-1|P=0) ≈ {prob_q_given_p:.4f} (1/2).")

# Run the check
result = check_quantum_measurement_probability()
print(result)