import numpy as np

def check_quantum_measurement_probability():
    """
    This function verifies the calculation for a quantum measurement problem.
    It checks each step of the provided solution to ensure its correctness.
    """
    try:
        # Define the initial state and operators as given in the problem
        initial_state = np.array([[-1], [2], [1]], dtype=complex)
        
        P_operator = np.array([
            [0, 1/np.sqrt(2), 0],
            [1/np.sqrt(2), 0, 1/np.sqrt(2)],
            [0, 1/np.sqrt(2), 0]
        ], dtype=float)
        
        Q_operator = np.array([
            [1, 0, 0],
            [0, 0, 0],
            [0, 0, -1]
        ], dtype=float)

        # --- Step 1: Normalize the initial state ---
        norm = np.linalg.norm(initial_state)
        if not np.isclose(norm, np.sqrt(6)):
            return f"Constraint failed: The norm of the initial state should be sqrt(6), but was calculated as {norm}."
        
        normalized_state = initial_state / norm

        # --- Step 2: First Measurement (P=0) ---
        # The solution claims the eigenvector for P=0 is (1/sqrt(2), 0, -1/sqrt(2))
        eigenvalue_P = 0
        eigenvector_P_0 = np.array([[1/np.sqrt(2)], [0], [-1/np.sqrt(2)]], dtype=complex)

        # Verify it's a correct eigenvector
        if not np.allclose(P_operator @ eigenvector_P_0, eigenvalue_P * eigenvector_P_0):
            return f"Constraint failed: The vector (1/sqrt(2), 0, -1/sqrt(2)) is not an eigenvector of P for eigenvalue 0."

        # Calculate the probability of measuring P=0
        # Prob(P=0) = |<eigenvector_P_0 | normalized_state>|^2
        # np.vdot performs the conjugate transpose and dot product
        prob_P_0 = np.abs(np.vdot(eigenvector_P_0, normalized_state))**2
        
        expected_prob_P_0 = 1/3
        if not np.isclose(prob_P_0, expected_prob_P_0):
            return f"Calculation mismatch: The probability of measuring P=0 is {prob_P_0:.4f}, but the solution states it is {expected_prob_P_0:.4f} (1/3)."

        # --- Step 3: State Collapse ---
        # After measuring P=0, the state collapses to the corresponding eigenvector
        state_after_P = eigenvector_P_0

        # --- Step 4: Second Measurement (Q=-1) ---
        # The solution claims the eigenvector for Q=-1 is (0, 0, 1)
        eigenvalue_Q = -1
        eigenvector_Q_neg1 = np.array([[0], [0], [1]], dtype=complex)

        # Verify it's a correct eigenvector
        if not np.allclose(Q_operator @ eigenvector_Q_neg1, eigenvalue_Q * eigenvector_Q_neg1):
            return f"Constraint failed: The vector (0, 0, 1) is not an eigenvector of Q for eigenvalue -1."

        # Calculate the probability of measuring Q=-1 in the new state
        # Prob(Q=-1 after P=0) = |<eigenvector_Q_-1 | state_after_P>|^2
        prob_Q_neg1_after_P = np.abs(np.vdot(eigenvector_Q_neg1, state_after_P))**2
        
        expected_prob_Q_neg1_after_P = 1/2
        if not np.isclose(prob_Q_neg1_after_P, expected_prob_Q_neg1_after_P):
            return f"Calculation mismatch: The probability of measuring Q=-1 after P=0 is {prob_Q_neg1_after_P:.4f}, but the solution states it is {expected_prob_Q_neg1_after_P:.4f} (1/2)."

        # --- Step 5: Joint Probability ---
        joint_probability = prob_P_0 * prob_Q_neg1_after_P
        
        expected_final_answer = 1/6
        if not np.isclose(joint_probability, expected_final_answer):
            return f"Final Answer Incorrect: The calculated joint probability is {joint_probability:.4f}, which does not match the expected answer {expected_final_answer:.4f} (1/6)."

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_quantum_measurement_probability()
print(result)