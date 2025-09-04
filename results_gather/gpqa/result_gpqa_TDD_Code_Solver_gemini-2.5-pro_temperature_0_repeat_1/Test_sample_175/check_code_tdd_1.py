import numpy as np

def check_answer():
    """
    This function checks the correctness of the given quantum mechanics problem solution.
    It calculates the probability of a sequence of measurements and compares it to the expected answer.
    """
    try:
        # --- Problem Setup ---
        # Initial state vector |ψ⟩
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

        # Expected final probability from the answer
        expected_total_prob = 1/6

        # --- Step 1: Normalize the initial state ---
        norm_psi = np.linalg.norm(psi_initial)
        if np.isclose(norm_psi, 0):
            return "Initial state vector cannot be a zero vector."
        psi_normalized = psi_initial / norm_psi
        
        # Verification of the norm from the provided answer (√6)
        if not np.isclose(norm_psi, np.sqrt(6)):
            return f"Constraint check failed: The norm of the initial state is {norm_psi:.4f}, not sqrt(6) as stated in the explanation."

        # --- Step 2: First Measurement (P=0) ---
        # Find eigenvalues and eigenvectors of P. eigh is used for Hermitian matrices.
        eigenvalues_P, eigenvectors_P = np.linalg.eigh(P)

        # Find the eigenvector corresponding to the eigenvalue p=0
        p_val = 0
        try:
            # Find the index of the eigenvalue that is close to 0
            p_index = np.where(np.isclose(eigenvalues_P, p_val))[0][0]
        except IndexError:
            return f"Constraint check failed: The operator P does not have an eigenvalue of {p_val}."
        
        p_eigenvector = eigenvectors_P[:, p_index]

        # Calculate the probability of measuring P=0
        # This is |⟨p=0|ψ_norm⟩|²
        amplitude_p = np.vdot(p_eigenvector, psi_normalized)
        prob_p = np.abs(amplitude_p)**2

        # Verification of the first probability (1/3)
        if not np.isclose(prob_p, 1/3):
            return f"Calculation Error: The probability of measuring P=0 is {prob_p:.4f}, but the explanation calculates it as 1/3."

        # The state collapses to the eigenvector of P after measurement
        state_after_P = p_eigenvector

        # --- Step 3: Second Measurement (Q=-1) ---
        # Find eigenvalues and eigenvectors of Q
        eigenvalues_Q, eigenvectors_Q = np.linalg.eigh(Q)

        # Find the eigenvector corresponding to the eigenvalue q=-1
        q_val = -1
        try:
            q_index = np.where(np.isclose(eigenvalues_Q, q_val))[0][0]
        except IndexError:
            return f"Constraint check failed: The operator Q does not have an eigenvalue of {q_val}."
            
        q_eigenvector = eigenvectors_Q[:, q_index]

        # Calculate the probability of measuring Q=-1, given the new state
        # This is |⟨q=-1|state_after_P⟩|²
        amplitude_q = np.vdot(q_eigenvector, state_after_P)
        prob_q_given_p = np.abs(amplitude_q)**2

        # Verification of the conditional probability (1/2)
        if not np.isclose(prob_q_given_p, 1/2):
            return f"Calculation Error: The conditional probability of measuring Q=-1 is {prob_q_given_p:.4f}, but the explanation calculates it as 1/2."

        # --- Step 4: Total Probability ---
        total_prob = prob_p * prob_q_given_p

        # --- Final Verification ---
        if np.isclose(total_prob, expected_total_prob):
            return "Correct"
        else:
            return f"Incorrect final answer. The calculated total probability is {total_prob:.4f} ({total_prob.as_integer_ratio()[0]}/{total_prob.as_integer_ratio()[1]}), which does not match the expected value of 1/6."

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Run the check and print the result
result = check_answer()
print(result)