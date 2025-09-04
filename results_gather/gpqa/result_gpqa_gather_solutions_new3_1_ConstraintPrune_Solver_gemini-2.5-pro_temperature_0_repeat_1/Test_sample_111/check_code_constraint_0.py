import numpy as np

def check_correctness():
    """
    This function checks the correctness of the answer to the quantum mechanics problem.
    It follows these steps:
    1. Normalizes the initial quantum state.
    2. Defines the operator A and finds its eigenvalues and eigenvectors.
    3. Calculates the probability of measuring the system in each eigenstate.
    4. Calculates the average (expectation) value of the operator.
    5. Compares the calculated results with the values given in option A.
    """
    # For calculation purposes, we can set hbar = 1 and add it back in the final result.
    hbar = 1.0

    # --- Step 1: Normalize the initial state |alpha> ---
    # The state is proportional to (1+i)|up> + (2-i)|down>
    # In vector form: [1+1j, 2-1j]
    alpha_unnormalized = np.array([1 + 1j, 2 - 1j])

    # The squared norm is <alpha|alpha> = |1+i|^2 + |2-i|^2
    norm_squared = np.vdot(alpha_unnormalized, alpha_unnormalized).real
    
    # The norm should be sqrt(7)
    if not np.isclose(norm_squared, 7.0):
        return f"Constraint check failed: The squared norm of the state should be 7, but was calculated as {norm_squared}."
    
    alpha_normalized = alpha_unnormalized / np.sqrt(norm_squared)

    # --- Step 2: Find eigenvalues and eigenvectors of operator A ---
    # Aij = hbar/2 if i!=j, and 0 otherwise.
    # A = (hbar/2) * [[0, 1], [1, 0]]
    A = (hbar / 2) * np.array([[0, 1], [1, 0]])

    # np.linalg.eigh is used for Hermitian matrices, guaranteeing real eigenvalues
    # and orthonormal eigenvectors.
    eigenvalues, eigenvectors = np.linalg.eigh(A)

    # The eigenvalues should be +hbar/2 and -hbar/2
    expected_eigenvalues = sorted([-hbar / 2, hbar / 2])
    if not np.allclose(sorted(eigenvalues), expected_eigenvalues):
        return f"Constraint check failed: The eigenvalues should be {expected_eigenvalues}, but were calculated as {sorted(eigenvalues)}."

    # --- Step 3: Calculate the probabilities ---
    # The probability of measuring an eigenvalue is the squared magnitude of the
    # projection of the state onto the corresponding eigenvector.
    # P(lambda_i) = |<v_i|alpha>|^2
    
    # The eigenvectors are the columns of the 'eigenvectors' matrix.
    prob_1 = abs(np.vdot(eigenvectors[:, 0], alpha_normalized))**2
    prob_2 = abs(np.vdot(eigenvectors[:, 1], alpha_normalized))**2
    
    # The exact probabilities are 9/14 and 5/14.
    exact_probs = sorted([9/14, 5/14])
    calculated_probs = sorted([prob_1, prob_2])

    if not np.allclose(calculated_probs, exact_probs):
        return f"Constraint check failed: The calculated probabilities {calculated_probs} do not match the exact probabilities {exact_probs}."

    # Check if the rounded probabilities match option A (0.64, 0.36)
    option_A_probs = sorted([0.64, 0.36])
    if not np.allclose(calculated_probs, option_A_probs, atol=0.005):
        return f"Constraint check failed: The calculated probabilities (rounded to two decimals) {[round(p, 2) for p in calculated_probs]} do not match option A's probabilities {option_A_probs}."

    # --- Step 4: Calculate the average value of A ---
    # Can be calculated as <A> = <alpha|A|alpha>
    average_value = np.vdot(alpha_normalized, A @ alpha_normalized).real

    # The exact average value is hbar/7
    exact_avg_val = hbar / 7
    if not np.isclose(average_value, exact_avg_val):
        return f"Constraint check failed: The calculated average value {average_value}*hbar does not match the exact value {exact_avg_val}*hbar."

    # --- Conclusion ---
    # All calculations match the values in option A.
    return "Correct"

# Run the check
result = check_correctness()
print(result)