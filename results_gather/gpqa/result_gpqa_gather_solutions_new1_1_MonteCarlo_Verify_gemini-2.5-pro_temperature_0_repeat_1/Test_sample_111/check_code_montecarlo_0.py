import numpy as np

def check_correctness():
    """
    This function checks the correctness of the answer to the quantum mechanics problem.

    The problem asks for:
    1. The probability of measuring the particle in each of the eigenstates of operator A.
    2. The average value of operator A.

    The steps are:
    1. Define the initial state and normalize it.
    2. Define the operator A.
    3. Find the eigenvalues and eigenvectors of A.
    4. Calculate the probabilities by projecting the state onto the eigenvectors.
    5. Calculate the average value of A for the given state.
    6. Compare the calculated results with the values in the proposed correct option C.
    """
    # For numerical calculations, we can set hbar = 1 and check the coefficient of the average value.
    hbar = 1.0

    # 1. Define and normalize the initial state |alpha>
    # |alpha> is proportional to (1+i)|up> + (2-i)|down>
    # In vector form, this is proportional to [1+1j, 2-1j]
    psi_unnormalized = np.array([1 + 1j, 2 - 1j])

    # Calculate the norm squared: <psi|psi> = |c1|^2 + |c2|^2
    norm_squared = np.vdot(psi_unnormalized, psi_unnormalized).real
    # The normalization constant is 1/sqrt(norm_squared)
    norm = np.sqrt(norm_squared)
    psi_normalized = psi_unnormalized / norm

    # Check if the normalization factor is correct (sqrt(7))
    if not np.isclose(norm, np.sqrt(7)):
        return f"Incorrect normalization constant. Expected sqrt(7), but calculation yields {norm:.4f}."

    # 2. Define the operator A
    # Aij = hbar/2 if i != j, and 0 otherwise.
    A = (hbar / 2) * np.array([[0, 1], [1, 0]])

    # 3. Find eigenvalues and eigenvectors of A
    eigenvalues, eigenvectors = np.linalg.eig(A)
    # The eigenvectors are the columns of the returned matrix.
    # We sort them by eigenvalue to ensure a consistent order for checking.
    sort_indices = np.argsort(eigenvalues)[::-1]  # Sort in descending order
    eigenvalues = eigenvalues[sort_indices]
    eigenvectors = eigenvectors[:, sort_indices]

    # The expected eigenvalues are +hbar/2 and -hbar/2
    expected_eigenvalues = np.array([hbar / 2, -hbar / 2])
    if not np.allclose(eigenvalues, expected_eigenvalues):
        return f"Calculated eigenvalues {eigenvalues} do not match expected eigenvalues {expected_eigenvalues}."

    # 4. Calculate the probabilities
    # Probability P_i = |<eigenvector_i | psi_normalized>|^2
    # The inner product <v|psi> is calculated by np.vdot(v, psi)
    prob1 = np.abs(np.vdot(eigenvectors[:, 0], psi_normalized))**2
    prob2 = np.abs(np.vdot(eigenvectors[:, 1], psi_normalized))**2

    # 5. Calculate the average value of A
    # Method 1: Sum of (probability * eigenvalue)
    avg_value_p = prob1 * eigenvalues[0] + prob2 * eigenvalues[1]
    # Method 2: <psi|A|psi> (for cross-checking)
    avg_value_direct = np.vdot(psi_normalized, A @ psi_normalized).real
    
    if not np.isclose(avg_value_p, avg_value_direct):
        return "Internal check failed: Average value calculations do not match."

    # 6. Compare with the proposed answer (Option C)
    # Option C: Probabilities 0.64, 0.36 and average value hbar/7.
    
    # The exact calculated probabilities are 9/14 and 5/14.
    expected_probs = {9.0 / 14.0, 5.0 / 14.0}
    calculated_probs = {prob1, prob2}

    # Check if the calculated probabilities match the exact values.
    # We check if the sets are close element-wise.
    is_prob_match = all(any(np.isclose(c, e) for e in expected_probs) for c in calculated_probs)
    if not is_prob_match:
        return f"Calculated probabilities {calculated_probs} do not match expected exact probabilities {expected_probs} (9/14 and 5/14)."

    # The exact average value is hbar/7.
    expected_avg_value = hbar / 7.0
    if not np.isclose(avg_value_direct, expected_avg_value):
        return f"Calculated average value {avg_value_direct:.4f} does not match expected value {expected_avg_value:.4f} (hbar/7)."

    # The values in option C are rounded. Let's check if our results match the rounded values.
    # Probs: 0.64, 0.36. Avg value coeff: 1/7.
    option_c_probs = {0.64, 0.36}
    is_rounded_prob_match = all(any(np.isclose(c, o, atol=0.005) for o in option_c_probs) for c in calculated_probs)
    if not is_rounded_prob_match:
        return f"Calculated probabilities {calculated_probs} do not match the rounded values in option C {option_c_probs}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)