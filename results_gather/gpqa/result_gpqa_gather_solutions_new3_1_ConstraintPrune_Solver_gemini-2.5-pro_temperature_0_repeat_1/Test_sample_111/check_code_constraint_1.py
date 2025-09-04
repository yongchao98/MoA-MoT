import numpy as np

def check_answer():
    """
    Checks the correctness of the provided answer for the quantum mechanics problem.
    
    The problem involves:
    1. Normalizing a quantum state |alpha>.
    2. Finding eigenvalues and eigenvectors of an operator A.
    3. Calculating the probability of measuring the system in each eigenstate.
    4. Calculating the average value of the operator A.
    """
    # Let hbar = 1 for calculation purposes. We can add it back symbolically at the end.
    hbar = 1.0

    # --- Step 1: Normalize the initial state |alpha> ---
    # The state is proportional to (1+i)|up> + (2-i)|down>
    # In vector form, where |up> = [1, 0] and |down> = [0, 1]:
    psi_unnormalized = np.array([1 + 1j, 2 - 1j])

    # Calculate the norm squared
    norm_sq = np.sum(np.abs(psi_unnormalized)**2)
    
    # The normalization constant is 1/sqrt(norm_sq)
    norm = np.sqrt(norm_sq)
    
    # The normalized state vector
    psi = psi_unnormalized / norm

    # --- Step 2: Define the operator A and find its eigenstates/eigenvalues ---
    # Aij = hbar/2 if i != j, and 0 otherwise.
    A = (hbar / 2) * np.array([[0, 1], [1, 0]])

    # Find eigenvalues and eigenvectors using numpy
    # The eigenvalues are the possible measurement outcomes.
    # The eigenvectors are the corresponding states.
    eigenvalues, eigenvectors = np.linalg.eig(A)

    # Note: np.linalg.eig returns eigenvectors as columns of the matrix.
    # Let's sort them for consistent ordering.
    sort_indices = np.argsort(eigenvalues)[::-1] # Sort from largest to smallest eigenvalue
    eigenvalues = eigenvalues[sort_indices]
    eigenvectors = eigenvectors[:, sort_indices]

    eigenvector_1 = eigenvectors[:, 0] # Corresponds to eigenvalue +hbar/2
    eigenvector_2 = eigenvectors[:, 1] # Corresponds to eigenvalue -hbar/2

    # --- Step 3: Calculate the probabilities ---
    # The probability of measuring an eigenvalue is the squared magnitude of the
    # inner product of the corresponding eigenvector with the state |alpha>.
    
    # Probability for the first eigenvalue (+hbar/2)
    # Inner product: <eigenvector_1 | psi>
    # In numpy, this is vdot(eigenvector_1, psi)
    inner_product_1 = np.vdot(eigenvector_1, psi)
    prob_1 = np.abs(inner_product_1)**2

    # Probability for the second eigenvalue (-hbar/2)
    inner_product_2 = np.vdot(eigenvector_2, psi)
    prob_2 = np.abs(inner_product_2)**2

    # --- Step 4: Calculate the average value of the operator ---
    # Method 1: Sum of (probability * eigenvalue)
    avg_value_1 = prob_1 * eigenvalues[0] + prob_2 * eigenvalues[1]
    
    # Method 2 (cross-check): <psi|A|psi>
    # A_psi = A @ psi
    # avg_value_2 = <psi|A_psi> = np.vdot(psi, A @ psi)
    # Both methods should yield the same result.

    # --- Step 5: Compare with the candidate answer A ---
    # Option A: Probabilities 0.64, 0.36 and average value hbar/7
    
    expected_prob_1 = 9.0 / 14.0  # ~0.6428
    expected_prob_2 = 5.0 / 14.0  # ~0.3571
    expected_avg_value = hbar / 7.0

    # Check probabilities
    # We check against the more precise fractional values
    if not np.isclose(prob_1, expected_prob_1):
        return f"Incorrect: The first probability is incorrect. Expected ~{expected_prob_1:.2f}, but calculated {prob_1:.4f}."
    if not np.isclose(prob_2, expected_prob_2):
        return f"Incorrect: The second probability is incorrect. Expected ~{expected_prob_2:.2f}, but calculated {prob_2:.4f}."
    
    # Check that the sum of probabilities is 1
    if not np.isclose(prob_1 + prob_2, 1.0):
        return f"Incorrect: The calculated probabilities do not sum to 1. Sum is {prob_1 + prob_2:.4f}."

    # Check average value
    if not np.isclose(avg_value_1, expected_avg_value):
        return f"Incorrect: The average value is incorrect. Expected hbar/7 ({expected_avg_value:.4f}*hbar), but calculated {avg_value_1:.4f}*hbar."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)