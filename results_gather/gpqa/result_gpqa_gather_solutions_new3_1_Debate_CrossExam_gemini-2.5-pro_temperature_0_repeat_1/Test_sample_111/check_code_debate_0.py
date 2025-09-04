import numpy as np

def check_correctness():
    """
    This function verifies the solution to the quantum mechanics problem.
    It calculates the probabilities and the average value from first principles
    and compares them to the values provided in option A.
    """
    # For numerical calculations, we can set hbar = 1 and add it back symbolically at the end.
    hbar = 1.0

    # --- Step 1: Normalize the initial state ---
    # The state |alpha> is proportional to (1+i)|up> + (2-i)|down>.
    # In vector form, with |up> = [1, 0] and |down> = [0, 1]:
    psi_unnormalized = np.array([1 + 1j, 2 - 1j])

    # Calculate the squared norm: <psi|psi>
    norm_squared = np.vdot(psi_unnormalized, psi_unnormalized).real
    
    # The norm should be sqrt(7)
    if not np.isclose(norm_squared, 7.0):
        return f"Constraint check failed: The squared norm of the initial state should be 7, but was calculated as {norm_squared}."

    # Normalize the state vector
    psi_normalized = psi_unnormalized / np.sqrt(norm_squared)

    # --- Step 2: Define the operator A and find its properties ---
    # Aij = hbar/2 if i != j, and 0 otherwise.
    # A = (hbar/2) * [[0, 1], [1, 0]]
    A_matrix = (hbar / 2) * np.array([[0, 1], [1, 0]])

    # Find eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(A_matrix)

    # Sort the results for consistent comparison. We expect eigenvalues +hbar/2 and -hbar/2.
    # We sort in descending order so lambda1 = +hbar/2 and lambda2 = -hbar/2.
    sort_indices = np.argsort(eigenvalues)[::-1]
    lambda1 = eigenvalues[sort_indices[0]]
    lambda2 = eigenvalues[sort_indices[1]]
    v1 = eigenvectors[:, sort_indices[0]]
    v2 = eigenvectors[:, sort_indices[1]]

    # --- Step 3: Calculate the probabilities ---
    # The probability of measuring lambda_i is |<v_i|alpha>|^2
    # np.vdot(a, b) calculates the inner product a*.b
    prob1 = np.abs(np.vdot(v1, psi_normalized))**2
    prob2 = np.abs(np.vdot(v2, psi_normalized))**2

    # The probabilities must sum to 1
    if not np.isclose(prob1 + prob2, 1.0):
        return f"Constraint check failed: Probabilities do not sum to 1. Calculated sum: {prob1 + prob2}."

    # --- Step 4: Calculate the average value ---
    # <A> = sum(P_i * lambda_i)
    avg_A = prob1 * lambda1 + prob2 * lambda2
    
    # --- Step 5: Compare with the provided answer (Option A) ---
    # Option A: Probabilities 0.64, 0.36 and average value hbar/7
    expected_probs = {0.64, 0.36}
    expected_avg_A = hbar / 7.0

    # Check if the calculated probabilities match the expected ones (within a tolerance)
    # We check if the set of calculated probabilities matches the set of expected probabilities.
    # Exact values are 9/14 (~0.6428) and 5/14 (~0.3571)
    if not ( (np.isclose(prob1, 0.64, atol=0.005) and np.isclose(prob2, 0.36, atol=0.005)) or \
             (np.isclose(prob1, 0.36, atol=0.005) and np.isclose(prob2, 0.64, atol=0.005)) ):
        return f"Probabilities do not match. Expected {{0.64, 0.36}}, but calculated {{{prob1:.4f}, {prob2:.4f}}}."

    # Check if the calculated average value matches the expected one
    if not np.isclose(avg_A, expected_avg_A):
        return f"Average value does not match. Expected {expected_avg_A:.4f}*hbar, but calculated {avg_A:.4f}*hbar."

    return "Correct"

# Run the check
result = check_correctness()
print(result)