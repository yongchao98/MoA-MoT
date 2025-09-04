import numpy as np
import sympy

def check_quantum_measurement_answer():
    """
    Checks the correctness of the answer to the quantum mechanics problem.

    The problem involves:
    1. Normalizing a quantum state |alpha> proportional to (1+i)|up> + (2-i)|down>.
    2. Finding the eigenvalues and eigenvectors of an operator A where Aij = hbar/2 for i!=j and 0 otherwise.
    3. Calculating the probabilities of measuring the system in the eigenstates of A.
    4. Calculating the average value of A.
    """
    # Use sympy for symbolic representation of hbar
    hbar = sympy.Symbol('hbar')

    # --- Step 1: Normalize the initial state |alpha> ---
    # The unnormalized state vector in the {|up>, |down>} basis
    alpha_unnormalized = np.array([1 + 1j, 2 - 1j])

    # Calculate the squared norm
    norm_sq = np.vdot(alpha_unnormalized, alpha_unnormalized).real
    if not np.isclose(norm_sq, 7.0):
        return f"Incorrect normalization: The squared norm should be 7, but was calculated as {norm_sq}."

    # The normalization constant is 1/sqrt(norm_sq)
    norm = np.sqrt(norm_sq)
    alpha_normalized = alpha_unnormalized / norm

    # --- Step 2: Define the operator A and find its eigenstates/eigenvalues ---
    # We can perform calculations with hbar=1 and add it back symbolically at the end.
    A_matrix = 0.5 * np.array([[0, 1], [1, 0]])

    # Find eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(A_matrix)
    
    # Sort eigenvalues and corresponding eigenvectors for consistent ordering
    sort_indices = np.argsort(eigenvalues)[::-1] # Sort in descending order (+0.5, -0.5)
    eigenvalues = eigenvalues[sort_indices]
    eigenvectors = eigenvectors[:, sort_indices]

    # --- Step 3: Calculate the probabilities ---
    # Probability P_i = |<eigenvector_i | alpha_normalized>|^2
    # np.vdot computes the conjugate dot product, which is the correct inner product <v|a>
    prob1 = np.abs(np.vdot(eigenvectors[:, 0], alpha_normalized))**2
    prob2 = np.abs(np.vdot(eigenvectors[:, 1], alpha_normalized))**2

    # --- Step 4: Calculate the average value of A ---
    # Method 1: Sum of (probability * eigenvalue)
    avg_val_1 = prob1 * eigenvalues[0] + prob2 * eigenvalues[1]
    
    # Method 2: <alpha|A|alpha>
    avg_val_2 = np.vdot(alpha_normalized, A_matrix @ alpha_normalized).real

    # Check if both methods give the same result
    if not np.isclose(avg_val_1, avg_val_2):
        return f"Internal calculation error: Average value methods disagree. Method 1: {avg_val_1}, Method 2: {avg_val_2}."

    # The calculated average value is avg_val_1 (in units of hbar)
    calculated_avg_value_coeff = avg_val_1

    # --- Step 5: Compare with the provided answer (Option A) ---
    # Option A: Probabilities = 0.64, 0.36; Average value = hbar / 7
    
    # The exact probabilities are 9/14 and 5/14
    expected_prob1 = 9/14
    expected_prob2 = 5/14
    
    # The exact average value coefficient is 1/7
    expected_avg_value_coeff = 1/7

    # Check probabilities
    # The order might be swapped, so check both possibilities
    calc_probs = sorted([prob1, prob2], reverse=True)
    expected_probs = sorted([expected_prob1, expected_prob2], reverse=True)

    if not np.allclose(calc_probs, expected_probs, atol=1e-3):
        return (f"Incorrect probabilities. "
                f"Expected: {expected_probs[0]:.2f}, {expected_probs[1]:.2f}. "
                f"Calculated: {calc_probs[0]:.2f}, {calc_probs[1]:.2f}.")

    # Check average value
    if not np.isclose(calculated_avg_value_coeff, expected_avg_value_coeff, atol=1e-6):
        return (f"Incorrect average value. "
                f"Expected coefficient of hbar: {expected_avg_value_coeff}. "
                f"Calculated coefficient of hbar: {calculated_avg_value_coeff}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_quantum_measurement_answer()
print(result)