import numpy as np

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    """
    # For calculation, we can set hbar = 1 and add it back symbolically at the end.
    hbar = 1.0
    
    # --- Step 1: Normalize the initial state |alpha> ---
    # Unnormalized state vector: (1+i)|up> + (2-i)|down>
    psi_un = np.array([1 + 1j, 2 - 1j])
    
    # Calculate the norm
    norm_sq = np.vdot(psi_un, psi_un).real
    if not np.isclose(norm_sq, 7.0):
        return f"Incorrect normalization: <psi|psi> should be 7, but calculated as {norm_sq}."
        
    norm = np.sqrt(norm_sq)
    
    # Normalized state vector
    alpha = psi_un / norm
    
    # --- Step 2: Define the operator A and find its eigenstates/eigenvalues ---
    # A = (hbar/2) * [[0, 1], [1, 0]]
    A = (hbar / 2) * np.array([[0, 1], [1, 0]])
    
    # Find eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(A)
    
    # Sort eigenvalues and corresponding eigenvectors for consistent comparison
    sort_indices = np.argsort(eigenvalues)[::-1] # Sort in descending order
    eigenvalues = eigenvalues[sort_indices]
    eigenvectors = eigenvectors[:, sort_indices]
    
    eigenstate_1 = eigenvectors[:, 0]
    eigenstate_2 = eigenvectors[:, 1]
    
    # --- Step 3: Calculate the probabilities ---
    # Probability P_i = |<eigenstate_i | alpha>|^2
    # <v|a> is the inner product, which is v*.dot(a) in numpy
    prob_1 = np.abs(np.vdot(eigenstate_1, alpha))**2
    prob_2 = np.abs(np.vdot(eigenstate_2, alpha))**2
    
    # --- Step 4: Calculate the average value <A> ---
    # Method 1: Sum of P_i * lambda_i
    avg_val_1 = prob_1 * eigenvalues[0] + prob_2 * eigenvalues[1]
    
    # Method 2: <alpha|A|alpha>
    avg_val_2 = np.vdot(alpha, A @ alpha).real
    
    if not np.isclose(avg_val_1, avg_val_2):
        return f"Inconsistency in average value calculation: Method 1 gives {avg_val_1}, Method 2 gives {avg_val_2}."

    # --- Step 5: Compare with the chosen answer (D) ---
    # Option D provides: Probs = 0.64, 0.36 and Avg_val = hbar/7
    
    # The calculated probabilities are 9/14 and 5/14
    expected_probs = np.array([9/14, 5/14])
    calculated_probs = np.array([prob_1, prob_2])
    
    # Check probabilities
    if not np.allclose(np.sort(calculated_probs), np.sort(expected_probs)):
        return f"Probability check failed. Expected [~0.357, ~0.643], but got {np.sort(calculated_probs)}."

    # Check average value
    expected_avg_val = hbar / 7.0
    if not np.isclose(avg_val_1, expected_avg_val):
        return f"Average value check failed. Expected {expected_avg_val}*hbar, but got {avg_val_1}*hbar."

    # Check against the rounded values in option D
    option_d_probs = np.array([0.36, 0.64])
    if not np.allclose(np.sort(calculated_probs), np.sort(option_d_probs), atol=0.005):
         return f"Probability check against option D failed. Option D has {option_d_probs}, but calculated values are {np.sort(calculated_probs)}."

    return "Correct"

# Run the check
result = check_answer()
print(result)