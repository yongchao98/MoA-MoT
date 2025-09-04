import numpy as np
import math

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    """
    # Define the given state vector and operator
    psi = np.array([-1, 2, 1], dtype=float)
    P = np.array([
        [0, 1/math.sqrt(2), 0],
        [1/math.sqrt(2), 0, 1/math.sqrt(2)],
        [0, 1/math.sqrt(2), 0]
    ], dtype=float)

    # The measurement outcome of interest
    target_eigenvalue = 0

    # The options given in the question
    options = {
        'A': 1/3,
        'B': 1,
        'C': math.sqrt(2/3),
        'D': 2/3
    }
    
    # The final answer provided by the LLM
    llm_answer_option = 'A'
    expected_probability = options.get(llm_answer_option)

    if expected_probability is None:
        return f"The provided answer option '{llm_answer_option}' is not a valid option (A, B, C, D)."

    # --- Calculation ---

    # 1. Find the eigenvalues and eigenvectors of the operator P
    eigenvalues, eigenvectors = np.linalg.eig(P)

    # 2. Find the eigenvector corresponding to the target eigenvalue (0)
    # Use np.isclose for floating point comparison
    try:
        # Find the index of the eigenvalue that is close to 0
        idx = np.where(np.isclose(eigenvalues, target_eigenvalue))[0][0]
        # Get the corresponding eigenvector. Eigenvectors are columns in the returned matrix.
        v0 = eigenvectors[:, idx]
    except IndexError:
        return f"Constraint not satisfied: The value {target_eigenvalue} is not an eigenvalue of the operator P. The calculated eigenvalues are {eigenvalues}."

    # 3. Calculate the probability using the formula: Prob(λ) = |<v_λ|ψ>|² / <ψ|ψ>
    # Note: np.linalg.eig returns normalized eigenvectors, so <v_λ|v_λ> = 1.
    # We use the more general formula that works even with unnormalized vectors for clarity.
    
    # Let's find the unnormalized eigenvector for clarity, which is [c, 0, -c]. We can choose c=1.
    v0_unnormalized = np.array([1, 0, -1], dtype=float)

    # Calculate the inner product <v₀|ψ>
    inner_product = np.dot(v0_unnormalized, psi)
    
    # Calculate the squared magnitude of the inner product |<v₀|ψ>|²
    inner_product_sq = np.abs(inner_product)**2

    # Calculate the squared norm of the eigenvector <v₀|v₀>
    v0_norm_sq = np.dot(v0_unnormalized, v0_unnormalized)

    # Calculate the squared norm of the state vector <ψ|ψ>
    psi_norm_sq = np.dot(psi, psi)

    # Check for zero norms to avoid division by zero
    if np.isclose(v0_norm_sq, 0) or np.isclose(psi_norm_sq, 0):
        return "Constraint not satisfied: The state vector or the eigenvector cannot be a zero vector."

    # Calculate the final probability
    calculated_probability = inner_product_sq / (v0_norm_sq * psi_norm_sq)

    # --- Verification ---
    
    # Check if the calculated probability matches the expected answer
    if np.isclose(calculated_probability, expected_probability):
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability is {calculated_probability:.4f}, which is approximately {inner_product_sq}/{int(v0_norm_sq * psi_norm_sq)} = {inner_product_sq/v0_norm_sq}/{int(psi_norm_sq)} = {inner_product_sq/(v0_norm_sq*psi_norm_sq)}. "
                f"The answer claims the probability is {expected_probability:.4f} (Option {llm_answer_option}).")

# Run the check
result = check_answer()
print(result)