import numpy as np

def check_quantum_measurement_answer():
    """
    This function verifies the solution to the given quantum mechanics problem.
    It calculates the probabilities and the average value from first principles
    and compares them to the values given in the proposed answer (Option D).
    """
    # --- Problem Setup ---
    # Let hbar = 1 for numerical calculations. We will handle it symbolically in the final comparison.
    hbar = 1.0

    # The answer to check is D: Probabilities 0.64, 0.36 and average value hbar/7.
    expected_probs_list = sorted([0.64, 0.36])
    expected_avg_val = hbar / 7.0
    expected_avg_val_str = "hbar/7"

    # --- Step 1: Normalize the initial state ---
    # The state |α> is proportional to (1+i)|up> + (2-i)|down>.
    # In vector form, |ψ> = [1+i, 2-i].
    psi_unnormalized = np.array([1 + 1j, 2 - 1j], dtype=complex)

    # Calculate the squared norm: <ψ|ψ> = |1+i|^2 + |2-i|^2 = 2 + 5 = 7
    norm_sq = np.vdot(psi_unnormalized, psi_unnormalized).real
    if not np.isclose(norm_sq, 7.0):
        return f"Incorrect: The squared norm of the initial state should be 7, but was calculated as {norm_sq}."

    # The normalized state is |α> = (1/sqrt(7)) * |ψ>.
    alpha = psi_unnormalized / np.sqrt(norm_sq)

    # --- Step 2: Define the operator A and find its eigenstates/eigenvalues ---
    # Aij = hbar/2 if i!=j, and 0 otherwise.
    # A = (hbar/2) * [[0, 1], [1, 0]]
    A = (hbar / 2.0) * np.array([[0, 1], [1, 0]], dtype=float)

    # Find eigenvalues and eigenvectors using numpy.
    # The eigenvalues are the possible measurement outcomes.
    eigenvalues, eigenvectors = np.linalg.eig(A)

    # Sort eigenvalues and corresponding eigenvectors for consistent ordering.
    # This ensures we compare the correct probability with the correct eigenvalue.
    sort_indices = np.argsort(eigenvalues)[::-1]  # Descending order
    eigenvalues = eigenvalues[sort_indices]
    eigenvectors = eigenvectors[:, sort_indices]

    # Eigenvalues should be +hbar/2 and -hbar/2
    if not np.allclose(eigenvalues, [hbar/2.0, -hbar/2.0]):
        return f"Incorrect: The calculated eigenvalues {eigenvalues} do not match the expected [+hbar/2, -hbar/2]."

    # --- Step 3: Calculate the probabilities ---
    # The probability of measuring an eigenvalue λ_i is P(λ_i) = |<v_i|α>|^2,
    # where |v_i> is the corresponding normalized eigenvector.

    # Probability for eigenvalue λ_1 = +hbar/2
    v1 = eigenvectors[:, 0]
    inner_product_1 = np.vdot(v1, alpha)
    prob_1 = np.abs(inner_product_1)**2

    # Probability for eigenvalue λ_2 = -hbar/2
    v2 = eigenvectors[:, 1]
    inner_product_2 = np.vdot(v2, alpha)
    prob_2 = np.abs(inner_product_2)**2
    
    calculated_probs_list = sorted([prob_1, prob_2])

    # --- Step 4: Calculate the average value of A ---
    # Method 1: <A> = Σ P_i * λ_i
    calculated_avg_val = (prob_1 * eigenvalues[0] + prob_2 * eigenvalues[1]).real
    
    # Method 2 (cross-check): <A> = <α|A|α>
    # calculated_avg_val_check = np.vdot(alpha, A @ alpha).real
    # if not np.isclose(calculated_avg_val, calculated_avg_val_check):
    #     return "Internal check failed: Average value calculation methods disagree."

    # --- Step 5: Compare calculated results with the given answer ---
    # The given answer uses 2 decimal places, so we use a suitable tolerance.
    tolerance = 0.005

    # Check probabilities
    if not np.allclose(calculated_probs_list, expected_probs_list, atol=tolerance):
        return (f"Incorrect: The calculated probabilities [{calculated_probs_list[0]:.2f}, {calculated_probs_list[1]:.2f}] "
                f"do not match the given probabilities {expected_probs_list}.")

    # Check average value
    if not np.isclose(calculated_avg_val, expected_avg_val, atol=tolerance):
        # Format the calculated value as a fraction of hbar
        # calculated_avg_val / hbar = 1/7 = 0.1428...
        # We can express this as hbar * (calculated_avg_val / hbar)
        return (f"Incorrect: The calculated average value is approximately {calculated_avg_val/hbar:.4f}*hbar, "
                f"which does not match the given average value of {expected_avg_val_str}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
if __name__ == "__main__":
    # The exact fractional results are:
    # Probabilities: 9/14 (~0.643) and 5/14 (~0.357)
    # Average value: hbar/7
    # The provided answer D (0.64, 0.36, hbar/7) is a correctly rounded version of these results.
    # Our code should confirm this.
    result = check_quantum_measurement_answer()
    print(result)