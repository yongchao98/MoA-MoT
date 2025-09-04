import numpy as np

def check_quantum_mechanics_problem():
    """
    This function checks the correctness of the provided solution to a quantum mechanics problem.

    The problem involves:
    1. A quantum state |alpha> proportional to (1+i)|up> + (2-i)|down>.
    2. An operator A with matrix elements Aij = hbar/2 for i!=j and 0 for i=j.
    3. Calculating the probabilities of measuring the system in the eigenstates of A.
    4. Calculating the average value of A.

    The function verifies the calculations against the provided answer D:
    Probabilities: 0.64, 0.36
    Average value: hbar/7
    """
    try:
        # For numerical calculations, we can set hbar = 1 and add it back symbolically at the end.
        hbar = 1.0

        # --- Step 1: Normalize the initial state |alpha> ---
        # The state is proportional to (1+i)|up> + (2-i)|down>.
        # In vector form with |up> = [1, 0] and |down> = [0, 1]:
        alpha_unnormalized = np.array([1 + 1j, 2 - 1j], dtype=complex)

        # The squared norm is <alpha_un|alpha_un> = |1+i|^2 + |2-i|^2 = (1^2+1^2) + (2^2+(-1)^2) = 2 + 5 = 7.
        norm_squared = np.vdot(alpha_unnormalized, alpha_unnormalized).real
        if not np.isclose(norm_squared, 7.0):
            return f"Constraint check failed: The squared norm of the initial state should be 7, but was calculated as {norm_squared}."

        # The normalized state vector |alpha>
        alpha = alpha_unnormalized / np.sqrt(norm_squared)

        # --- Step 2: Define the operator A ---
        # Aij = hbar/2 if i!=j, and 0 otherwise.
        A = (hbar / 2) * np.array([[0, 1], [1, 0]], dtype=float)

        # --- Step 3: Find eigenvalues and eigenstates of A ---
        # The characteristic equation det(A - lambda*I) = 0 gives lambda^2 - (hbar/2)^2 = 0.
        # So, the eigenvalues are lambda = +/- hbar/2.
        eigenvalues, eigenvectors = np.linalg.eig(A)

        # Sort for consistent comparison
        eigenvalues.sort()
        expected_eigenvalues = np.array([-hbar/2, hbar/2])
        if not np.allclose(eigenvalues, expected_eigenvalues):
            return f"Constraint check failed: The calculated eigenvalues {eigenvalues} do not match the expected eigenvalues {expected_eigenvalues} (with hbar=1)."

        # Get the eigenvectors corresponding to +hbar/2 and -hbar/2
        # np.linalg.eig returns eigenvectors as columns.
        # Let's re-run eig to get the original correspondence before sorting.
        eigenvalues, eigenvectors = np.linalg.eig(A)
        idx_plus = np.where(np.isclose(eigenvalues, hbar/2))[0][0]
        idx_minus = np.where(np.isclose(eigenvalues, -hbar/2))[0][0]
        
        eigenstate_plus = eigenvectors[:, idx_plus]
        eigenstate_minus = eigenvectors[:, idx_minus]

        # --- Step 4: Calculate the probabilities ---
        # The probability of measuring an eigenvalue lambda_i is P(lambda_i) = |<eigenstate_i|alpha>|^2.
        # np.vdot(a, b) computes the inner product <a|b> = a*.b
        prob_plus = np.abs(np.vdot(eigenstate_plus, alpha))**2
        prob_minus = np.abs(np.vdot(eigenstate_minus, alpha))**2

        # The exact probabilities from the solution are 9/14 and 5/14.
        expected_prob_1 = 9.0 / 14.0  # approx 0.6428
        expected_prob_2 = 5.0 / 14.0  # approx 0.3571

        # Check if the calculated probabilities match the exact values, regardless of order.
        if not ((np.isclose(prob_plus, expected_prob_1) and np.isclose(prob_minus, expected_prob_2)) or
                (np.isclose(prob_plus, expected_prob_2) and np.isclose(prob_minus, expected_prob_1))):
            return f"Probability calculation is incorrect. Expected probabilities [9/14, 5/14], but got [{prob_plus:.4f}, {prob_minus:.4f}]."

        # Check if probabilities sum to 1
        if not np.isclose(prob_plus + prob_minus, 1.0):
            return f"Constraint check failed: The sum of probabilities is {prob_plus + prob_minus}, which is not 1."

        # --- Step 5: Calculate the average value of A ---
        # The average value <A> can be calculated as <alpha|A|alpha>.
        avg_A = np.vdot(alpha, A @ alpha).real

        # The expected average value from the solution is hbar/7.
        expected_avg_A = hbar / 7.0
        if not np.isclose(avg_A, expected_avg_A):
            return f"Average value calculation is incorrect. Expected hbar/7 ({expected_avg_A:.4f} with hbar=1), but got {avg_A:.4f}."

        # --- Final check against option D ---
        # Option D provides probabilities 0.64, 0.36 and average value hbar/7.
        # Our calculated probabilities (0.6428..., 0.3571...) round to the values in option D.
        # Our calculated average value matches hbar/7.
        # All calculations align with the provided solution and option D.

        return "Correct"

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_quantum_mechanics_problem()
print(result)