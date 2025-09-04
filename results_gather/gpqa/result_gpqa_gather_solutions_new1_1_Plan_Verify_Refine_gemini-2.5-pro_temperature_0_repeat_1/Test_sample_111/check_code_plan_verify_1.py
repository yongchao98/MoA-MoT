import numpy as np

def check_answer():
    """
    This function performs the quantum mechanics calculations and verifies the answer.
    """
    # Define hbar as 1 for numerical calculations. The final average value will be a multiple of hbar.
    hbar = 1.0

    # 1. Define the unnormalized state vector |alpha> ~ (1+i)|up> + (2-i)|down>
    # In the basis where |up> = [1, 0] and |down> = [0, 1]
    alpha_unnormalized = np.array([1 + 1j, 2 - 1j], dtype=complex)

    # 2. Normalize the state vector
    # The norm squared is <alpha|alpha> = (1-i)(1+i) + (2+i)(2-i) = 2 + 5 = 7
    norm = np.sqrt(np.vdot(alpha_unnormalized, alpha_unnormalized).real)
    if not np.isclose(norm, np.sqrt(7)):
        return f"Incorrect normalization constant. Expected sqrt(7), but got {norm}."
    alpha_normalized = alpha_unnormalized / norm

    # 3. Define the operator matrix A
    # Aij = hbar/2 if i != j, and 0 otherwise.
    A = (hbar / 2) * np.array([[0, 1], [1, 0]], dtype=float)

    # 4. Find the eigenvalues and eigenvectors of A
    # The eigenvalues are +/- hbar/2. The eigenvectors are (1/sqrt(2))*[1, 1] and (1/sqrt(2))*[1, -1]
    eigenvalues, eigenvectors = np.linalg.eigh(A)

    # 5. Calculate the probabilities of measuring each eigenvalue
    # P(lambda_i) = |<eigenvector_i|alpha>|^2
    # The order of eigenvalues/vectors from eigh() is ascending.
    # Eigenvalue -hbar/2 corresponds to eigenvector eigenvectors[:, 0]
    # Eigenvalue +hbar/2 corresponds to eigenvector eigenvectors[:, 1]
    prob_neg = np.abs(np.vdot(eigenvectors[:, 0], alpha_normalized))**2
    prob_pos = np.abs(np.vdot(eigenvectors[:, 1], alpha_normalized))**2
    
    # The exact probabilities are 5/14 and 9/14
    if not (np.isclose(prob_neg, 5/14) and np.isclose(prob_pos, 9/14)):
        return f"Calculated probabilities {prob_neg:.4f}, {prob_pos:.4f} do not match the expected exact values 5/14, 9/14."

    # 6. Calculate the average (expectation) value of A
    # <A> = <alpha|A|alpha>
    avg_value = np.vdot(alpha_normalized, A @ alpha_normalized).real
    
    # The exact average value is hbar/7
    if not np.isclose(avg_value, hbar / 7):
        return f"Calculated average value {avg_value:.4f}*hbar does not match the expected exact value hbar/7."

    # 7. Check against the provided answer (Option A)
    # Option A: Probabilities 0.64, 0.36 and average value hbar/7
    ans_probs = {0.64, 0.36}
    ans_avg_val_coeff = 1/7

    # Check if the set of calculated probabilities (rounded to 2 decimal places) matches the answer's probabilities
    calculated_probs_rounded = {round(prob_neg, 2), round(prob_pos, 2)}
    if calculated_probs_rounded != ans_probs:
        return f"Probability constraint not satisfied. Expected probabilities {ans_probs}, but calculated {calculated_probs_rounded}."

    # Check if the calculated average value coefficient matches the answer's coefficient
    if not np.isclose(avg_value / hbar, ans_avg_val_coeff):
        return f"Average value constraint not satisfied. Expected {ans_avg_val_coeff}*hbar, but calculated {avg_value/hbar}*hbar."

    return "Correct"

# Run the check
result = check_answer()
print(result)