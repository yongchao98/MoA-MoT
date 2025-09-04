import numpy as np

def check_quantum_measurement_answer():
    """
    This function verifies the solution to the given quantum mechanics problem.
    It calculates the probabilities and the expectation value from scratch and
    compares them to the values provided in the proposed answer.
    """
    # Let hbar = 1 for numerical calculations. We will add it back symbolically at the end.
    hbar = 1.0

    # --- Step 1: Define and normalize the initial state |alpha> ---
    # The state is proportional to (1+i)|up> + (2-i)|down>
    # In vector notation, |up> = [1, 0] and |down> = [0, 1]
    alpha_unnormalized = np.array([1 + 1j, 2 - 1j])

    # The squared norm is <alpha_unnormalized|alpha_unnormalized>
    # = (1-i)(1+i) + (2+i)(2-i) = (1^2 + 1^2) + (2^2 + 1^2) = 2 + 5 = 7
    norm_squared = np.vdot(alpha_unnormalized, alpha_unnormalized).real
    
    # The normalized state is |alpha> = (1/sqrt(norm)) * |alpha_unnormalized>
    alpha_normalized = alpha_unnormalized / np.sqrt(norm_squared)

    # --- Step 2: Define the operator A ---
    # The matrix A has elements Aij = hbar/2 if i != j, and 0 otherwise.
    A = (hbar / 2) * np.array([[0, 1], [1, 0]])

    # --- Step 3: Find the eigenstates and eigenvalues of A ---
    # These represent the possible states after measurement and the measured values.
    # We use np.linalg.eigh for Hermitian matrices, which A is.
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    # Eigenvalues are lambda_1 = -hbar/2, lambda_2 = +hbar/2
    # Eigenvectors are the columns of the 'eigenvectors' matrix.

    # --- Step 4: Calculate the measurement probabilities ---
    # The probability of measuring an eigenvalue lambda_i is P(i) = |<eigenstate_i|alpha>|^2
    probabilities = [np.abs(np.vdot(vec, alpha_normalized))**2 for vec in eigenvectors.T]
    
    # Sort probabilities in descending order to match the format of the answer options.
    calculated_probs = sorted(probabilities, reverse=True)
    
    # The exact probabilities are 9/14 and 5/14.
    exact_probs = [9/14, 5/14]

    # --- Step 5: Calculate the average (expectation) value of A ---
    # <A> = <alpha|A|alpha>
    expectation_value = np.vdot(alpha_normalized, A @ alpha_normalized).real
    
    # The exact expectation value is hbar/7. For hbar=1, this is 1/7.
    exact_exp_value_coeff = 1/7

    # --- Step 6: Verify the provided answer (Option B) ---
    # Option B states: Probabilities are 0.64, 0.36 and average value is hbar/7.
    answer_probs = [0.64, 0.36]
    answer_exp_value_coeff = 1/7

    # Check 1: Probabilities
    # The answer provides rounded values, so we check if our exact values round to the answer's values.
    # 9/14 = 0.6428... -> 0.64
    # 5/14 = 0.3571... -> 0.36
    # We use np.allclose with a tolerance that accounts for this rounding.
    if not np.allclose(calculated_probs, answer_probs, atol=0.005):
        return (f"Incorrect probabilities. "
                f"The answer provides probabilities {answer_probs}, but the calculated exact probabilities are "
                f"[{exact_probs[0]:.4f}, {exact_probs[1]:.4f}] which are approximately {calculated_probs}.")

    # Check 2: Average Value
    if not np.isclose(expectation_value, answer_exp_value_coeff):
        return (f"Incorrect average value. "
                f"The answer provides an average value of hbar/7, but the calculated value is hbar*{expectation_value:.4f}. "
                f"The exact coefficient for hbar should be {exact_exp_value_coeff:.4f}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_quantum_measurement_answer()
print(result)