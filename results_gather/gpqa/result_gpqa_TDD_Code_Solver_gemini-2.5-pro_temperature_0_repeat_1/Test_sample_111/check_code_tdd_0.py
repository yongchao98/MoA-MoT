import numpy as np

def check_answer():
    """
    This function calculates the solution to the quantum mechanics problem
    and checks if it matches the provided answer C.
    """
    # Let hbar = 1 for numerical calculations. We can add it back symbolically later.
    hbar = 1.0

    # --- Step 1: Define and normalize the initial state |alpha> ---
    # |alpha> is proportional to (1+i)|up> + (2-i)|down>
    # In vector form, this is [1+i, 2-i]
    alpha_unnormalized = np.array([1 + 1j, 2 - 1j])

    # The squared norm is the sum of the squared magnitudes of the components.
    # |1+i|^2 = 1^2 + 1^2 = 2
    # |2-i|^2 = 2^2 + (-1)^2 = 5
    # Total squared norm = 2 + 5 = 7
    norm = np.sqrt(np.vdot(alpha_unnormalized, alpha_unnormalized).real)
    if not np.isclose(norm**2, 7.0):
        return f"Incorrect: The squared norm of the state vector should be 7, but was calculated as {norm**2}."
    
    alpha_normalized = alpha_unnormalized / norm

    # --- Step 2: Define the operator A and find its eigenstates ---
    # Aij = hbar/2 if i != j, and 0 otherwise.
    # A = (hbar/2) * [[0, 1], [1, 0]]
    A = (hbar / 2) * np.array([[0, 1], [1, 0]])

    # Find eigenvalues and eigenvectors. eigh is used for Hermitian matrices.
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    # Eigenvalues will be [-0.5, 0.5] (or [-hbar/2, +hbar/2])
    # Eigenvectors will be the columns of the 'eigenvectors' matrix.
    # e.g., for +hbar/2: (1/sqrt(2))*[1, 1]
    # e.g., for -hbar/2: (1/sqrt(2))*[1, -1]

    # --- Step 3: Calculate the probabilities ---
    # The probability of measuring an eigenstate |e_i> is P_i = |<e_i|alpha>|^2
    calculated_probs = []
    for i in range(len(eigenvalues)):
        eigenvector_i = eigenvectors[:, i]
        # <e_i|alpha> is the inner product
        projection = np.vdot(eigenvector_i, alpha_normalized)
        prob = np.abs(projection)**2
        calculated_probs.append(prob)

    # The exact probabilities are 9/14 and 5/14.
    # 9/14 ~= 0.642857
    # 5/14 ~= 0.357143
    exact_probs = [9/14, 5/14]

    # Sort both lists in descending order for consistent comparison.
    calculated_probs.sort(reverse=True)
    exact_probs.sort(reverse=True)

    if not np.allclose(calculated_probs, exact_probs):
        return f"Incorrect: Calculated probabilities {calculated_probs} do not match the exact probabilities {exact_probs}."

    # --- Step 4: Calculate the average value <A> ---
    # <A> = <alpha|A|alpha>
    A_alpha = A @ alpha_normalized
    calculated_avg_value = np.vdot(alpha_normalized, A_alpha).real
    
    # The exact average value is hbar/7. With hbar=1, this is 1/7.
    exact_avg_value = hbar / 7.0
    
    if not np.isclose(calculated_avg_value, exact_avg_value):
        return f"Incorrect: Calculated average value {calculated_avg_value}*hbar does not match the exact average value {exact_avg_value}*hbar."

    # --- Step 5: Compare with the given answer C ---
    # Answer C: Probabilities 0.64, 0.36 and average value hbar/7
    answer_c_probs = [0.64, 0.36]
    answer_c_avg_value_coeff = 1/7.0

    # Check if the rounded exact probabilities match the answer
    rounded_exact_probs = np.round(exact_probs, 2)
    answer_c_probs.sort(reverse=True) # Ensure order is the same
    if not np.allclose(rounded_exact_probs, answer_c_probs):
        return f"Incorrect: The probabilities from answer C {answer_c_probs} are not the correctly rounded values of the exact probabilities {exact_probs} (which round to {rounded_exact_probs})."

    # Check if the average value from the answer matches the exact one
    if not np.isclose(answer_c_avg_value_coeff, exact_avg_value):
        return f"Incorrect: The average value from answer C ({answer_c_avg_value_coeff}*hbar) does not match the exact average value ({exact_avg_value}*hbar)."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)