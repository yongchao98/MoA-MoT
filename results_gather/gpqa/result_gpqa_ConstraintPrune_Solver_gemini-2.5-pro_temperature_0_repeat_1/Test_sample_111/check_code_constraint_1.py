import numpy as np

def check_quantum_measurement_answer():
    """
    Checks the correctness of the given answer to the quantum mechanics problem.
    """
    # Let hbar = 1 for numerical calculations. The final answer for the average value
    # will be scaled by hbar.
    hbar = 1.0

    # --- Parse the given answer (C) ---
    # C) 0.64, 0.36 and hbar / 7
    try:
        ans_probs = {0.64, 0.36}
        ans_avg_val_numerical = hbar / 7.0
    except (ValueError, TypeError):
        return "Invalid format in the provided answer."

    # --- Step 1: Define and normalize the state vector |alpha> ---
    # The state is proportional to (1+i)|up> + (2-i)|down>
    # In vector form: [1+i, 2-i]
    alpha_unnormalized = np.array([1 + 1j, 2 - 1j])

    # The norm squared is <alpha|alpha> = (1-i)(1+i) + (2+i)(2-i) = 2 + 5 = 7
    # The normalization constant is 1/sqrt(7)
    norm = np.linalg.norm(alpha_unnormalized)
    if not np.isclose(norm, np.sqrt(7)):
        return f"Internal check failed: The norm of the unnormalized state vector is {norm}, but it should be sqrt(7)."
    
    alpha = alpha_unnormalized / norm

    # --- Step 2: Define the operator A ---
    # Aij = hbar/2 if i != j, and 0 otherwise. This is (hbar/2) * Pauli-X matrix.
    A = (hbar / 2) * np.array([[0, 1], 
                               [1, 0]])

    # --- Step 3: Find eigenvalues and eigenstates of A ---
    # We use np.linalg.eigh for Hermitian matrices. It returns sorted eigenvalues
    # and corresponding eigenvectors as columns.
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    # Eigenvalues will be [-hbar/2, +hbar/2]
    # Eigenvectors will be the columns of the 'eigenvectors' matrix.

    # --- Step 4: Calculate the measurement probabilities ---
    # Probability P_i = |<eigenstate_i | alpha>|^2
    # The inner product in numpy for complex vectors is np.vdot(v1, v2) = v1* . v2
    
    prob_list = []
    for i in range(len(eigenvalues)):
        eigenstate = eigenvectors[:, i]
        projection = np.vdot(eigenstate, alpha)
        probability = np.abs(projection)**2
        prob_list.append(probability)

    # Check if probabilities sum to 1 (a sanity check)
    if not np.isclose(sum(prob_list), 1.0):
        return f"Constraint not satisfied: Calculated probabilities {prob_list} do not sum to 1."

    # Compare the set of calculated probabilities with the answer's probabilities.
    # We sort them to perform a consistent comparison.
    sorted_calc_probs = sorted(prob_list)
    sorted_ans_probs = sorted(list(ans_probs))
    
    # The answer is rounded to two decimal places, so we use an appropriate tolerance.
    if not np.allclose(sorted_calc_probs, sorted_ans_probs, atol=0.005):
        return (f"Incorrect probabilities. "
                f"Expected (approximately): {sorted_ans_probs}, "
                f"Calculated: {[round(p, 4) for p in sorted_calc_probs]}. "
                f"The calculated exact probabilities are 5/14 and 9/14.")

    # --- Step 5: Calculate the average value of A ---
    # <A> = <alpha| A |alpha>
    # A @ alpha computes the matrix-vector product A|alpha>
    avg_val_calculated = np.vdot(alpha, A @ alpha)

    # Since A is Hermitian and |alpha> is a physical state, the expectation value must be real.
    if not np.isclose(avg_val_calculated.imag, 0):
        return f"Constraint not satisfied: The calculated average value {avg_val_calculated} is not a real number."
    
    avg_val_calculated = avg_val_calculated.real

    # Compare the calculated average value with the one from the answer.
    if not np.isclose(avg_val_calculated, ans_avg_val_numerical):
        return (f"Incorrect average value. "
                f"Expected: {ans_avg_val_numerical:.4f} (i.e., hbar/7), "
                f"Calculated: {avg_val_calculated:.4f}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_quantum_measurement_answer()
print(result)