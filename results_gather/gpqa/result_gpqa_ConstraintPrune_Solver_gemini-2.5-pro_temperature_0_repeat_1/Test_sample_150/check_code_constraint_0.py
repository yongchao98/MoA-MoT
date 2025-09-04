import numpy as np

def check_answer():
    """
    Checks the correctness of the given answer to the quantum mechanics problem.
    """
    # Define the state vector and the observable operator from the question
    psi = np.array([-1, 2, 1], dtype=float)
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=float)

    # The value to be measured
    measured_value = 0
    
    # The expected answer from option A
    expected_probability = 1/3

    # --- Step 1: Find eigenvalues and eigenvectors of the operator P ---
    try:
        eigenvalues, eigenvectors = np.linalg.eig(P)
    except np.linalg.LinAlgError:
        return "Error: The operator matrix P is singular or could not be processed."

    # --- Step 2: Check if the measured value is a possible outcome (an eigenvalue) ---
    # Find the index of the eigenvalue that is close to the measured_value
    eigenvalue_indices = np.where(np.isclose(eigenvalues, measured_value))[0]

    if len(eigenvalue_indices) == 0:
        return f"Incorrect. The value {measured_value} is not an eigenvalue of the operator P. The possible measurement outcomes are {np.round(eigenvalues, 5)}."

    # Get the normalized eigenvector corresponding to the measured value.
    # np.linalg.eig returns normalized eigenvectors.
    # If there are multiple eigenvectors for the same eigenvalue (degeneracy), we must project onto the subspace.
    # Here, the eigenvalue 0 is not degenerate.
    v_0 = eigenvectors[:, eigenvalue_indices[0]]

    # --- Step 3: Normalize the initial state vector psi ---
    norm_psi = np.linalg.norm(psi)
    if np.isclose(norm_psi, 0):
        return "Error: The initial state vector is a zero vector."
    psi_norm = psi / norm_psi

    # --- Step 4: Calculate the probability ---
    # Probability is the squared magnitude of the inner product (projection)
    # of the normalized state vector onto the normalized eigenvector.
    # For real vectors, conj() has no effect but is good practice for the general case.
    amplitude = np.dot(v_0.conj(), psi_norm)
    calculated_probability = np.abs(amplitude)**2

    # --- Step 5: Compare the calculated probability with the expected answer ---
    if np.isclose(calculated_probability, expected_probability):
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability is {calculated_probability:.4f}, "
                f"but the answer 'A' corresponds to a probability of {expected_probability:.4f}.")

# Run the check
result = check_answer()
print(result)