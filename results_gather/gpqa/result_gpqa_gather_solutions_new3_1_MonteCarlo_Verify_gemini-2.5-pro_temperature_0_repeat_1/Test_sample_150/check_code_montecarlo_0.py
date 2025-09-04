import numpy as np

def check_correctness():
    """
    Checks the correctness of the calculated probability for the given quantum mechanics problem.

    The function performs the following steps:
    1. Defines the state vector `psi` and the observable operator `P`.
    2. Calculates the eigenvalues and eigenvectors of `P`.
    3. Verifies that the target value (0) is indeed an eigenvalue.
    4. Selects the eigenvector corresponding to the eigenvalue 0.
    5. Calculates the probability using Born's rule: Prob(λ) = |<v_λ|ψ>|^2 / <ψ|ψ>.
    6. Compares the calculated probability with the expected answer (1/3).
    """
    # Define the state vector and the observable operator from the question
    psi = np.array([-1, 2, 1], dtype=complex)
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ])
    
    # The target outcome is the eigenvalue 0
    target_eigenvalue = 0
    
    # The expected probability from the answer is 1/3
    expected_probability = 1/3

    # Step 1: Find eigenvalues and eigenvectors of the operator P
    try:
        eigenvalues, eigenvectors = np.linalg.eig(P)
    except np.linalg.LinAlgError:
        return "Error: Could not compute eigenvalues and eigenvectors for the operator P."

    # Step 2: Check if the target value is an eigenvalue and get its eigenvector
    # Use a tolerance for floating-point comparison
    matching_indices = np.where(np.isclose(eigenvalues, target_eigenvalue))
    
    if len(matching_indices[0]) == 0:
        return (f"Constraint not satisfied: The target value {target_eigenvalue} is not an eigenvalue of the operator P. "
                f"The calculated eigenvalues are {eigenvalues}.")
    
    # The eigenvector corresponding to the target eigenvalue.
    # np.linalg.eig returns normalized eigenvectors.
    v_lambda = eigenvectors[:, matching_indices[0][0]]

    # Step 3: Calculate the squared norm of the state vector (denominator)
    psi_norm_sq = np.dot(psi.conj(), psi).real
    if np.isclose(psi_norm_sq, 0):
        return "Constraint not satisfied: The state vector cannot be a zero vector."

    # Step 4: Calculate the squared magnitude of the inner product (numerator)
    inner_product = np.dot(v_lambda.conj(), psi)
    inner_product_sq_mag = np.abs(inner_product)**2

    # Step 5: Calculate the final probability
    calculated_probability = inner_product_sq_mag / psi_norm_sq

    # Step 6: Compare the calculated probability with the expected answer
    if np.isclose(calculated_probability, expected_probability):
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability is {calculated_probability:.4f}, "
                f"but the expected answer is 1/3 ({expected_probability:.4f}).")

# Run the check
result = check_correctness()
print(result)