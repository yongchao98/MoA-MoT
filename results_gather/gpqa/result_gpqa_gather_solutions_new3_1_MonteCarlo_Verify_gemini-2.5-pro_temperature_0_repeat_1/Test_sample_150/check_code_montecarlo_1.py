import numpy as np

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    """
    # 1. Define the given state vector and operator matrix
    psi = np.array([-1, 2, 1], dtype=complex)
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=complex)
    
    target_eigenvalue = 0

    # 2. Find the eigenvalues and eigenvectors of the operator P
    eigenvalues, eigenvectors = np.linalg.eig(P)

    # 3. Find the eigenvector corresponding to the target eigenvalue (0)
    # np.isclose is used for safe floating-point comparison
    try:
        # Find the index of the eigenvalue that is close to 0
        idx = np.where(np.isclose(eigenvalues, target_eigenvalue))[0][0]
        # The corresponding eigenvector is the column at that index
        v_0 = eigenvectors[:, idx]
    except IndexError:
        return f"Constraint not satisfied: The target value {target_eigenvalue} is not an eigenvalue of the operator P. The calculated eigenvalues are {eigenvalues.real}."

    # 4. Calculate the probability using the general formula:
    # Prob(λ) = |<v_λ|ψ>|^2 / (<v_λ|v_λ> * <ψ|ψ>)
    # np.vdot computes the dot product (inner product) with complex conjugation on the first vector.
    
    # Numerator: |<v_0|ψ>|^2
    inner_product = np.vdot(v_0, psi)
    numerator = np.abs(inner_product)**2

    # Denominator: <v_0|v_0> * <ψ|ψ>
    # Note: np.linalg.eig returns normalized eigenvectors, so np.vdot(v_0, v_0) should be 1.0
    norm_v0_sq = np.vdot(v_0, v_0) 
    norm_psi_sq = np.vdot(psi, psi)
    denominator = norm_v0_sq * norm_psi_sq

    if np.isclose(denominator, 0):
        return "Error: The norm of the state vector is zero, which is physically invalid."

    calculated_probability = numerator / denominator

    # 5. Compare the calculated probability with the provided answer.
    # The final answer given is <<<B>>>, which corresponds to 1/3.
    expected_probability = 1/3

    if np.isclose(calculated_probability, expected_probability):
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability is {calculated_probability.real:.4f}, "
                f"which is approximately {numerator.real:.4f} / {denominator.real:.4f}. "
                f"The expected answer is 1/3, which is approximately {expected_probability:.4f}.")

# Run the check
result = check_answer()
print(result)