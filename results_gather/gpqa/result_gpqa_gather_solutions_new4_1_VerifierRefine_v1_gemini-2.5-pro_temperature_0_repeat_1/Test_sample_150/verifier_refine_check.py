import numpy as np

def check_quantum_probability():
    """
    Checks the calculation for the probability of measuring an eigenvalue in a quantum system.

    The problem is:
    State vector |ψ⟩ = [-1, 2, 1]
    Observable P = [[0, 1/√2, 0], [1/√2, 0, 1/√2], [0, 1/√2, 0]]
    Calculate the probability of measuring the eigenvalue 0.

    The expected answer is C, which corresponds to 1/3.
    """
    # Define the state vector and the observable operator
    psi = np.array([-1, 2, 1], dtype=complex)
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=complex)

    # The measurement outcome we are interested in
    target_eigenvalue = 0

    # The expected probability from the answer C
    expected_probability = 1/3

    # Step 1: Find the eigenvalues and eigenvectors of the operator P
    try:
        eigenvalues, eigenvectors = np.linalg.eig(P)
    except np.linalg.LinAlgError as e:
        return f"Error during eigenvalue decomposition: {e}"

    # Step 2: Find the eigenvector corresponding to the target eigenvalue (0)
    # We find the index of the eigenvalue closest to our target
    try:
        zero_eig_index = np.argmin(np.abs(eigenvalues - target_eigenvalue))
    except ValueError:
        return "Could not find the target eigenvalue 0 among the calculated eigenvalues."

    # Check if the found eigenvalue is indeed close to 0
    if not np.isclose(eigenvalues[zero_eig_index], target_eigenvalue):
        return f"The value 0 is not an eigenvalue of the operator P. The calculated eigenvalues are {eigenvalues}."

    # The corresponding eigenvector (eigenvectors are the columns of the matrix)
    # The eigenvectors from np.linalg.eig are already normalized.
    v0 = eigenvectors[:, zero_eig_index]

    # Step 3: Calculate the probability using the formula: Prob(λ) = |⟨v_λ|ψ⟩|² / ⟨ψ|ψ⟩
    # Since v0 is normalized, ⟨v_0|v_0⟩ = 1. The formula simplifies to |⟨v_0|ψ⟩|² / ⟨ψ|ψ⟩

    # Calculate the squared norm of the state vector, ⟨ψ|ψ⟩
    psi_norm_sq = np.dot(psi.conj(), psi).real

    # Calculate the inner product ⟨v_0|ψ⟩
    inner_product = np.dot(v0.conj(), psi)

    # Calculate the squared magnitude of the inner product, |⟨v_0|ψ⟩|²
    inner_product_sq_mag = np.abs(inner_product)**2

    # Calculate the final probability
    calculated_probability = inner_product_sq_mag / psi_norm_sq

    # Step 4: Compare the calculated probability with the expected answer
    if np.isclose(calculated_probability, expected_probability):
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability is {calculated_probability:.6f}, "
                f"which is not equal to the expected probability of 1/3 ({expected_probability:.6f}). "
                f"The final answer 'C' corresponds to 1/3, but the calculation does not match.")

# Run the check
result = check_quantum_probability()
print(result)