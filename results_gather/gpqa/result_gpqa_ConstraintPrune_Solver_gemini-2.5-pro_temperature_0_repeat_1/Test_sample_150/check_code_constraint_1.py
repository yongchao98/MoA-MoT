import numpy as np

def check_answer():
    """
    Checks the correctness of the given answer to the quantum mechanics problem.
    """
    # 1. Define the state vector and the observable matrix
    # State vector |ψ⟩
    psi = np.array([-1, 2, 1], dtype=complex)
    
    # Observable P
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=float)

    # The eigenvalue we are interested in
    target_eigenvalue = 0
    
    # The proposed answer from the LLM
    llm_answer = 1/3

    # 2. Normalize the state vector
    norm_psi = np.linalg.norm(psi)
    if np.isclose(norm_psi, 0):
        return "The state vector is a zero vector and cannot be normalized."
    psi_normalized = psi / norm_psi

    # 3. Find the eigenvalues and eigenvectors of the observable P
    eigenvalues, eigenvectors = np.linalg.eig(P)

    # 4. Find the eigenvector corresponding to the eigenvalue 0
    # Find the index of the eigenvalue that is close to 0
    try:
        # Use np.isclose for robust floating-point comparison
        eigenvalue_index = np.where(np.isclose(eigenvalues, target_eigenvalue))[0][0]
    except IndexError:
        return f"The observable P does not have an eigenvalue of {target_eigenvalue}. The calculated eigenvalues are {eigenvalues}."

    # The corresponding eigenvector (eigenvectors are columns of the returned matrix)
    # np.linalg.eig already returns normalized eigenvectors
    v_0 = eigenvectors[:, eigenvalue_index]

    # 5. Calculate the probability
    # Probability = |⟨v_0|ψ_normalized⟩|²
    # np.vdot computes the dot product with the conjugate of the first vector
    inner_product = np.vdot(v_0, psi_normalized)
    probability = np.abs(inner_product)**2

    # 6. Compare the calculated probability with the given answer
    if np.isclose(probability, llm_answer):
        return "Correct"
    else:
        return (f"Incorrect. The calculated probability is {probability:.4f}, "
                f"which is not equal to the proposed answer of {llm_answer:.4f} (1/3).")

# Execute the check
result = check_answer()
print(result)