import numpy as np

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    """
    # 1. Define the given state vector and observable operator
    psi = np.array([-1, 2, 1])
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ])

    # 2. Find the eigenvalues and eigenvectors of the operator P
    eigenvalues, eigenvectors = np.linalg.eig(P)
    
    # The eigenvectors are the columns of the returned matrix
    eigenvectors = eigenvectors.T

    # 3. Find the eigenvector corresponding to the eigenvalue 0
    target_eigenvalue = 0
    # Use a tolerance for floating point comparison
    try:
        # Find the index of the eigenvalue closest to 0
        idx = np.argmin(np.abs(eigenvalues - target_eigenvalue))
        v_0 = eigenvectors[idx]
    except (ValueError, IndexError):
        return "Incorrect: The value 0 is not an eigenvalue of the operator P."

    # Check if the found eigenvalue is indeed close to 0
    if not np.isclose(eigenvalues[idx], target_eigenvalue):
        return f"Incorrect: The value 0 is not an eigenvalue of the operator P. Closest eigenvalue is {eigenvalues[idx]}."

    # 4. Calculate the terms for the probability formula
    # The eigenvectors from np.linalg.eig are already normalized, so <v_0|v_0> = 1.
    
    # Calculate the squared norm of the state vector, <ψ|ψ>
    # np.vdot handles complex conjugation, which is good practice.
    psi_norm_sq = np.vdot(psi, psi).real
    
    # Calculate the inner product, <v_0|ψ>
    inner_product = np.vdot(v_0, psi)
    
    # Calculate the squared magnitude of the inner product, |<v_0|ψ>|^2
    inner_product_sq_mag = np.abs(inner_product)**2
    
    # 5. Calculate the final probability
    calculated_prob = inner_product_sq_mag / psi_norm_sq

    # 6. Compare with the expected answer
    expected_prob_val = 1/3
    
    if not np.isclose(calculated_prob, expected_prob_val):
        return f"Incorrect: The calculated probability is {calculated_prob:.4f}, but the expected probability is 1/3 ({expected_prob_val:.4f})."

    # 7. Check if the final answer choice matches the calculated value
    # The final answer block states the options are A) 1/3, B) 2/3, C) 1, D) sqrt(2/3)
    # and the final answer is <<<A>>>.
    options = {'A': 1/3, 'B': 2/3, 'C': 1, 'D': np.sqrt(2/3)}
    final_answer_choice = 'A'
    
    if final_answer_choice not in options:
        return f"Incorrect: The final answer choice '{final_answer_choice}' is not a valid option."

    if not np.isclose(options[final_answer_choice], calculated_prob):
        return f"Incorrect: The final answer choice is '{final_answer_choice}', which corresponds to a value of {options[final_answer_choice]:.4f}. However, the correctly calculated probability is {calculated_prob:.4f}."

    return "Correct"

# Run the check
result = check_answer()
if result == "Correct":
    print("Correct")
else:
    print(result)
