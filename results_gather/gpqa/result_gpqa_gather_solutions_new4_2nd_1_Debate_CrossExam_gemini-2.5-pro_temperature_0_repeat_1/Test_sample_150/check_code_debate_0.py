import numpy as np

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the quantum mechanics problem.
    It recalculates the probability from first principles and compares it to the provided answer.
    """
    # 1. Define the problem's parameters from the question
    # State vector |ψ⟩
    psi = np.array([-1, 2, 1], dtype=float)

    # Observable operator P
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=float)

    # The eigenvalue of interest is 0
    target_eigenvalue = 0

    # The options as given in the question text
    options = {
        'A': np.sqrt(2)/3,
        'B': 2/3,
        'C': 1.0,
        'D': 1/3
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_value = 1/3
    llm_answer_option = 'D'

    # 2. Find the eigenvector for the target eigenvalue
    try:
        eigenvalues, eigenvectors = np.linalg.eig(P)
        # Find the index of the eigenvalue that is close to 0
        idx = np.where(np.isclose(eigenvalues, target_eigenvalue))[0][0]
        # Get the corresponding eigenvector |v₀⟩. Eigenvectors are columns.
        v0 = eigenvectors[:, idx]
    except (np.linalg.LinAlgError, IndexError) as e:
        return f"Error during eigenvalue/eigenvector calculation: {e}"

    # 3. Calculate the probability using the general formula:
    # Prob(λ) = |⟨v_λ|ψ⟩|² / (⟨v_λ|v_λ⟩ * ⟨ψ|ψ⟩)
    # Note: np.linalg.eig returns normalized eigenvectors, so ⟨v₀|v₀⟩ is 1.
    # We calculate all terms for clarity and robustness.

    # Numerator: |⟨v₀|ψ⟩|²
    inner_product = np.dot(v0.conj(), psi)
    numerator = np.abs(inner_product)**2

    # Denominator: ⟨v₀|v₀⟩ * ⟨ψ|ψ⟩
    norm_v0_sq = np.dot(v0.conj(), v0)
    norm_psi_sq = np.dot(psi.conj(), psi)
    
    if np.isclose(norm_psi_sq, 0):
        return "Error: The norm of the state vector is zero."

    denominator = norm_v0_sq * norm_psi_sq
    calculated_probability = numerator / denominator

    # 4. Check the correctness of the LLM's answer
    
    # Check 1: Is the numerical result correct?
    if not np.isclose(calculated_probability, llm_answer_value):
        return (f"Incorrect numerical result. The calculated probability is {calculated_probability:.5f}, "
                f"but the LLM's answer is {llm_answer_value:.5f}.")

    # Check 2: Does the selected option letter correspond to the correct value?
    if llm_answer_option not in options:
        return f"Invalid option '{llm_answer_option}'. Valid options are A, B, C, D."
        
    if not np.isclose(options[llm_answer_option], calculated_probability):
         return (f"Incorrect final answer. The LLM chose option '{llm_answer_option}', which corresponds to a value of "
                f"{options[llm_answer_option]:.5f}. However, the calculated correct probability is {calculated_probability:.5f}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)