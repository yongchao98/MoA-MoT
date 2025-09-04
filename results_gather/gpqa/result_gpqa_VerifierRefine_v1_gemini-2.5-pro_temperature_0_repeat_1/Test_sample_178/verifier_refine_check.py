import numpy as np
from scipy.linalg import expm

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer.
    It defines the matrices, checks their properties, evaluates each statement,
    and compares the findings with the given answer and its reasoning.
    """
    # Define a tolerance for floating point comparisons
    TOL = 1e-9

    # Helper functions to check matrix properties
    def is_hermitian(A):
        """Checks if a matrix is Hermitian (A = A†)."""
        return np.allclose(A, A.conj().T, atol=TOL)

    def is_anti_hermitian(A):
        """Checks if a matrix is anti-Hermitian (A† = -A)."""
        return np.allclose(A.conj().T, -A, atol=TOL)

    def is_unitary(A):
        """Checks if a matrix is unitary (A†A = AA† = I)."""
        n = A.shape[0]
        identity = np.identity(n)
        return np.allclose(A @ A.conj().T, identity, atol=TOL) and \
               np.allclose(A.conj().T @ A, identity, atol=TOL)

    def is_quantum_state(rho):
        """
        Checks if a matrix represents a quantum state (density matrix).
        Returns a tuple (bool, reason_string).
        """
        # 1. Must be Hermitian
        if not is_hermitian(rho):
            return False, "it is not Hermitian."
        
        # 2. Trace must be 1
        if not np.isclose(np.trace(rho).real, 1.0, atol=TOL):
            return False, f"its trace is {np.trace(rho).real}, not 1."
            
        # 3. Must be positive semi-definite (all eigenvalues >= 0)
        # Since it's Hermitian, eigenvalues are real.
        try:
            eigenvalues = np.linalg.eigvalsh(rho)
            if not np.all(eigenvalues >= -TOL):
                return False, f"it is not positive semi-definite (has negative eigenvalues: {eigenvalues[eigenvalues < 0]})."
        except np.linalg.LinAlgError:
            return False, "eigenvalue calculation failed."
            
        return True, "it is a valid quantum state."

    # --- Main analysis ---

    # Define the matrices from the problem
    try:
        W = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=float)
        X = np.array([[1j, -1, 2j], [1, 0, 1], [2j, -1, -1j]], dtype=complex)
        Y = np.array([[0.5, 0.1, 0.2], [0.1, 0.25, 0.1], [0.2, 0.1, 0.25]], dtype=float)
        Z = np.array([[3, 2j, 5], [-2j, -2, -4j], [5, 4j, 4]], dtype=complex)
    except Exception as e:
        return f"Failed to define matrices: {e}"

    # The provided answer from the LLM is D.
    llm_answer = 'D'

    # --- Verify the properties claimed in the LLM's reasoning ---
    reasoning_checks = {
        "W is Unitary": is_unitary(W),
        "X is not Unitary": not is_unitary(X),
        "X is anti-Hermitian": is_anti_hermitian(X),
        "e^X is Unitary": is_unitary(expm(X)),
        "Z is Hermitian": is_hermitian(Z),
        "X is not Hermitian": not is_hermitian(X),
        "Y is a quantum state": is_quantum_state(Y)[0],
        "(e^X)*Y*(e^{-X}) is a quantum state": is_quantum_state(expm(X) @ Y @ expm(-X))[0]
    }

    for claim, result in reasoning_checks.items():
        if not result:
            return f"Incorrect. The reasoning provided in the answer is flawed. The claim '{claim}' is false."

    # --- Evaluate each statement to find the correct one ---
    # A) W and X represent the evolution operator (unitary).
    statement_A_correct = is_unitary(W) and is_unitary(X)

    # B) e^X changes the norm of a vector (i.e., e^X is not unitary).
    statement_B_correct = not is_unitary(expm(X))

    # C) Z and X represent observables (Hermitian).
    statement_C_correct = is_hermitian(Z) and is_hermitian(X)

    # D) (e^X)*Y*(e^{-X}) represents a quantum state.
    y_is_state, _ = is_quantum_state(Y)
    transformed_matrix = expm(X) @ Y @ expm(-X)
    transformed_is_state, transformed_reason = is_quantum_state(transformed_matrix)
    statement_D_correct = y_is_state and transformed_is_state

    # Determine the correct option based on the analysis
    correct_option = None
    if statement_A_correct: correct_option = 'A'
    elif statement_B_correct: correct_option = 'B'
    elif statement_C_correct: correct_option = 'C'
    elif statement_D_correct: correct_option = 'D'
    
    # There might be multiple correct statements, let's find all
    all_correct_options = []
    if statement_A_correct: all_correct_options.append('A')
    if statement_B_correct: all_correct_options.append('B')
    if statement_C_correct: all_correct_options.append('C')
    if statement_D_correct: all_correct_options.append('D')

    if len(all_correct_options) > 1:
        return f"Incorrect. The question implies a single correct statement, but the code found multiple correct options: {all_correct_options}."
    if len(all_correct_options) == 0:
        return "Incorrect. The code found that none of the statements are correct."

    # Compare with the LLM's answer
    if correct_option == llm_answer:
        return "Correct"
    else:
        # Build a detailed reason for the failure
        reason = f"Incorrect. The provided answer is {llm_answer}, but the code determined the correct answer is {correct_option}.\n"
        reason += "Here is the breakdown:\n"
        reason += f"A) W and X are unitary? W is unitary: {is_unitary(W)}. X is unitary: {is_unitary(X)}. Statement is {statement_A_correct}.\n"
        reason += f"B) e^X is not unitary? e^X is unitary: {is_unitary(expm(X))}. Statement is {statement_B_correct}.\n"
        reason += f"C) Z and X are Hermitian? Z is Hermitian: {is_hermitian(Z)}. X is Hermitian: {not is_hermitian(X)}. Statement is {statement_C_correct}.\n"
        reason += f"D) (e^X)*Y*(e^{{-X}}) is a quantum state? This is {statement_D_correct}. Reason: {transformed_reason}"
        return reason

# Execute the check and print the result
print(check_correctness())