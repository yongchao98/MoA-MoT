import numpy as np
from scipy.linalg import expm

def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying the properties of the given matrices.
    """
    # Define helper functions for quantum mechanical properties
    def is_hermitian(M, tol=1e-9):
        """Checks if a matrix is Hermitian."""
        return np.allclose(M, M.conj().T, atol=tol)

    def is_anti_hermitian(M, tol=1e-9):
        """Checks if a matrix is anti-Hermitian (skew-Hermitian)."""
        return np.allclose(M, -M.conj().T, atol=tol)

    def is_unitary(M, tol=1e-9):
        """Checks if a matrix is unitary."""
        identity = np.identity(M.shape[0], dtype=complex)
        return np.allclose(M @ M.conj().T, identity, atol=tol)

    def is_density_matrix(M, tol=1e-9):
        """Checks if a matrix is a valid density matrix."""
        # 1. Must be Hermitian
        if not is_hermitian(M, tol):
            return False, "Matrix is not Hermitian."
        # 2. Trace must be 1
        if not np.isclose(np.trace(M).real, 1.0, atol=tol):
            return False, f"Trace is {np.trace(M)}, not 1."
        # 3. Must be positive semi-definite (all eigenvalues >= 0)
        # For a Hermitian matrix, eigenvalues are guaranteed to be real.
        eigenvalues = np.linalg.eigvalsh(M)
        if not np.all(eigenvalues >= -tol):
            return False, f"Matrix has negative eigenvalues: {eigenvalues}."
        return True, "Is a valid density matrix."

    # Define the matrices from the question
    W = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=complex)
    X = np.array([[1j, -1, 2j], [1, 0, 1], [2j, -1, -1j]], dtype=complex)
    Y = np.array([[0.5, 0.1, 0.2], [0.1, 0.25, 0.1], [0.2, 0.1, 0.25]], dtype=complex)
    Z = np.array([[3, 2j, 5], [-2j, -2, -4j], [5, 4j, 4]], dtype=complex)

    # --- Evaluate each statement ---
    statement_correctness = {}

    # A) There exists a vector to which if one multiplies e^X, the norm of the vector changes.
    # This is true iff e^X is NOT unitary. e^X is unitary iff X is anti-Hermitian.
    # A unitary operator preserves the norm of a vector.
    e_X = expm(X)
    statement_correctness['A'] = not is_unitary(e_X)

    # B) (e^X)*Y*(e^{-X}) represents a quantum state.
    # This is true if Y is a density matrix and the transformation UYU_dagger preserves its properties.
    # The transformation is UYU_dagger if U=e^X is unitary, which is true if X is anti-Hermitian.
    is_Y_dm, _ = is_density_matrix(Y)
    is_X_anti_h = is_anti_hermitian(X)
    statement_correctness['B'] = is_Y_dm and is_X_anti_h

    # C) W and X represent the evolution operator of some quantum system.
    # This is true iff both W and X are unitary.
    statement_correctness['C'] = is_unitary(W) and is_unitary(X)

    # D) Z and X represent observables.
    # This is true iff both Z and X are Hermitian.
    statement_correctness['D'] = is_hermitian(Z) and is_hermitian(X)

    # --- Compare with LLM's answer ---
    llm_answer = 'B'
    true_statements = [k for k, v in statement_correctness.items() if v]

    if len(true_statements) == 1 and true_statements[0] == llm_answer:
        return "Correct"
    else:
        # Construct a detailed reason for the incorrectness
        reason = f"The final answer is incorrect. The LLM chose '{llm_answer}', but the correct statement is '{true_statements[0] if true_statements else 'None'}'.\n\n"
        reason += "Here is a breakdown of each statement's validity:\n"
        
        # Reason for A
        reason += f"A) 'Norm changes for e^X': This is FALSE.\n"
        reason += f"   - A norm-preserving operator is unitary. The operator e^X is unitary because X is anti-Hermitian (X = -X†).\n"
        
        # Reason for B
        reason += f"B) '(e^X)*Y*(e^-X) is a quantum state': This is TRUE.\n"
        reason += f"   - Y is a valid density matrix (Hermitian, trace=1, positive semi-definite).\n"
        reason += f"   - The transformation (e^X)Y(e^-X) is a unitary transformation (since X is anti-Hermitian), which always maps a valid density matrix to another valid density matrix.\n"
        
        # Reason for C
        reason += f"C) 'W and X are evolution operators': This is FALSE.\n"
        reason += f"   - An evolution operator must be unitary. While W is unitary, X is not.\n"
        
        # Reason for D
        reason += f"D) 'Z and X are observables': This is FALSE.\n"
        reason += f"   - An observable must be Hermitian. While Z is Hermitian, X is not (it is anti-Hermitian).\n"
        
        reason += f"\nConclusion: The only true statement is B. The LLM's answer '{llm_answer}' is correct, but its reasoning led to the wrong letter in the final output. The provided final answer was <<<B>>>, which corresponds to the statement '(e^X)*Y*(e^{-X}) represents a quantum state'. This statement is indeed the correct one."
        
        # This is a special case where the LLM's reasoning is correct, and its final output letter is correct, but the prompt's final answer is different.
        # Let's re-read the prompt's final answer.
        # The prompt's final answer is <<<B>>>.
        # The analysis shows B is the correct statement.
        # So the prompt's final answer is correct.
        # The code should return "Correct".
        # The logic `if len(true_statements) == 1 and true_statements[0] == llm_answer:` handles this.
        # The `else` block is for when the LLM's answer is genuinely wrong.
        # Let's adjust the `else` block to be more general.
        
        final_reason = f"The LLM's answer '{llm_answer}' is incorrect.\n"
        final_reason += "Here is the evaluation of each statement:\n"
        for statement, is_correct in statement_correctness.items():
            final_reason += f" - Statement {statement} is {'TRUE' if is_correct else 'FALSE'}.\n"
        final_reason += f"The only correct statement is '{true_statements[0] if true_statements else 'None'}'."
        return final_reason

# Execute the check
print(check_answer())