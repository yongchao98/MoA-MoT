import numpy as np
from scipy.linalg import expm

def check_correctness():
    """
    This function verifies the correct statement about the given matrices based on quantum mechanics principles.
    
    It checks each statement (A, B, C, D) by:
    1. Defining the matrices W, X, Y, Z.
    2. Implementing helper functions to check for matrix properties:
       - is_hermitian: M == M.conj().T (for observables)
       - is_anti_hermitian: M.conj().T == -M (for generators of unitary evolution)
       - is_unitary: M.conj().T @ M == I (for evolution operators)
       - is_density_matrix: Hermitian, trace=1, positive semi-definite (for quantum states)
    3. Systematically evaluating each statement based on these properties.
    4. Comparing the logically derived correct statement with the provided answer ('C').
    """
    
    # 1. Define the matrices from the problem description
    try:
        W = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=float)
        X = np.array([[1j, -1, 2j], [1, 0, 1], [2j, -1, -1j]], dtype=complex)
        Y = np.array([[0.5, 0.1, 0.2], [0.1, 0.25, 0.1], [0.2, 0.1, 0.25]], dtype=float)
        Z = np.array([[3, 2j, 5], [-2j, -2, -4j], [5, 4j, 4]], dtype=complex)
    except Exception as e:
        return f"Failed to initialize matrices. Error: {e}"

    # 2. Helper functions to check matrix properties
    def is_hermitian(matrix):
        return np.allclose(matrix, matrix.conj().T)

    def is_anti_hermitian(matrix):
        return np.allclose(matrix.conj().T, -matrix)

    def is_unitary(matrix):
        identity = np.identity(matrix.shape[0])
        return np.allclose(matrix.conj().T @ matrix, identity)

    def is_density_matrix(matrix):
        # Must be Hermitian
        if not is_hermitian(matrix):
            return False, "it is not Hermitian"
        # Must have trace = 1
        if not np.isclose(np.trace(matrix).real, 1.0):
            return False, f"its trace is {np.trace(matrix).real:.2f}, not 1"
        # Must be positive semi-definite (all eigenvalues >= 0)
        # np.linalg.eigvalsh is used for Hermitian matrices for better stability and real eigenvalues
        eigenvalues = np.linalg.eigvalsh(matrix)
        if not np.all(eigenvalues >= -1e-9):  # Use tolerance for floating point errors
            return False, f"it is not positive semi-definite (eigenvalues: {np.round(eigenvalues, 3)})"
        return True, "it is a valid density matrix"

    # 3. Evaluate each statement from the question
    results = {}

    # Statement A: There exists a vector to which if one multiplies e^X, the norm of the vector changes.
    # This is true iff e^X is NOT unitary. e^X is unitary iff X is anti-Hermitian.
    if is_anti_hermitian(X):
        # e^X is unitary, so it preserves norms. The statement is false.
        results['A'] = (False, "X is anti-Hermitian, so e^X is a unitary operator which always preserves the norm of a vector.")
    else:
        results['A'] = (True, "X is not anti-Hermitian, so e^X is not unitary and would change a vector's norm.")

    # Statement B: Z and X represent observables.
    # This is true iff both Z and X are Hermitian.
    is_Z_herm = is_hermitian(Z)
    is_X_herm = is_hermitian(X)
    if is_Z_herm and is_X_herm:
        results['B'] = (True, "Both Z and X are Hermitian.")
    else:
        reason = []
        if not is_Z_herm: reason.append("Z is not Hermitian.")
        if not is_X_herm: reason.append("X is not Hermitian.")
        results['B'] = (False, " ".join(reason))

    # Statement C: (e^X)*Y*(e^{-X}) represents a quantum state.
    # This is true if Y is a density matrix and the transformation is unitary.
    # The transformation e^X is unitary because X is anti-Hermitian.
    # A unitary transformation of a valid density matrix results in another valid density matrix.
    is_Y_dm, y_reason = is_density_matrix(Y)
    is_X_anti_herm = is_anti_hermitian(X)
    if is_Y_dm and is_X_anti_herm:
        results['C'] = (True, "Y is a valid density matrix, and since X is anti-Hermitian, e^X is unitary. A unitary transformation of a density matrix results in another valid density matrix.")
    else:
        reason = []
        if not is_Y_dm: reason.append(f"Y is not a valid density matrix because {y_reason}.")
        if not is_X_anti_herm: reason.append("The transformation operator e^X is not unitary because X is not anti-Hermitian.")
        results['C'] = (False, " ".join(reason))

    # Statement D: W and X represent the evolution operator of some quantum system.
    # This is true iff both W and X are unitary.
    is_W_unitary = is_unitary(W)
    is_X_unitary = is_unitary(X)
    if is_W_unitary and is_X_unitary:
        results['D'] = (True, "Both W and X are unitary.")
    else:
        reason = []
        if not is_W_unitary: reason.append("W is not unitary.")
        if not is_X_unitary: reason.append("X is not unitary.")
        results['D'] = (False, " ".join(reason))

    # 4. Determine the correct statement and compare with the LLM's answer
    llm_answer = 'C'
    
    true_statements = [s for s, (is_true, _) in results.items() if is_true]

    if len(true_statements) == 0:
        return "Incorrect. The provided answer is C, but the code found that no statement is true."
    
    if len(true_statements) > 1:
        return f"Incorrect. The provided answer is C, but the code found multiple true statements: {', '.join(true_statements)}. This indicates a flaw in the question's design."

    correct_statement = true_statements[0]

    if correct_statement == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer}, but the code determined that statement {correct_statement} is the correct one.\n"
                f"Reason for {correct_statement} being correct: {results[correct_statement][1]}\n"
                f"Reason for {llm_answer} being incorrect: {results[llm_answer][1]}")

# Execute the check and print the result
result = check_correctness()
print(result)