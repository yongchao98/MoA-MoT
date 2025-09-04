import numpy as np
from scipy.linalg import expm

def check_correctness():
    """
    This function defines the matrices from the problem and checks the validity of each statement (A, B, C, D)
    to determine if the provided answer 'A' is correct.
    """
    # Define the matrices as complex numpy arrays
    W = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=complex)
    X = np.array([[1j, -1, 2j], [1, 0, 1], [2j, -1, -1j]], dtype=complex)
    Y = np.array([[0.5, 0.1, 0.2], [0.1, 0.25, 0.1], [0.2, 0.1, 0.25]], dtype=complex)
    Z = np.array([[3, 2j, 5], [-2j, -2, -4j], [5, 4j, 4]], dtype=complex)

    # --- Helper functions to check matrix properties ---
    def is_hermitian(matrix):
        # A matrix is Hermitian if it equals its conjugate transpose.
        return np.allclose(matrix, matrix.conj().T)

    def is_skew_hermitian(matrix):
        # A matrix is skew-Hermitian if it equals the negative of its conjugate transpose.
        return np.allclose(matrix, -matrix.conj().T)

    def is_unitary(matrix):
        # A matrix is unitary if its product with its conjugate transpose is the identity matrix.
        identity = np.identity(matrix.shape[0])
        return np.allclose(matrix @ matrix.conj().T, identity)

    def is_density_matrix(matrix):
        # Check the three conditions for a density matrix.
        # 1. Hermitian
        if not is_hermitian(matrix):
            return False, "is not Hermitian"
        # 2. Trace is 1
        if not np.isclose(np.trace(matrix), 1.0):
            return False, f"has a trace of {np.trace(matrix).real:.2f}, not 1"
        # 3. Positive semi-definite (non-negative eigenvalues)
        # Use eigvalsh for Hermitian matrices for better stability and real eigenvalues
        eigenvalues = np.linalg.eigvalsh(matrix)
        if not np.all(eigenvalues >= -1e-9): # Use tolerance for floating point errors
            return False, "is not positive semi-definite (has negative eigenvalues)"
        return True, ""

    # --- Evaluate each statement ---

    # B) Z and X represent observables.
    # This requires both Z and X to be Hermitian.
    z_is_hermitian = is_hermitian(Z)
    x_is_hermitian = is_hermitian(X)
    statement_B_is_true = z_is_hermitian and x_is_hermitian
    if statement_B_is_true and "A" != "B":
        return "Incorrect. Statement B is also true because both Z and X are Hermitian."
    if not x_is_hermitian and "A" == "B":
        return "Incorrect. Statement B is false because X is not a Hermitian matrix."

    # D) W and X represent the evolution operator of some quantum system.
    # This requires both W and X to be unitary.
    w_is_unitary = is_unitary(W)
    x_is_unitary = is_unitary(X)
    statement_D_is_true = w_is_unitary and x_is_unitary
    if statement_D_is_true and "A" != "D":
        return "Incorrect. Statement D is also true because both W and X are unitary."
    if not x_is_unitary and "A" == "D":
        return "Incorrect. Statement D is false because X is not a unitary matrix."

    # C) There exists a vector to which if one multiplies e^X, the norm of the vector changes.
    # This is equivalent to saying e^X is NOT unitary.
    # e^X is unitary if and only if X is skew-Hermitian.
    x_is_skew_hermitian = is_skew_hermitian(X)
    statement_C_is_true = not x_is_skew_hermitian
    if statement_C_is_true and "A" != "C":
        return "Incorrect. Statement C is also true because e^X is not unitary."
    if not statement_C_is_true and "A" == "C":
        return "Incorrect. Statement C is false. X is skew-Hermitian, which means e^X is unitary and preserves the norm of any vector."

    # A) (e^X)*Y*(e^{-X}) represents a quantum state.
    # If X is skew-Hermitian, then e^X is unitary and e^{-X} = (e^X)†.
    # The operation is then a unitary transformation of Y.
    # A unitary transformation of a valid density matrix results in another valid density matrix.
    # So, we must check if Y is a valid density matrix and if the transformation is unitary.
    y_is_density, y_reason = is_density_matrix(Y)
    
    # We already established X is skew-Hermitian, so e^X is unitary.
    # Therefore, statement A is true if and only if Y is a density matrix.
    statement_A_is_true = y_is_density and x_is_skew_hermitian
    
    # --- Final Verification ---
    llm_answer = "A"
    
    # Check if the LLM's answer is correct and unique
    is_correct = (llm_answer == "A" and statement_A_is_true) and \
                 not statement_B_is_true and \
                 not statement_C_is_true and \
                 not statement_D_is_true

    if is_correct:
        return "Correct"
    
    # If the LLM answer is 'A' but it's false, explain why.
    if llm_answer == "A" and not statement_A_is_true:
        if not y_is_density:
            return f"Incorrect. Statement A is false because Y is not a valid quantum state (density matrix); it {y_reason}."
        if not x_is_skew_hermitian:
            return "Incorrect. Statement A is false because the transformation is not unitary (X is not skew-Hermitian)."

    # If another statement is true, the answer 'A' is not uniquely correct.
    true_statements = []
    if statement_A_is_true: true_statements.append("A")
    if statement_B_is_true: true_statements.append("B")
    if statement_C_is_true: true_statements.append("C")
    if statement_D_is_true: true_statements.append("D")

    if len(true_statements) > 1:
        return f"Incorrect. The answer 'A' is not the only correct option. Statements {', '.join(true_statements)} are all true."
    
    # General fallback for any other error
    return f"Incorrect. The correct statement is {true_statements[0]}, but the provided answer was {llm_answer}."


# Run the check
result = check_correctness()
print(result)