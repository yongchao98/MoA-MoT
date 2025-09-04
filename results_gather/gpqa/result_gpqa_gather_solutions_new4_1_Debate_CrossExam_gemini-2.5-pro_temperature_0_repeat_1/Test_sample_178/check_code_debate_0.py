import numpy as np
from scipy.linalg import expm

def check_correctness():
    """
    This function checks the correctness of the provided answer to the quantum mechanics matrix question.
    It defines the matrices and the required quantum mechanical properties, evaluates each statement,
    and compares the result with the given answer 'C'.
    """

    # Helper function to parse matrix strings like "(c11, c12; c21, c22)"
    def parse_matrix(s: str) -> np.ndarray:
        """Parses the specific matrix string format into a numpy array."""
        # Replace 'i' with 'j' for Python's complex number notation
        s = s.replace('i', 'j')
        # Remove parentheses and spaces
        s = s.strip()[1:-1].replace(" ", "")
        rows = s.split(';')
        matrix = []
        for row_str in rows:
            if not row_str: continue
            elements = [complex(x) for x in row_str.split(',')]
            matrix.append(elements)
        return np.array(matrix, dtype=complex)

    # Helper functions to check matrix properties based on quantum mechanics definitions
    def is_hermitian(m, tol=1e-9):
        """Checks if a matrix is Hermitian (A = A†)."""
        return np.allclose(m, m.conj().T, atol=tol)

    def is_anti_hermitian(m, tol=1e-9):
        """Checks if a matrix is anti-Hermitian (A = -A†)."""
        return np.allclose(m, -m.conj().T, atol=tol)

    def is_unitary(m, tol=1e-9):
        """Checks if a matrix is unitary (A * A† = I)."""
        identity = np.identity(m.shape[0], dtype=complex)
        return np.allclose(m @ m.conj().T, identity, atol=tol)

    def is_positive_semidefinite(m, tol=1e-9):
        """Checks if a Hermitian matrix is positive semi-definite (all eigenvalues >= 0)."""
        if not is_hermitian(m, tol):
            return False
        # eigvalsh is efficient and numerically stable for Hermitian matrices
        eigenvalues = np.linalg.eigvalsh(m)
        return np.all(eigenvalues >= -tol)

    def is_density_matrix(m, tol=1e-9):
        """Checks if a matrix is a valid density matrix (Hermitian, trace=1, positive semi-definite)."""
        if not is_hermitian(m, tol):
            return (False, "it is not Hermitian")
        # The trace of a complex matrix can be complex, but for a Hermitian matrix it must be real.
        if not np.isclose(np.trace(m).real, 1.0, atol=tol) or not np.isclose(np.trace(m).imag, 0.0, atol=tol):
            return (False, f"its trace is {np.trace(m):.4f}, not 1")
        if not is_positive_semidefinite(m, tol):
            eigenvalues = np.linalg.eigvalsh(m)
            return (False, f"it is not positive semi-definite (eigenvalues: {[f'{v:.4f}' for v in eigenvalues]})")
        return (True, "it is a valid density matrix")

    # --- Matrix Definitions from the question ---
    W_str = "(0, 0, 1; 0, 1, 0; 1, 0, 0)"
    X_str = "(i, -1, 2i; 1, 0, 1;  2i, -1, -i)"
    Y_str = "(0.5, 0.1, 0.2; 0.1, 0.25, 0.1; 0.2, 0.1, 0.25)"
    Z_str = "(3, 2i, 5; -2i, -2, -4i; 5, 4i, 4)"
    
    # The final answer provided by the LLM to be checked
    provided_answer = 'C'

    try:
        W = parse_matrix(W_str)
        X = parse_matrix(X_str)
        Y = parse_matrix(Y_str)
        Z = parse_matrix(Z_str)
    except Exception as e:
        return f"Failed to parse matrices. Error: {e}"

    # --- Evaluate all statements based on quantum mechanics principles ---
    results = {}

    # A) Z and X represent observables. (Both must be Hermitian)
    results['A'] = is_hermitian(Z) and is_hermitian(X)
    
    # B) There exists a vector to which if one multiplies e^X, the norm of the vector changes. 
    # (This is true if e^X is NOT unitary. e^X is unitary iff X is anti-Hermitian)
    results['B'] = not is_anti_hermitian(X)

    # C) (e^X)*Y*(e^{-X}) represents a quantum state.
    # (This is true if Y is a density matrix and the transformation UYU† preserves its properties,
    # which happens if U=e^X is unitary, i.e., X is anti-Hermitian)
    y_is_dm, _ = is_density_matrix(Y)
    results['C'] = y_is_dm and is_anti_hermitian(X)

    # D) W and X represent the evolution operator of some quantum system. (Both must be unitary)
    results['D'] = is_unitary(W) and is_unitary(X)

    # --- Check the provided answer ('C') against the evaluation ---
    
    if not results[provided_answer]:
        # The provided answer corresponds to a statement that our code found to be false.
        # We need to explain why it's false.
        reason = f"Statement {provided_answer} is false because "
        if provided_answer == 'A':
            reason += "X is not Hermitian."
        elif provided_answer == 'B':
            reason += "X is anti-Hermitian, which means e^X is unitary and preserves vector norms."
        elif provided_answer == 'C':
            y_is_dm_check, y_reason = is_density_matrix(Y)
            x_is_anti_herm_check = is_anti_hermitian(X)
            if not y_is_dm_check: reason += f"Y is not a valid density matrix ({y_reason}) "
            if not x_is_anti_herm_check: reason += "and X is not anti-Hermitian."
        elif provided_answer == 'D':
            reason += "X is not unitary."
        
        return f"Incorrect. The provided answer is '{provided_answer}', but this statement is false. {reason.strip()}"

    # The provided answer corresponds to a true statement.
    # Now, check if any other statement is also true, which would make the question ambiguous.
    true_statements = [s for s, is_true in results.items() if is_true]
    
    if len(true_statements) > 1:
        return f"Incorrect. The provided answer '{provided_answer}' corresponds to a true statement, but other statements are also true: {true_statements}. The question is ambiguous as it asks to choose one correct statement."

    if len(true_statements) == 1 and true_statements[0] == provided_answer:
        return "Correct"
    else:
        # This case means the provided answer is not the single correct one.
        correct_answer = true_statements[0] if true_statements else "None"
        return f"Incorrect. The provided answer is '{provided_answer}', but the only correct statement is '{correct_answer}'."

# Execute the checker function and print the result
print(check_correctness())