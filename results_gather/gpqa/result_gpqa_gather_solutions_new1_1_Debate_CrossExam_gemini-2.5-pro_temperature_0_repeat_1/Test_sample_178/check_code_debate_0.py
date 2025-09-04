import numpy as np
from scipy.linalg import expm

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics matrix question.
    """
    # Define the matrices from the question
    W = np.array([[0, 0, 1], 
                  [0, 1, 0], 
                  [1, 0, 0]], dtype=complex)

    X = np.array([[1j, -1, 2j], 
                  [1, 0, 1],  
                  [2j, -1, -1j]], dtype=complex)

    Y = np.array([[0.5, 0.1, 0.2], 
                  [0.1, 0.25, 0.1], 
                  [0.2, 0.1, 0.25]], dtype=complex)

    Z = np.array([[3, 2j, 5], 
                  [-2j, -2, -4j], 
                  [5, 4j, 4]], dtype=complex)

    # --- Helper functions for quantum mechanical properties ---
    def is_hermitian(M, tol=1e-9):
        """Checks if a matrix is Hermitian."""
        return np.allclose(M, M.conj().T, atol=tol)

    def is_skew_hermitian(M, tol=1e-9):
        """Checks if a matrix is skew-Hermitian."""
        return np.allclose(M, -M.conj().T, atol=tol)

    def is_unitary(M, tol=1e-9):
        """Checks if a matrix is unitary."""
        identity = np.identity(M.shape[0])
        return np.allclose(M.conj().T @ M, identity, atol=tol)

    def is_density_matrix(M, tol=1e-9):
        """Checks if a matrix is a valid density matrix."""
        # 1. Must be Hermitian
        if not is_hermitian(M, tol):
            return False, "is not Hermitian"
        
        # 2. Trace must be 1
        if not np.isclose(np.trace(M).real, 1.0, atol=tol):
            return False, f"has a trace of {np.trace(M).real}, not 1"
        
        # 3. Must be positive semi-definite (all eigenvalues >= 0)
        # Use eigvalsh for Hermitian matrices to get real eigenvalues
        eigenvalues = np.linalg.eigvalsh(M)
        if not np.all(eigenvalues >= -tol):
            return False, f"is not positive semi-definite (eigenvalues: {eigenvalues})"
            
        return True, "is a valid density matrix"

    # --- Evaluate each statement from the question ---
    
    # Statement A: W and X represent the evolution operator of some quantum system.
    # Constraint: Both W and X must be unitary.
    w_is_unitary = is_unitary(W)
    x_is_unitary = is_unitary(X)
    statement_A_is_true = w_is_unitary and x_is_unitary
    
    # Statement B: Z and X represent observables.
    # Constraint: Both Z and X must be Hermitian.
    z_is_hermitian = is_hermitian(Z)
    x_is_hermitian = is_hermitian(X)
    statement_B_is_true = z_is_hermitian and x_is_hermitian

    # Statement C: There exists a vector to which if one multiplies e^X, the norm of the vector changes.
    # Constraint: e^X must NOT be unitary. This is true if X is NOT skew-Hermitian.
    x_is_skew_hermitian = is_skew_hermitian(X)
    # If X is skew-Hermitian, e^X is unitary, and the statement is false.
    statement_C_is_true = not x_is_skew_hermitian

    # Statement D: (e^X)*Y*(e^{-X}) represents a quantum state.
    # Constraint 1: Y must be a valid density matrix.
    # Constraint 2: The transformation must preserve the density matrix properties.
    # This is true if the transformation is unitary, i.e., e^X is unitary.
    # e^X is unitary if X is skew-Hermitian.
    y_is_density, y_reason = is_density_matrix(Y)
    # We already checked if X is skew-Hermitian for statement C.
    statement_D_is_true = y_is_density and x_is_skew_hermitian

    # --- Final check against the provided answer "D" ---
    correct_statements = []
    if statement_A_is_true: correct_statements.append("A")
    if statement_B_is_true: correct_statements.append("B")
    if statement_C_is_true: correct_statements.append("C")
    if statement_D_is_true: correct_statements.append("D")

    # The provided answer is D. We check if D is the *only* correct statement.
    if len(correct_statements) == 1 and correct_statements[0] == "D":
        return "Correct"
    else:
        # Provide a detailed reason for the error.
        reasons = []
        if not statement_D_is_true:
            if not y_is_density:
                reasons.append(f"Statement D is false because Y {y_reason}.")
            if not x_is_skew_hermitian:
                reasons.append("Statement D is false because the evolution operator e^X is not unitary (since X is not skew-Hermitian).")
        
        if statement_A_is_true:
            reasons.append("Statement A is true, but the answer was D.")
        if statement_B_is_true:
            reasons.append("Statement B is true, but the answer was D.")
        if statement_C_is_true:
            reasons.append("Statement C is true, but the answer was D.")
            
        if len(correct_statements) > 1:
             reasons.append(f"Multiple statements are correct: {correct_statements}. The answer should be unique.")
        elif not correct_statements:
             reasons.append("No statement is correct.")

        return f"Incorrect. The provided answer 'D' is wrong. Reason(s): {' '.join(reasons)}"

# Run the check
print(check_answer())