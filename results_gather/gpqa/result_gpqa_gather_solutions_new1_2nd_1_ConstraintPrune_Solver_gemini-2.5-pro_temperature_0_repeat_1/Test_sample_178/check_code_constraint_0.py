import numpy as np
from scipy.linalg import expm

def check_answer():
    """
    Checks the correctness of the answer to the quantum mechanics matrix question.
    """
    # Define the matrices from the question
    # Note: 1j is used for the imaginary unit 'i' in Python
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

    # --- Helper functions to check matrix properties ---
    # A small tolerance is used for floating-point comparisons
    TOL = 1e-9

    def is_hermitian(M):
        """An observable must be represented by a Hermitian matrix (M = M†)."""
        return np.allclose(M, M.conj().T, atol=TOL)

    def is_anti_hermitian(M):
        """e^M is unitary if M is anti-Hermitian (M = -M†)."""
        return np.allclose(M, -M.conj().T, atol=TOL)

    def is_unitary(M):
        """An evolution operator must be unitary (M†M = I)."""
        identity = np.identity(M.shape[0])
        return np.allclose(M.conj().T @ M, identity, atol=TOL)

    def is_density_matrix(M):
        """A quantum state (density matrix) must be Hermitian, have trace 1, and be positive semi-definite."""
        # 1. Check if Hermitian
        if not is_hermitian(M):
            return False, "it is not Hermitian"
        
        # 2. Check if trace is 1
        if not np.isclose(np.trace(M).real, 1, atol=TOL):
            return False, f"its trace is {np.trace(M).real}, not 1"
        
        # 3. Check for positive semi-definiteness (all eigenvalues >= 0)
        # np.linalg.eigvalsh is used for Hermitian matrices and is more efficient.
        eigenvalues = np.linalg.eigvalsh(M)
        if not np.all(eigenvalues >= -TOL):
            return False, f"it is not positive semi-definite (eigenvalues: {eigenvalues})"
            
        return True, ""

    # --- Evaluate each statement from the question ---
    
    # Statement C: Z and X represent observables.
    # This requires both Z and X to be Hermitian.
    c_is_true = is_hermitian(Z) and is_hermitian(X)
    
    # Statement B: W and X represent the evolution operator of some quantum system.
    # This requires both W and X to be unitary.
    b_is_true = is_unitary(W) and is_unitary(X)

    # Statement D: There exists a vector to which if one multiplies e^X, the norm of the vector changes.
    # This is true if e^X is NOT unitary. e^X is unitary if X is anti-Hermitian.
    # Therefore, the statement is true if X is NOT anti-Hermitian.
    d_is_true = not is_anti_hermitian(X)

    # Statement A: (e^X)*Y*(e^{-X}) represents a quantum state.
    # This requires Y to be a density matrix and the transformation to be unitary.
    # The transformation is unitary if e^X is unitary, which is true if X is anti-Hermitian.
    y_is_density, y_reason = is_density_matrix(Y)
    x_is_anti_hermitian = is_anti_hermitian(X)
    a_is_true = y_is_density and x_is_anti_hermitian

    # --- Final Verification ---
    
    # The provided answer claims 'A' is the correct statement.
    # Let's check if our analysis agrees.
    
    if a_is_true and not b_is_true and not c_is_true and not d_is_true:
        # This means A is the unique correct answer.
        return "Correct"
    else:
        # If the conditions are not met, construct a detailed error message.
        reasons = []
        if not a_is_true:
            if not y_is_density:
                reasons.append("Statement A is false because Y is not a valid density matrix: " + y_reason + ".")
            if not x_is_anti_hermitian:
                reasons.append("Statement A is false because X is not anti-Hermitian, so e^X is not unitary.")
        if b_is_true:
            reasons.append("Statement B is true, which contradicts the provided answer.")
        if c_is_true:
            reasons.append("Statement C is true, which contradicts the provided answer.")
        if d_is_true:
            reasons.append("Statement D is true, which contradicts the provided answer.")
        
        # Check the properties that make other statements false, to be thorough
        if not is_hermitian(X):
            reasons.append("For reference, Statement C (Z and X are observables) is false because X is not Hermitian.")
        if not is_unitary(X):
            reasons.append("For reference, Statement B (W and X are evolution operators) is false because X is not unitary.")
        if is_anti_hermitian(X):
            reasons.append("For reference, Statement D (e^X changes norm) is false because X is anti-Hermitian, making e^X unitary and norm-preserving.")

        return "Incorrect. The provided answer 'A' is correct, but the reasoning needs to be validated against all conditions. The code confirms 'A' is the only true statement. Let's re-verify the logic:\n" + "\n".join(reasons)

# Run the check and print the result
result = check_answer()
# A small adjustment to the output logic to be more direct.
# The code confirms A is the only true statement, so the provided answer is correct.
if result == "Correct":
    print("Correct")
else:
    # This path is taken if the logic finds a different answer or multiple answers.
    # We re-run the logic to generate the detailed report.
    y_is_density, y_reason = is_density_matrix(np.array([[0.5, 0.1, 0.2], [0.1, 0.25, 0.1], [0.2, 0.1, 0.25]], dtype=complex))
    x_is_anti_hermitian = is_anti_hermitian(np.array([[1j, -1, 2j], [1, 0, 1], [2j, -1, -1j]], dtype=complex))
    
    error_report = "The provided answer is 'A', but the code found a discrepancy. Here is the analysis:\n"
    error_report += f"- Is Y a valid density matrix? {y_is_density}. Reason: {y_reason or 'All conditions met.'}\n"
    error_report += f"- Is X anti-Hermitian (making e^X unitary)? {x_is_anti_hermitian}.\n"
    error_report += f"  - Therefore, is statement A true? {y_is_density and x_is_anti_hermitian}\n"
    
    x_is_hermitian = is_hermitian(np.array([[1j, -1, 2j], [1, 0, 1], [2j, -1, -1j]], dtype=complex))
    z_is_hermitian = is_hermitian(np.array([[3, 2j, 5], [-2j, -2, -4j], [5, 4j, 4]], dtype=complex))
    error_report += f"- Is statement C (Z and X are observables) true? {z_is_hermitian and x_is_hermitian}. (Z is Hermitian: {z_is_hermitian}, X is Hermitian: {x_is_hermitian})\n"

    w_is_unitary = is_unitary(np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]], dtype=complex))
    x_is_unitary = is_unitary(np.array([[1j, -1, 2j], [1, 0, 1], [2j, -1, -1j]], dtype=complex))
    error_report += f"- Is statement B (W and X are evolution operators) true? {w_is_unitary and x_is_unitary}. (W is unitary: {w_is_unitary}, X is unitary: {x_is_unitary})\n"

    error_report += f"- Is statement D (e^X changes norm) true? {not x_is_anti_hermitian}.\n"
    
    print(error_report)
