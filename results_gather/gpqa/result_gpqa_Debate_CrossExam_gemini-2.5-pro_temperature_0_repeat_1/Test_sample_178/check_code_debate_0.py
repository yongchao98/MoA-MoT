import numpy as np
from scipy.linalg import expm

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer by:
    1. Defining the matrices from the problem statement.
    2. Defining helper functions to check for quantum mechanical properties (Hermitian, Unitary, Density Matrix).
    3. Systematically evaluating each statement (A, B, C, D).
    4. Comparing the evaluation with the LLM's reasoning and final answer.
    5. Returning "Correct" if the answer is correct, or a reason if it is incorrect.
    """

    # --- 1. Define Matrices ---
    # Helper to parse matrix strings, handling 'i' for imaginary unit
    def parse_matrix(s):
        s = s.replace('i', 'j')
        rows = s.split(';')
        matrix_list = [[complex(x) for x in row.strip().split(',')] for row in rows]
        return np.array(matrix_list)

    W = parse_matrix("0, 0, 1; 0, 1, 0; 1, 0, 0")
    X = parse_matrix("i, -1, 2i; 1, 0, 1; 2i, -1, -i")
    Y = parse_matrix("0.5, 0.1, 0.2; 0.1, 0.25, 0.1; 0.2, 0.1, 0.25")
    Z = parse_matrix("3, 2i, 5; -2i, -2, -4i; 5, 4i, 4")

    # --- 2. Define Property Checkers ---
    TOL = 1e-9  # Tolerance for floating point comparisons

    def is_hermitian(matrix):
        return np.allclose(matrix, matrix.conj().T, atol=TOL)

    def is_anti_hermitian(matrix):
        return np.allclose(matrix, -matrix.conj().T, atol=TOL)

    def is_unitary(matrix):
        identity = np.identity(matrix.shape[0])
        return np.allclose(matrix.conj().T @ matrix, identity, atol=TOL)

    def is_positive_semidefinite(matrix):
        if not is_hermitian(matrix):
            return False
        # Use eigvalsh for Hermitian matrices for stability and real eigenvalues
        eigenvalues = np.linalg.eigvalsh(matrix)
        return np.all(eigenvalues >= -TOL)

    def is_density_matrix(matrix):
        if not is_hermitian(matrix):
            return f"it is not Hermitian."
        if not np.isclose(np.trace(matrix).real, 1.0, atol=TOL):
            return f"its trace is {np.trace(matrix).real}, not 1."
        if not is_positive_semidefinite(matrix):
            eigenvalues = np.linalg.eigvalsh(matrix)
            return f"it is not positive semi-definite (eigenvalues are {np.round(eigenvalues, 3)})."
        return None # No error

    # --- 3. Evaluate Statements ---

    # Statement B: W and X represent the evolution operator (must be unitary).
    w_is_unitary = is_unitary(W)
    x_is_unitary = is_unitary(X)
    statement_b_is_true = w_is_unitary and x_is_unitary
    
    # Check LLM reasoning for B
    if not w_is_unitary:
        return "Incorrect. The reasoning is flawed because W is not a unitary matrix, contrary to the LLM's finding."
    if x_is_unitary:
        return "Incorrect. The reasoning is flawed because X is a unitary matrix, contrary to the LLM's finding."
    if statement_b_is_true:
        return "Incorrect. Statement B is true, so A cannot be the unique correct answer."

    # Statement C: Z and X represent observables (must be Hermitian).
    z_is_hermitian = is_hermitian(Z)
    x_is_hermitian = is_hermitian(X)
    statement_c_is_true = z_is_hermitian and x_is_hermitian

    # Check LLM reasoning for C
    if not z_is_hermitian:
        return "Incorrect. The reasoning is flawed because Z is not a Hermitian matrix, contrary to the LLM's finding."
    if x_is_hermitian:
        return "Incorrect. The reasoning is flawed because X is a Hermitian matrix, contrary to the LLM's finding."
    if statement_c_is_true:
        return "Incorrect. Statement C is true, so A cannot be the unique correct answer."

    # Statement D: e^X can change a vector's norm (i.e., e^X is NOT unitary).
    # This relies on X being anti-Hermitian.
    x_is_anti_herm = is_anti_hermitian(X)
    if not x_is_anti_herm:
        return "Incorrect. The reasoning is flawed because X is not an anti-Hermitian matrix, which is a key step in the LLM's logic for statement D."
    
    e_X = expm(X)
    e_X_is_unitary = is_unitary(e_X)
    statement_d_is_true = not e_X_is_unitary

    # Check LLM reasoning for D
    if not e_X_is_unitary:
        return "Incorrect. The reasoning is flawed. Since X is anti-Hermitian, e^X must be unitary. However, the calculation shows it is not, indicating a potential numerical issue or a flaw in the premise."
    if statement_d_is_true:
        return "Incorrect. Statement D is true, so A cannot be the unique correct answer."

    # Statement A: (e^X)*Y*(e^{-X}) represents a quantum state (density matrix).
    # This is a unitary transformation of Y. It holds if Y is a density matrix.
    y_check_result = is_density_matrix(Y)
    if y_check_result is not None:
        return f"Incorrect. The answer A is wrong because the expression does not represent a quantum state. The matrix Y is not a valid density matrix because {y_check_result}"
    
    # If Y is a density matrix and the transformation is unitary, A is true.
    statement_a_is_true = (y_check_result is None) and e_X_is_unitary

    # --- 4. Final Verdict ---
    # The LLM chose A. This is correct if A is true and B, C, D are false.
    if statement_a_is_true and not statement_b_is_true and not statement_c_is_true and not statement_d_is_true:
        return "Correct"
    else:
        # This case handles if the final conclusion is wrong for any reason.
        return f"Incorrect. The evaluation of the statements is: A={statement_a_is_true}, B={statement_b_is_true}, C={statement_c_is_true}, D={statement_d_is_true}. The provided answer 'A' is not the sole correct statement."

# Run the check and print the result.
result = check_correctness()
print(result)