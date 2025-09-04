import numpy as np

def check_correctness():
    """
    This function checks the correctness of the given answer to the quantum mechanics problem.

    The problem involves:
    1. Normalizing the state |alpha> proportional to (1+i)|up> + (2-i)|down>.
    2. Finding the eigenvalues and eigenstates of the operator A, where Aij = hbar/2 for i!=j and 0 otherwise.
    3. Calculating the probability of measuring the system in each eigenstate of A.
    4. Calculating the expectation value of A.
    5. Comparing the results with the provided option C: Probs = [0.64, 0.36], Avg = hbar/7.
    """
    # For numerical calculations, we can set hbar = 1 and add it back at the end.
    hbar = 1.0

    # --- 1. Normalize the initial state ---
    # The unnormalized state vector is v = [1+i, 2-i]^T
    alpha_unnormalized = np.array([1 + 1j, 2 - 1j])

    # The squared norm is <v|v> = (1-i)(1+i) + (2+i)(2-i) = (1+1) + (4+1) = 7
    norm_squared = np.vdot(alpha_unnormalized, alpha_unnormalized).real
    
    # The normalization constant is 1/sqrt(7)
    normalization_factor = 1 / np.sqrt(norm_squared)
    
    # The normalized state |alpha>
    alpha = normalization_factor * alpha_unnormalized

    # --- 2. Define the operator A and find its eigen-system ---
    # A = (hbar/2) * [[0, 1], [1, 0]]
    A = (hbar / 2) * np.array([[0, 1], [1, 0]])

    # Find eigenvalues and eigenvectors of A
    # The function np.linalg.eigh is suitable for Hermitian matrices and returns
    # real eigenvalues and a matrix of corresponding orthonormal eigenvectors.
    eigenvalues, eigenvectors = np.linalg.eigh(A)

    # Eigenvalues are lambda_1 = +hbar/2 and lambda_2 = -hbar/2
    # Eigenvectors are |a_1> and |a_2>
    
    # --- 3. Calculate measurement probabilities ---
    # The probability of measuring eigenvalue lambda_k is |<a_k|alpha>|^2
    
    # The eigenvectors are the columns of the 'eigenvectors' matrix.
    eigenstate_1 = eigenvectors[:, 0]
    eigenstate_2 = eigenvectors[:, 1]

    # Calculate the projection of |alpha> onto each eigenstate
    projection_1 = np.vdot(eigenstate_1, alpha)
    projection_2 = np.vdot(eigenstate_2, alpha)

    # The probabilities are the squared magnitudes of these projections
    prob_1 = np.abs(projection_1)**2
    prob_2 = np.abs(projection_2)**2

    # --- 4. Calculate the average (expectation) value ---
    # <A> = <alpha|A|alpha>
    avg_A = np.vdot(alpha, A @ alpha).real

    # --- 5. Check against the provided answer (Option C) ---
    # Option C states: Probabilities are 0.64, 0.36 and average value is hbar/7.
    
    # Check probabilities
    # The calculated probabilities should be compared with the given ones.
    # We sort both lists to ensure a consistent comparison.
    calculated_probs = sorted([prob_1, prob_2], reverse=True)
    answer_probs = sorted([0.64, 0.36], reverse=True)

    # We use np.allclose to account for floating point inaccuracies and rounding in the answer.
    # The answer seems rounded to two decimal places, so a tolerance of 0.005 is appropriate.
    if not np.allclose(calculated_probs, answer_probs, atol=0.005):
        # The exact probabilities are 9/14 (~0.6428) and 5/14 (~0.3571)
        exact_probs = sorted([9/14, 5/14], reverse=True)
        return (f"Incorrect probabilities. The answer provides {answer_probs}, but the "
                f"calculated probabilities are approximately {[round(p, 4) for p in calculated_probs]}. "
                f"The exact values are {[f'{p.numerator}/{p.denominator}' for p in [Fraction(9,14), Fraction(5,14)]]}.")

    # Check average value
    # The answer states the average value is hbar/7.
    expected_avg_A = hbar / 7.0
    
    if not np.isclose(avg_A, expected_avg_A):
        return (f"Incorrect average value. The answer states hbar/7, which is approximately "
                f"{expected_avg_A:.4f}*hbar. The calculated value is {avg_A:.4f}*hbar.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
# We need Fraction for a more descriptive error message, so import it inside the main block.
from fractions import Fraction
result = check_correctness()
print(result)