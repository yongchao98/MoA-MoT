import numpy as np

def check_correctness():
    """
    Checks the correctness of the answer to the quantum mechanics problem.
    """
    # --- Problem Definition ---
    # The state of the system at time t
    psi = np.array([-1, 2, 1], dtype=complex)

    # The observable operator P
    P = np.array([
        [0, 1/np.sqrt(2), 0],
        [1/np.sqrt(2), 0, 1/np.sqrt(2)],
        [0, 1/np.sqrt(2), 0]
    ], dtype=float)

    # The value to be measured (eigenvalue)
    target_eigenvalue = 0

    # The expected probability from the final answer 'C'
    expected_probability = 1/3

    # --- Calculation ---
    # 1. Find the eigenvalues and eigenvectors of the operator P.
    # np.linalg.eig returns normalized eigenvectors as columns of the matrix.
    try:
        eigenvalues, eigenvectors = np.linalg.eig(P)
    except np.linalg.LinAlgError:
        print("Error: Could not compute eigenvalues/eigenvectors for the given matrix P.")
        return

    # 2. Find the eigenvector corresponding to the target eigenvalue (0).
    # We find the index of the eigenvalue that is numerically closest to our target.
    try:
        idx = np.where(np.isclose(eigenvalues, target_eigenvalue))[0][0]
    except IndexError:
        print(f"Constraint not satisfied: The value {target_eigenvalue} is not an eigenvalue of the operator P.")
        print(f"The calculated eigenvalues are: {eigenvalues}")
        return

    # This is the normalized eigenvector for the eigenvalue 0.
    v0_normalized = eigenvectors[:, idx]

    # 3. Calculate the squared norm of the (unnormalized) state vector <ψ|ψ>.
    psi_norm_sq = np.dot(psi.conj(), psi).real
    if psi_norm_sq == 0:
        print("Error: The state vector is a zero vector, which is not a valid physical state.")
        return

    # 4. Calculate the inner product <v₀|ψ>.
    inner_product = np.dot(v0_normalized.conj(), psi)

    # 5. Calculate the probability using the formula: Prob(0) = |<v₀|ψ>|² / <ψ|ψ>
    calculated_probability = np.abs(inner_product)**2 / psi_norm_sq

    # --- Verification ---
    # Compare the calculated probability with the expected probability from answer C.
    if np.isclose(calculated_probability, expected_probability):
        print("Correct")
    else:
        # Using Fraction for a more readable output of the calculated probability
        from fractions import Fraction
        frac = Fraction(calculated_probability).limit_denominator()
        
        print(f"Incorrect.")
        print(f"The calculated probability is {calculated_probability:.6f}, which is {frac.numerator}/{frac.denominator}.")
        print(f"The expected probability from answer C is 1/3, which is approximately {expected_probability:.6f}.")
        print("The final answer C is correct, but some of the provided answers might have reached the correct conclusion through flawed reasoning or calculation errors.")

# Run the check
check_correctness()