import numpy as np

def solve_quantum_fidelity():
    """
    Solves the quantum fidelity problem based on the hypothetical universe's rules.
    """
    # Step 1: Define the quantum system
    O = np.array([[3, 1], [1, 2]])
    psi_initial = np.array([np.sqrt(3)/2, 1/2])

    print("Observable Operator O:\n", O)
    print("\nInitial State |ψ⟩:\n", psi_initial, "\n")

    # Step 2: Find the eigenvalues and eigenvectors of the observable O.
    # np.linalg.eigh sorts eigenvalues in ascending order for Hermitian matrices.
    eigenvalues, eigenvectors = np.linalg.eigh(O)

    # In this 2D case, the second-largest eigenvalue is the smaller one.
    lambda_2 = eigenvalues[0]
    v_2 = eigenvectors[:, 0]  # This is the target state from Rule 1.

    # The largest eigenvalue.
    lambda_1 = eigenvalues[1]
    v_1 = eigenvectors[:, 1]

    # Step 3: Project the initial state onto the eigenbasis to find the coefficients.
    # |ψ⟩ = c₁|v₁⟩ + c₂|v₂⟩
    # Since eigenvectors from eigh are orthonormal and all our vectors are real,
    # c₁ = ⟨v₁|ψ⟩ = v₁ᵀ • |ψ⟩ and c₂ = ⟨v₂|ψ⟩ = v₂ᵀ • |ψ⟩.
    c_1 = np.dot(v_1, psi_initial)
    c_2 = np.dot(v_2, psi_initial)

    # Step 4: Calculate the components of the fidelity formula.
    # The fidelity F between the resulting state (from Rule 2) and the target
    # state |v₂⟩ (from Rule 1) is given by:
    # F = |c₂λ₂³|² / (|c₁λ₁³|² + |c₂λ₂³|²) = (c₂² * λ₂⁶) / (c₁² * λ₁⁶ + c₂² * λ₂⁶)

    # We will now calculate each term of this equation.
    c1_sq_lambda1_pow6 = (c_1**2) * (lambda_1**6)
    c2_sq_lambda2_pow6 = (c_2**2) * (lambda_2**6)

    # Step 5: Compute the final fidelity.
    denominator = c1_sq_lambda1_pow6 + c2_sq_lambda2_pow6
    if denominator == 0:
        fidelity = 0.0 # Avoid division by zero, though unlikely here.
    else:
        fidelity = c2_sq_lambda2_pow6 / denominator

    # Step 6: Print the results and the breakdown of the final equation.
    print("The fidelity F is calculated using the formula: F = (c₂² * λ₂⁶) / (c₁² * λ₁⁶ + c₂² * λ₂⁶)\n")

    print(f"The eigenvalues are λ₁ (largest) = {lambda_1:.5f} and λ₂ (second-largest) = {lambda_2:.5f}.")
    print(f"The initial state projected onto the eigenbasis has coefficients c₁ = {c_1:.5f} and c₂ = {c_2:.5f}.\n")

    print("The terms of the final equation are:")
    print(f"Numerator term (c₂² * λ₂⁶) = ({c_2:.5f})² * ({lambda_2:.5f})⁶ = {c2_sq_lambda2_pow6:.8f}")
    print(f"First denominator term (c₁² * λ₁⁶) = ({c_1:.5f})² * ({lambda_1:.5f})⁶ = {c1_sq_lambda1_pow6:.8f}")
    print(f"Second denominator term is the same as the numerator: {c2_sq_lambda2_pow6:.8f}\n")

    print("Final fidelity calculation:")
    print(f"F = {c2_sq_lambda2_pow6:.8f} / ({c1_sq_lambda1_pow6:.8f} + {c2_sq_lambda2_pow6:.8f})")
    print(f"F = {c2_sq_lambda2_pow6:.8f} / {denominator:.8f}\n")

    print("The final fidelity is:")
    print(f"<<<{fidelity}>>>")


if __name__ == '__main__':
    solve_quantum_fidelity()
