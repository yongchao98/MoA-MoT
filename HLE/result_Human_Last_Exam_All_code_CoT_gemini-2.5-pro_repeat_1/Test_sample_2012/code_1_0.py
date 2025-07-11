import numpy as np

def solve_quantum_fidelity():
    """
    Calculates the fidelity based on the custom physical laws of Universe U.
    """
    # The observable operator
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 1: Find the eigenvalues of the observable operator.
    # Eigenvectors are also calculated but not directly needed for the final formula.
    eigenvalues, eigenvectors = np.linalg.eig(O)

    # Sort eigenvalues in descending order to identify the largest and second-largest.
    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_indices]

    lambda_1 = sorted_eigenvalues[0]  # Largest eigenvalue
    lambda_2 = sorted_eigenvalues[1]  # Second-largest eigenvalue

    # Both eigenvalues are positive, so their magnitudes are the values themselves.
    # The fidelity F = |λ₂|⁶ / (|λ₁|⁶ + |λ₂|⁶).
    lambda_1_pow_6 = lambda_1**6
    lambda_2_pow_6 = lambda_2**6

    fidelity = lambda_2_pow_6 / (lambda_1_pow_6 + lambda_2_pow_6)

    print("The observable operator is:")
    print(O)
    print("\nStep 1: Find the eigenvalues of the operator.")
    print(f"The eigenvalues are λ₁ (largest) = {lambda_1:.8f} and λ₂ (second-largest) = {lambda_2:.8f}.")

    print("\nStep 2: Determine the fidelity formula based on the laws of universe U.")
    print("The post-measurement state is a superposition with coefficients proportional to the cube of the eigenvalues' magnitudes.")
    print("The fidelity with the eigenstate of the second-largest eigenvalue simplifies to F = |λ₂|⁶ / (|λ₁|⁶ + |λ₂|⁶).")

    print("\nStep 3: Calculate the numerical values for the formula.")
    print(f"|λ₁|⁶ = {lambda_1:.8f}⁶ = {lambda_1_pow_6:.8f}")
    print(f"|λ₂|⁶ = {lambda_2:.8f}⁶ = {lambda_2_pow_6:.8f}")
    
    print("\nFinal Equation:")
    print(f"Fidelity = {lambda_2_pow_6:.8f} / ({lambda_1_pow_6:.8f} + {lambda_2_pow_6:.8f})")

    print(f"\nThe calculated fidelity is: {fidelity:.8f}")
    
    # Final answer in the specified format
    print(f"<<<{fidelity}>>>")

solve_quantum_fidelity()