import numpy as np

def solve_fidelity_puzzle():
    """
    Solves the quantum mechanics puzzle based on the hypothetical universe's rules.
    """
    # Define the observable operator O
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 1: Find the eigenvalues of the operator O.
    # For a Hermitian matrix, eigvalsh returns eigenvalues in ascending order.
    eigenvalues = np.linalg.eigvalsh(O)
    
    # Step 2: Identify the second-largest (smallest) and largest eigenvalues.
    # eigvalsh sorts them, so lambda_2 is the first element and lambda_1 is the second.
    lambda_2 = eigenvalues[0]  # Second-largest eigenvalue
    lambda_1 = eigenvalues[1]  # Largest eigenvalue

    # Step 3: Calculate the terms for the fidelity formula F = λ₂⁶ / (λ₁⁶ + λ₂⁶).
    lambda_2_pow_6 = lambda_2**6
    lambda_1_pow_6 = lambda_1**6
    denominator = lambda_1_pow_6 + lambda_2_pow_6

    # Step 4: Calculate the final fidelity.
    fidelity = lambda_2_pow_6 / denominator

    # Output the explanation and the result as requested.
    print("The observable operator is:")
    print(O)
    print("\nStep 1: Find the eigenvalues of the operator.")
    print(f"The eigenvalues are λ₁ = {lambda_1:.4f} (largest) and λ₂ = {lambda_2:.4f} (second-largest).")
    
    print("\nStep 2: Determine the fidelity formula based on the rules.")
    print("The resulting state |ψ'⟩ is proportional to λ₁³|v₁⟩ + λ₂³|v₂⟩.")
    print("The target state is |v₂⟩.")
    print("The fidelity F is given by the formula: F = |⟨v₂|ψ'⟩|² = λ₂⁶ / (λ₁⁶ + λ₂⁶)")

    print("\nStep 3: Calculate the numerical values for the equation.")
    print(f"λ₁⁶ = ({lambda_1:.4f})⁶ ≈ {lambda_1_pow_6:.4f}")
    print(f"λ₂⁶ = ({lambda_2:.4f})⁶ ≈ {lambda_2_pow_6:.4f}")
    print(f"F = {lambda_2_pow_6:.4f} / ({lambda_1_pow_6:.4f} + {lambda_2_pow_6:.4f})")

    print("\nFinal Answer:")
    print(f"The fidelity is {fidelity:.6f}")
    
    # Return the raw numeric answer for the final capture
    return fidelity

# Execute the function
final_answer = solve_fidelity_puzzle()
# The <<<...>>> format is for the final answer extraction.
# The numeric value is printed above for the user to see clearly.
print(f"\n<<<{final_answer}>>>")
