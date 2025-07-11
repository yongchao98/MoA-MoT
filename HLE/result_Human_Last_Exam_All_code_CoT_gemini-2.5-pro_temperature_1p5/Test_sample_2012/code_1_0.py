import numpy as np

# This script calculates the fidelity based on the described hypothetical physics.

# Define the observable operator matrix
O = np.array([[3, 1],
              [1, 2]])

# Step 1: Find the eigenvalues of the operator O.
eigenvalues = np.linalg.eigvals(O)

# Step 2: Sort the eigenvalues to identify the largest and second-largest.
# np.linalg.eigvals does not guarantee order, so we sort.
eigenvalues.sort()  # Sorts in ascending order
# For a 2D system, the second-largest is the smallest, and the largest is the largest.
lambda2 = eigenvalues[0]
lambda1 = eigenvalues[1]

# Step 3: Calculate the fidelity.
# According to Rule 2, the post-measurement state |ψ_final⟩ is a superposition
# with coefficients proportional to the cube of the eigenvalues.
# |ψ_final⟩ is proportional to (λ₁³|v₁⟩ + λ₂³|v₂⟩).
# The fidelity F with the eigenstate |v₂⟩ (associated with the second-largest eigenvalue λ₂) is:
# F = |⟨v₂ | ψ_final⟩|²
# Due to the orthonormality of eigenvectors, this simplifies to:
# F = (λ₂⁶) / (λ₁⁶ + λ₂⁶)

# Calculate the terms of the equation
lambda1_pow6 = lambda1**6
lambda2_pow6 = lambda2**6
denominator = lambda1_pow6 + lambda2_pow6
fidelity = lambda2_pow6 / denominator

# Step 4: Print the results, showing each number in the final equation.
print(f"The eigenvalues of the operator are λ₁ = {lambda1:.6f} (largest) and λ₂ = {lambda2:.6f} (second-largest).")
print("\nThe formula for the fidelity (F) of the post-measurement state with the second eigenstate is F = λ₂⁶ / (λ₁⁶ + λ₂⁶).")
print("\nCalculating the values for the formula:")
print(f"λ₁⁶ = ({lambda1:.6f})⁶ = {lambda1_pow6:.6f}")
print(f"λ₂⁶ = ({lambda2:.6f})⁶ = {lambda2_pow6:.6f}")

print("\nFinal Fidelity Equation:")
print(f"F = {lambda2_pow6:.6f} / ({lambda1_pow6:.6f} + {lambda2_pow6:.6f})")
print(f"F = {lambda2_pow6:.6f} / {denominator:.6f}")
print(f"\nThe calculated fidelity is: {fidelity:.8f}")
