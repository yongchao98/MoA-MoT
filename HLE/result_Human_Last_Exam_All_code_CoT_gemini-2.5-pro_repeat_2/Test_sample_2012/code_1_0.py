import numpy as np

# Define the observable operator matrix
O = np.array([[3, 1],
              [1, 2]])

# Calculate the eigenvalues and eigenvectors
eigenvalues, eigenvectors = np.linalg.eig(O)

# The problem is defined by the eigenvalues. We need to identify the
# largest and second-largest. In a 2D system, this is just the larger
# and the smaller eigenvalue. We sort them to ensure consistency.
sorted_indices = np.argsort(eigenvalues)[::-1]  # Sort in descending order
lambda_1 = eigenvalues[sorted_indices[0]]  # Largest eigenvalue
lambda_2 = eigenvalues[sorted_indices[1]]  # Second-largest eigenvalue

# According to Rule 2, the post-measurement state is a superposition
# with coefficients proportional to the cube of the eigenvalues.
# Let the final unnormalized state be |psi_un> = (lambda_1^3)|phi_1> + (lambda_2^3)|phi_2>
# The fidelity with the eigenstate |phi_2> is given by:
# F = |<phi_2 | psi_final>| = | (lambda_2^3) / sqrt((lambda_1^3)^2 + (lambda_2^3)^2) |
# F = |lambda_2^3| / sqrt(lambda_1^6 + lambda_2^6)

# Since lambda_2 = (5 - sqrt(5))/2 is positive, the absolute value is not needed.
lambda_2_cubed = lambda_2**3
lambda_1_to_the_6 = lambda_1**6
lambda_2_to_the_6 = lambda_2**6
sum_of_powers = lambda_1_to_the_6 + lambda_2_to_the_6

# Calculate the final fidelity
fidelity = lambda_2_cubed / np.sqrt(sum_of_powers)

# Print the results as requested by the prompt
print("This script calculates the fidelity based on the rules of universe U.")
print("\nStep 1: Find the eigenvalues of the operator O.")
print(f"The eigenvalues are lambda_1 (largest) = {lambda_1:.6f} and lambda_2 (second-largest) = {lambda_2:.6f}")

print("\nStep 2: Determine the fidelity formula from the problem's rules.")
print("Fidelity = lambda_2^3 / sqrt(lambda_1^6 + lambda_2^6)")

print("\nStep 3: Calculate the components of the formula.")
print(f"Numerator: lambda_2^3 = {lambda_2_cubed:.6f}")
print(f"Term in denominator: lambda_1^6 = {lambda_1_to_the_6:.6f}")
print(f"Term in denominator: lambda_2^6 = {lambda_2_to_the_6:.6f}")
print(f"Sum in denominator: lambda_1^6 + lambda_2^6 = {sum_of_powers:.6f}")

print("\nStep 4: Compute the final fidelity.")
print(f"F = {lambda_2_cubed:.6f} / sqrt({sum_of_powers:.6f})")
print(f"The final fidelity is: {fidelity:.8f}")
