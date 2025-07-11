import numpy as np

# Step 1: Define the observable matrix O and the initial state psi
O = np.array([[3, 1],
              [1, 2]])
psi = np.array([np.sqrt(3)/2, 1/2])

# Step 2: Find the eigenvalues and eigenvectors of the matrix O
eigenvalues, eigenvectors = np.linalg.eig(O)

# Step 3: Identify the largest and second-largest (i.e., smallest) eigenvalues
# and their corresponding eigenvectors. We sort them to ensure correct assignment.
sort_indices = np.argsort(eigenvalues)
lambda_small, lambda_large = eigenvalues[sort_indices]
# eigenvectors are columns, so we sort them and transpose for easier access
v_small, v_large = eigenvectors[:, sort_indices].T

# Step 4: Calculate the coefficients of the initial state in the eigenbasis.
# This is done by projecting the initial state onto each eigenvector.
c_small = np.dot(v_small.conj(), psi)
c_large = np.dot(v_large.conj(), psi)

# Step 5: Calculate the squared magnitudes of the coefficients and the 6th power of the eigenvalues
c_small_sq = np.abs(c_small)**2
c_large_sq = np.abs(c_large)**2

lambda_small_6 = lambda_small**6
lambda_large_6 = lambda_large**6

# Step 6: Calculate the components of the fidelity formula
# F = |c_small|^2 * lambda_small^6 / ( |c_large|^2 * lambda_large^6 + |c_small|^2 * lambda_small^6 )
numerator = c_small_sq * lambda_small_6
denominator_part1 = c_large_sq * lambda_large_6
denominator_part2 = c_small_sq * lambda_small_6 # This is the same as the numerator

# Step 7: Calculate the final fidelity
fidelity = numerator / (denominator_part1 + denominator_part2)

# Print the breakdown of the calculation and the final result
print("The fidelity F is calculated as: |c_small|^2 * lambda_small^6 / (|c_large|^2 * lambda_large^6 + |c_small|^2 * lambda_small^6)\n")
print("Component values:")
print(f"  lambda_small = {lambda_small:.5f}, so lambda_small^6 = {lambda_small_6:.5f}")
print(f"  lambda_large = {lambda_large:.5f}, so lambda_large^6 = {lambda_large_6:.5f}")
print(f"  |c_small|^2 = {c_small_sq:.5f}")
print(f"  |c_large|^2 = {c_large_sq:.5f}\n")

print("Final equation with numbers:")
print(f"F = {numerator:.5f} / ({denominator_part1:.5f} + {denominator_part2:.5f})")
print(f"F = {numerator:.5f} / ({denominator_part1 + denominator_part2:.5f})")
print(f"F = {fidelity:.8f}")

print(f"\n<<<The final answer is {fidelity:.8f}>>>")
<<<0.00282668>>>