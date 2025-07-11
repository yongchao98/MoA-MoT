import numpy as np

# Define the observable operator and the initial state
O = np.array([[3., 1.], [1., 2.]])
psi = np.array([np.sqrt(3)/2., 1./2.])

# Step 1: Find the eigenvalues and eigenvectors of the operator O.
# np.linalg.eigh sorts eigenvalues in ascending order.
eigenvalues, eigenvectors = np.linalg.eigh(O)

# The second-largest eigenvalue is the first one in the sorted list.
lambda_2 = eigenvalues[0]
v2 = eigenvectors[:, 0]

# The largest eigenvalue is the second one.
lambda_1 = eigenvalues[1]
v1 = eigenvectors[:, 1]

# Step 2: Decompose the initial state into the eigenbasis.
# c_i = <v_i|psi>
# np.vdot calculates the inner product (v_i^dagger * psi)
c1 = np.vdot(v1, psi)
c2 = np.vdot(v2, psi)

# Step 3 & 4: Calculate the components for the fidelity formula.
# F = |c2|^2 * lambda_2^6 / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)
c1_sq_mag = np.abs(c1)**2
c2_sq_mag = np.abs(c2)**2
lambda_1_pow6 = lambda_1**6
lambda_2_pow6 = lambda_2**6

numerator = c2_sq_mag * lambda_2_pow6
denominator = c1_sq_mag * lambda_1_pow6 + numerator

fidelity = numerator / denominator

# Print the components of the final equation and the result
print("Fidelity is calculated using the formula: F = |c2|^2 * λ2^6 / (|c1|^2 * λ1^6 + |c2|^2 * λ2^6)")
print("\nIntermediate values for the equation:")
print(f"|c1|^2 = {c1_sq_mag}")
print(f"|c2|^2 = {c2_sq_mag}")
print(f"λ1^6 (largest eigenvalue to the power of 6) = {lambda_1_pow6}")
print(f"λ2^6 (second-largest eigenvalue to the power of 6) = {lambda_2_pow6}")
print(f"\nFinal Equation:")
print(f"F = ({c2_sq_mag}) * ({lambda_2_pow6}) / (({c1_sq_mag}) * ({lambda_1_pow6}) + ({c2_sq_mag}) * ({lambda_2_pow6}))")
print(f"F = {numerator} / {denominator}")
print(f"\nThe calculated fidelity is: {fidelity}")

# Final answer in the required format
print(f"\n<<<{fidelity}>>>")