import numpy as np

# Define the observable operator and the initial state
O = np.array([[3, 1],
              [1, 2]])

psi = np.array([np.sqrt(3)/2, 1/2])

# Step 1: Find the eigenvalues and eigenvectors of the observable O
# np.linalg.eigh sorts eigenvalues in ascending order.
# evals[0] is the smaller (second-largest) eigenvalue, lambda_2
# evals[1] is the larger (largest) eigenvalue, lambda_1
evals, evecs = np.linalg.eigh(O)

lambda_2 = evals[0]
lambda_1 = evals[1]

# v2 is the eigenvector for lambda_2, v1 is the eigenvector for lambda_1
v2 = evecs[:, 0]
v1 = evecs[:, 1]

# Step 2: Decompose the initial state into the eigenbasis of O
# We calculate the coefficients c1 and c2. np.vdot is used for complex-safe inner product.
c1 = np.vdot(v1, psi)
c2 = np.vdot(v2, psi)

# Get the squared magnitudes of the coefficients
c1_sq = np.abs(c1)**2
c2_sq = np.abs(c2)**2

# Step 3: Calculate the 6th power of the eigenvalues
lambda1_6 = lambda_1**6
lambda2_6 = lambda_2**6

# Step 4: Calculate the fidelity using the derived formula
# F = |c2|^2 * lambda_2^6 / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)
numerator = c2_sq * lambda2_6
denominator = c1_sq * lambda1_6 + c2_sq * lambda2_6
fidelity = numerator / denominator

# Print the components of the calculation as requested
print("The fidelity F is calculated as: F = (|c2|^2 * lambda_2^6) / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)\n")
print(f"Eigenvalue lambda_1 (largest): {lambda_1:.4f}")
print(f"Eigenvalue lambda_2 (second-largest): {lambda_2:.4f}\n")
print(f"|c1|^2 (projection on v1): {c1_sq:.4f}")
print(f"|c2|^2 (projection on v2): {c2_sq:.4f}\n")
print(f"lambda_1^6: {lambda1_6:.4f}")
print(f"lambda_2^6: {lambda2_6:.4f}\n")

print("Substituting these values into the formula:")
print(f"F = ({c2_sq:.4f} * {lambda2_6:.4f}) / (({c1_sq:.4f} * {lambda1_6:.4f}) + ({c2_sq:.4f} * {lambda2_6:.4f}))")
print(f"F = {numerator:.4f} / ({c1_sq * lambda1_6:.4f} + {c2_sq * lambda2_6:.4f})")
print(f"F = {numerator:.4f} / {denominator:.4f}")
print(f"\nThe final fidelity is: {fidelity:.8f}")

<<<0.00230283>>>