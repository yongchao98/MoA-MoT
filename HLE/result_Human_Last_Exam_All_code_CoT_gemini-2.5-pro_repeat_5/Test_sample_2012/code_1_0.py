import numpy as np

# Step 1: Define the quantum system
# The observable operator O
O = np.array([[3, 1],
              [1, 2]], dtype=float)

# The initial state |psi>
psi_initial = np.array([np.sqrt(3)/2, 1/2], dtype=float)

print("This script calculates the fidelity of a post-measurement state in a fictional universe.")
print("-" * 40)
print(f"The observable operator O is:\n{O}")
print(f"The initial state |psi> is: {psi_initial}")
print("-" * 40)

# Step 2: Find the eigenvalues and eigenvectors of the observable O.
# numpy.linalg.eigh is used for Hermitian matrices and returns eigenvalues in ascending order.
eigenvalues, eigenvectors = np.linalg.eigh(O)

# In a 2D system, the "second-largest" eigenvalue is the smaller one.
lambda_2 = eigenvalues[0]  # Second-largest (smaller) eigenvalue
v_2 = eigenvectors[:, 0]   # Corresponding eigenvector

# The largest eigenvalue is the second one in the sorted list.
lambda_1 = eigenvalues[1]  # Largest eigenvalue
v_1 = eigenvectors[:, 1]   # Corresponding eigenvector

print("Step 2: Find eigenvalues and eigenvectors of O.")
print(f"The largest eigenvalue is lambda_1 = {lambda_1:.6f}")
print(f"The second-largest eigenvalue is lambda_2 = {lambda_2:.6f}")
print(f"The eigenstate for the second-largest eigenvalue, |v_2>, is: {v_2}")
print("-" * 40)

# Step 3: Express the initial state in the eigenbasis of O.
# The initial state can be written as |psi> = c1*|v1> + c2*|v2>.
# Since the eigenvectors |v1> and |v2> are orthonormal, the coefficients are found by projection.
c_1 = np.vdot(v_1, psi_initial)
c_2 = np.vdot(v_2, psi_initial)

print("Step 3: Express the initial state in the eigenbasis.")
print(f"The coefficient for the largest eigenstate is c1 = <v_1|psi> = {c_1:.6f}")
print(f"The coefficient for the second-largest eigenstate is c2 = <v_2|psi> = {c_2:.6f}")
print("-" * 40)

# Step 4 & 5: Calculate the fidelity based on the laws of universe U.
# Rule 2 states the new coefficients are proportional to the cube of the eigenvalues.
# This means the post-measurement state |psi_post> is proportional to:
# c1*(lambda_1^3)*|v1> + c2*(lambda_2^3)*|v2>
# The fidelity (F) with the state |v_2> is |<v_2|psi_post>|^2.
# After normalization, this leads to the formula:
# F = (|c2|^2 * lambda_2^6) / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)

c1_sq = np.abs(c_1)**2
c2_sq = np.abs(c_2)**2
l1_pow6 = lambda_1**6
l2_pow6 = lambda_2**6

numerator = c2_sq * l2_pow6
denominator = (c1_sq * l1_pow6) + (c2_sq * l2_pow6)
fidelity = numerator / denominator

print("Step 4 & 5: Calculate the final fidelity.")
print("The fidelity F is calculated using the formula:")
print("F = (|c2|^2 * lambda_2^6) / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)")
print("\nPlugging in the calculated values:")
print(f"F = ({c2_sq:.6f} * {l2_pow6:.6f}) / (({c1_sq:.6f} * {l1_pow6:.6f}) + ({c2_sq:.6f} * {l2_pow6:.6f}))")
print(f"F = {numerator:.6f} / ({c1_sq * l1_pow6:.6f} + {c2_sq * l2_pow6:.6f})")
print(f"F = {numerator:.6f} / {denominator:.6f}")
print("\n" + "-"*40)
print(f"The final fidelity is: {fidelity}")
print("-" * 40)

<<<0.0000028198590315206253>>>