import numpy as np

# This script calculates the fidelity of a quantum state in a hypothetical universe
# based on a set of modified physical laws.

# Define the observable operator and the initial state
O = np.array([
    [3, 1],
    [1, 2]
])
psi = np.array([np.sqrt(3)/2, 1/2])

print("--- Step 1: Find Eigenvalues and Eigenvectors of the Operator O ---")
# Find the eigenvalues and eigenvectors of the operator O
eigvals, eigvecs = np.linalg.eig(O)

# Sort the eigenvalues and corresponding eigenvectors in descending order
sorted_indices = np.argsort(eigvals)[::-1]
sorted_eigvals = eigvals[sorted_indices]
sorted_eigvecs = eigvecs[:, sorted_indices]

# Assign eigenvalues and eigenvectors to named variables
lambda_1 = sorted_eigvals[0]  # Largest eigenvalue
lambda_2 = sorted_eigvals[1]  # Second-largest eigenvalue
v_1 = sorted_eigvecs[:, 0]
v_2 = sorted_eigvecs[:, 1]  # Target eigenstate for fidelity calculation

print(f"Largest eigenvalue (lambda_1): {lambda_1}")
print(f"Second-largest eigenvalue (lambda_2): {lambda_2}\n")

print("--- Step 2: Express the initial state |psi> in the eigenbasis ---")
# Project the initial state onto the eigenvectors to find the coefficients c1 and c2
c1 = np.dot(v_1.conj().T, psi)
c2 = np.dot(v_2.conj().T, psi)

# For the fidelity formula, we need the squared magnitudes of the coefficients
c1_sq = np.abs(c1)**2
c2_sq = np.abs(c2)**2

print(f"Squared magnitude of coefficient c1, |<v1|psi>|^2: {c1_sq}")
print(f"Squared magnitude of coefficient c2, |<v2|psi>|^2: {c2_sq}\n")
# Sanity check: the sum of squared magnitudes should be 1 since |psi> is normalized.
# print(f"Check: |c1|^2 + |c2|^2 = {c1_sq + c2_sq}\n")


print("--- Step 3: Calculate the Fidelity ---")
# According to the problem, the fidelity F of the post-measurement state
# with respect to |v2> is given by the formula:
# F = (|c2|^2 * lambda_2^6) / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)

# Calculate the 6th power of the eigenvalues
lambda_1_pow6 = lambda_1**6
lambda_2_pow6 = lambda_2**6

# Calculate the numerator and denominator of the fidelity formula
numerator = c2_sq * lambda_2_pow6
denominator = (c1_sq * lambda_1_pow6) + (c2_sq * lambda_2_pow6)

# Calculate the final fidelity
fidelity = numerator / denominator

print("The final equation for fidelity is F = (A * B) / (C * D + A * B), where:")
print(f"  A = |c2|^2    = {c2_sq}")
print(f"  B = lambda_2^6 = {lambda_2_pow6}")
print(f"  C = |c1|^2    = {c1_sq}")
print(f"  D = lambda_1^6 = {lambda_1_pow6}")
print("\nSubstituting the numbers:")
print(f"Numerator = {c2_sq} * {lambda_2_pow6} = {numerator}")
print(f"Denominator = ({c1_sq} * {lambda_1_pow6}) + ({c2_sq} * {lambda_2_pow6}) = {denominator}")

print("\nFinal Fidelity F = Numerator / Denominator")
print(f"F = {numerator} / {denominator}")
print(f"F = {fidelity}")

# <<<Final Answer>>>
print(f"\n<<<{fidelity}>>>")
