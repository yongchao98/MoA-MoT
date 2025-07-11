import numpy as np

# Step 1: Define the quantum system
# Observable operator O and initial state psi
O = np.array([[3, 1],
              [1, 2]])
psi = np.array([np.sqrt(3)/2, 1/2])

# Step 2: Find eigenvalues and eigenvectors of the observable O
# We use np.linalg.eigh for a Hermitian matrix.
eigenvalues, eigenvectors = np.linalg.eigh(O)

# Sort eigenvalues and corresponding eigenvectors in descending order
idx = eigenvalues.argsort()[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# The largest eigenvalue is lambda1, the second-largest is lambda2
lambda1 = eigenvalues[0]
lambda2 = eigenvalues[1]

# The corresponding eigenvectors are v1 and v2
v1 = eigenvectors[:, 0]
v2 = eigenvectors[:, 1] # This is the target state for the fidelity calculation

# Step 3: Express the initial state in the eigenbasis of O
# The coefficients are the inner products of the eigenvectors with the initial state.
# c1 = <v1|psi> and c2 = <v2|psi>
c1 = np.dot(v1, psi)
c2 = np.dot(v2, psi)

# Step 4: Calculate the components for the fidelity formula
# According to the problem's rules, the post-measurement state |psi_final> is
# proportional to O^3|psi>. The fidelity with the target state |v2> is given by:
# F = |c2 * lambda2^3| / sqrt(|c1|^2 * lambda1^6 + |c2|^2 * lambda2^6)

# Calculate the numerator of the fidelity formula
numerator = np.abs(c2 * (lambda2**3))

# Calculate the denominator of the fidelity formula
denominator = np.sqrt(np.abs(c1)**2 * (lambda1**6) + np.abs(c2)**2 * (lambda2**6))

# Step 5: Calculate the final fidelity
fidelity = numerator / denominator

# Print the numbers in the final equation as requested
print("The fidelity is calculated as the ratio of two numbers.")
print(f"Numerator |c2 * lambda2^3|: {numerator}")
print(f"Denominator sqrt(|c1|^2*lambda1^6 + |c2|^2*lambda2^6): {denominator}")
print(f"Fidelity = {numerator} / {denominator}")
print(f"The final fidelity is: {fidelity}")

# Output the final answer in the specified format
print(f"<<<{fidelity}>>>")