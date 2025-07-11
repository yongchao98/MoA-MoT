import numpy as np

# Step 1: Define the operator and the initial state
O = np.array([[3, 1],
              [1, 2]])

# Initial state |psi> = [sqrt(3)/2, 1/2]
psi = np.array([np.sqrt(3)/2, 1/2])

# Step 2: Find eigenvalues and eigenvectors of the operator O
# np.linalg.eigh returns eigenvalues in ascending order and corresponding eigenvectors
eigenvalues, eigenvectors = np.linalg.eigh(O)

# The second-largest eigenvalue is the first element, and the largest is the second.
lambda_2 = eigenvalues[0]
lambda_1 = eigenvalues[1]

# The corresponding eigenvectors
v2 = eigenvectors[:, 0]
v1 = eigenvectors[:, 1]

# Step 3: Decompose the initial state into the eigenbasis
# Calculate coefficients c1 and c2 by projecting |psi> onto the eigenvectors
# c1 = <v1|psi>, c2 = <v2|psi>
# .conj() is used for correctness with complex vectors, though not needed here.
c1 = np.dot(v1.conj(), psi)
c2 = np.dot(v2.conj(), psi)

# Step 4 & 5: Calculate the fidelity
# According to the derived formula: F = |c2|^2 * lambda_2^6 / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)
# Let's calculate the numerator and denominator
numerator = (np.abs(c2)**2) * (lambda_2**6)
denominator = (np.abs(c1)**2) * (lambda_1**6) + (np.abs(c2)**2) * (lambda_2**6)

fidelity = numerator / denominator

# Print the result
print(f"The observable operator is:\n{O}\n")
print(f"The initial state |psi> is: {psi}\n")
print(f"The largest eigenvalue (lambda_1) is: {lambda_1}")
print(f"The corresponding eigenvector |v1> is: {v1}\n")
print(f"The second-largest eigenvalue (lambda_2) is: {lambda_2}")
print(f"The corresponding eigenvector |v2> is: {v2}\n")
print(f"The projection of |psi> onto |v1> (c1) is: {c1}")
print(f"The projection of |psi> onto |v2> (c2) is: {c2}\n")
print("The fidelity is calculated using the formula F = |c2|^2 * lambda_2^6 / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)")
print(f"F = ({np.abs(c2)**2:.6f} * {lambda_2**6:.6f}) / (({np.abs(c1)**2:.6f} * {lambda_1**6:.6f}) + ({np.abs(c2)**2:.6f} * {lambda_2**6:.6f}))")
print(f"F = {numerator:.6f} / ({denominator:.6f})")
print(f"The final fidelity is: {fidelity}")

<<<2.8131333734020084e-06>>>