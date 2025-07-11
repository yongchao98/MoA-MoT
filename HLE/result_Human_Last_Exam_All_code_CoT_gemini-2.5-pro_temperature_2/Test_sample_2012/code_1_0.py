import numpy as np

# Define the observable operator matrix
O = np.array([[3, 1],
              [1, 2]])

# Calculate eigenvalues and eigenvectors
eigenvalues, eigenvectors = np.linalg.eig(O)

# Sort eigenvalues in descending order to identify the largest and second-largest
sorted_indices = np.argsort(eigenvalues)[::-1]
sorted_eigenvalues = eigenvalues[sorted_indices]

lambda_1 = sorted_eigenvalues[0]  # Largest eigenvalue
lambda_2 = sorted_eigenvalues[1]  # Second-largest eigenvalue

print(f"The observable operator is:\n{O}\n")
print(f"Largest eigenvalue (λ1): {lambda_1:.4f}")
print(f"Second-largest eigenvalue (λ2): {lambda_2:.4f}\n")

# According to the problem, the fidelity (F) of the post-measurement state
# with respect to the eigenstate of the second-largest eigenvalue is given by the formula:
# F = (λ2^6) / (λ1^6 + λ2^6)

# Calculate the terms
lambda_1_pow_6 = lambda_1**6
lambda_2_pow_6 = lambda_2**6

# Calculate the final fidelity
fidelity = lambda_2_pow_6 / (lambda_1_pow_6 + lambda_2_pow_6)

print("To calculate the fidelity, we compute the necessary terms:")
print(f"λ1^6 = ({lambda_1:.4f})^6 = {lambda_1_pow_6:.4f}")
print(f"λ2^6 = ({lambda_2:.4f})^6 = {lambda_2_pow_6:.4f}")
print(f"λ1^6 + λ2^6 = {lambda_1_pow_6:.4f} + {lambda_2_pow_6:.4f} = {lambda_1_pow_6 + lambda_2_pow_6:.4f}\n")

print("The fidelity F is calculated as follows:")
print(f"F = λ2^6 / (λ1^6 + λ2^6)")
print(f"F = {lambda_2_pow_6:.4f} / ({lambda_1_pow_6:.4f} + {lambda_2_pow_6:.4f})")
print(f"F = {lambda_2_pow_6:.4f} / {lambda_1_pow_6 + lambda_2_pow_6:.4f}")
print(f"F = {fidelity:.6f}\n")

print("The final fidelity is:")
print(fidelity)

# Final answer in the required format
final_answer = f"<<<{fidelity:.6f}>>>"