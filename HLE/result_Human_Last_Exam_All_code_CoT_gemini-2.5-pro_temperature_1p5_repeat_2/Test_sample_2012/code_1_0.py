import numpy as np

# Define the observable operator matrix
O = np.array([
    [3, 1],
    [1, 2]
])

# Step 1: Find the eigenvalues of the operator O.
# np.linalg.eigh returns eigenvalues in ascending order for Hermitian matrices.
eigenvalues, eigenvectors = np.linalg.eigh(O)

# Step 2: Identify the second-largest (smaller) and largest eigenvalues.
# For a 2D system, the second-largest is simply the smallest.
lambda_S = eigenvalues[0] # Second-largest (smaller) eigenvalue
lambda_L = eigenvalues[1] # Largest eigenvalue

print(f"The observable operator is:\n{O}\n")
print(f"The largest eigenvalue (lambda_L) is: {lambda_L}")
print(f"The second-largest eigenvalue (lambda_S) is: {lambda_S}\n")

# Step 3: Calculate the terms for the fidelity formula F = (lambda_S^6) / (lambda_L^6 + lambda_S^6)
lambda_S_pow6 = lambda_S**6
lambda_L_pow6 = lambda_L**6
sum_of_powers = lambda_L_pow6 + lambda_S_pow6

# Step 4: Calculate the final fidelity
fidelity = lambda_S_pow6 / sum_of_powers

# Step 5: Print the final equation with all its components
print("The fidelity F is calculated as lambda_S^6 / (lambda_L^6 + lambda_S^6).")
print("\n--- Calculation Breakdown ---")
print(f"Value of lambda_S^6: {lambda_S_pow6}")
print(f"Value of lambda_L^6: {lambda_L_pow6}")
print(f"Value of (lambda_L^6 + lambda_S^6): {sum_of_powers}")
print("\n--- Final Equation ---")
print(f"Fidelity = {lambda_S_pow6} / {sum_of_powers}")
print(f"Fidelity = {fidelity}")
print("--------------------------\n")
print(f"The final fidelity is {fidelity}.")
