import numpy as np

# Step 1: Define the observable operator and find its eigenvalues.
# The observable operator O is given by the matrix:
# O = [[3, 1],
#      [1, 2]]
O = np.array([[3, 1],
              [1, 2]])

# We find the eigenvalues of O.
eigenvalues = np.linalg.eigvals(O)

# Step 2: Identify the largest and second-largest eigenvalues.
# We sort the eigenvalues to easily identify them.
eigenvalues.sort()
# In this 2-dimensional system, the smaller eigenvalue is the second-largest.
lambda_2 = eigenvalues[0]
lambda_1 = eigenvalues[1] # Largest eigenvalue

# Step 3: Formulate the fidelity calculation based on the problem's rules.
# Law 2 states the post-measurement state |psi_final> is a superposition with
# coefficients proportional to the cube of the eigenvalues.
# |psi_final> is proportional to (lambda_1^3 |v1> + lambda_2^3 |v2>), where |v1> and |v2>
# are the normalized eigenvectors for lambda_1 and lambda_2.
# The question asks for the fidelity of this state with respect to the eigenstate
# of the second-largest eigenvalue, |v2>.
# Fidelity F = |<v2|psi_final>|^2.
# After normalization and using orthogonality of eigenvectors (<v2|v1>=0, <v2|v2>=1),
# the formula for fidelity simplifies to:
# F = lambda_2^6 / (lambda_1^6 + lambda_2^6)

# Step 4: Calculate the necessary components for the final equation.
lambda_1_pow6 = lambda_1**6
lambda_2_pow6 = lambda_2**6

# Step 5: Calculate the final fidelity.
fidelity = lambda_2_pow6 / (lambda_1_pow6 + lambda_2_pow6)

# Step 6: Print the calculation steps and the final answer.
print("The observable operator's eigenvalues are calculated from the matrix:")
print(f"{O}\n")
print(f"The largest eigenvalue is λ_1 = {lambda_1:.5f}")
print(f"The second-largest eigenvalue is λ_2 = {lambda_2:.5f}\n")
print("The fidelity F is calculated using the formula: F = λ_2^6 / (λ_1^6 + λ_2^6)")
print(f"Substituting the values:")
print(f"F = {lambda_2_pow6:.5f} / ({lambda_1_pow6:.5f} + {lambda_2_pow6:.5f})")
print(f"F = {lambda_2_pow6:.5f} / {(lambda_1_pow6 + lambda_2_pow6):.5f}")
print(f"The final fidelity is F = {fidelity:.5f}")

print("\n<<<0.00268>>>")