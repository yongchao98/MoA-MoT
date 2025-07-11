import numpy as np

# Step 1: Define the K-matrix for the nu=2 Bosonic Integer Quantum Hall state
K_b = np.array([[0, 1],
                [1, 0]])

print("The K-matrix for the final state of bosons is K_b:")
print(K_b)
print("-" * 30)

# Step 2: Define the topological correction term Sigma from the gapped fermions
# We assume they form a nu=1 IQH state for each species.
Sigma = np.array([[1, 0],
                  [0, 1]])

print("The correction matrix from the gapped fermion sector is Sigma:")
print(Sigma)
print("-" * 30)

# Step 3: Calculate the K-matrix for the Composite Fermions (K_cf)
# The equation is K_cf = 4 * K_b + Sigma
K_cf_term1 = 4 * K_b
K_cf = K_cf_term1 + Sigma

print("Calculating the Composite Fermion K-matrix using K_cf = 4*K_b + Sigma:")
print(f"4 * K_b = \n{K_cf_term1}")
print(f"K_cf = \n{K_cf_term1} + \n{Sigma} = \n{K_cf}")
print("-" * 30)


# Step 4: Calculate the K-matrix for the original fermions (K_f)
# The relation is K_f = K_cf - 2*I
I = np.identity(2, dtype=int)
K_f_term2 = 2 * I
K_f = K_cf - K_f_term2

print("Calculating the final K-matrix for the original fermions using K_f = K_cf - 2*I:")
print(f"K_f = \n{K_cf} - \n{K_f_term2}")

# Final Answer Output
final_k_matrix = K_f
print("\nThe final K-matrix of the resulting fractional state is:")
print(final_k_matrix)
print("\nWhich is:")
print(f"[[{final_k_matrix[0,0]}, {final_k_matrix[0,1]}],")
print(f" [{final_k_matrix[1,0]}, {final_k_matrix[1,1]}]]")
