import numpy as np

def solve_cohomology_dim():
    """
    Calculates the dimension of the second cohomology group H^2(G,M).
    The problem is reduced to calculating the dimension of H^2(C_n, M),
    which is dim(ker(T-I)) - dim(Im(N)).
    """
    # n is the dimension of the vector space M, which is 128.
    # This also corresponds to the order of the cyclic group action.
    n = 128

    print(f"The problem reduces to computing the dimension of H^2(C_{n}, M), where n={n}.")
    print("This dimension is given by dim(M^C_n) - dim(N(M)).")
    print("-" * 30)

    # Construct the permutation matrix T for a cyclic permutation on a basis of size n.
    # T sends basis vector e_i to e_{i+1} (with e_n going to e_0).
    # Using programming indices 0 to n-1: e_i -> e_{i+1}, e_{n-1} -> e_0.
    T = np.zeros((n, n), dtype=np.float64)
    for i in range(n - 1):
        T[i + 1, i] = 1
    T[0, n - 1] = 1

    # Identity matrix of size n.
    I = np.identity(n, dtype=np.float64)

    # 1. Calculate dim(M^C_n) = dim(ker(T-I)).
    # This is given by n - rank(T-I).
    A = T - I
    rank_A = np.linalg.matrix_rank(A)
    dim_ker_T_minus_I = n - rank_A
    
    print(f"First, we compute the dimension of the space of invariants, M^C_n = ker(T-I).")
    print(f"dim(ker(T-I)) = n - rank(T-I)")
    print(f"The dimension n = {n}.")
    print(f"The rank of (T-I) is {int(rank_A)}.")
    print(f"So, dim(M^C_n) = {n} - {int(rank_A)} = {int(dim_ker_T_minus_I)}.")
    print("-" * 30)

    # 2. Calculate dim(N(M)), which is the rank of the norm operator N.
    # N is the sum of all powers of T from 0 to n-1.
    # The columns of N are all identical, being the sum of all basis vectors.
    N = np.ones((n, n), dtype=np.float64)
    rank_N = np.linalg.matrix_rank(N)
    dim_Im_N = rank_N

    print(f"Next, we compute the dimension of the image of the norm map, N(M).")
    print(f"dim(N(M)) = rank(N)")
    print(f"The rank of the norm matrix N is {int(rank_N)}.")
    print(f"So, dim(N(M)) = {int(dim_Im_N)}.")
    print("-" * 30)

    # 3. The dimension of the cohomology group H^2(G,M).
    dim_H2 = dim_ker_T_minus_I - dim_Im_N
    
    print("The final dimension of the cohomology group H^2(G,M) is:")
    print(f"dim(H^2(G,M)) = dim(M^C_n) - dim(N(M)) = {int(dim_ker_T_minus_I)} - {int(dim_Im_N)} = {int(dim_H2)}")

solve_cohomology_dim()