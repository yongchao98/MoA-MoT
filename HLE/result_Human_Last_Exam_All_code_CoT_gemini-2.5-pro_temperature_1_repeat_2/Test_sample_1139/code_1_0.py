import math

def calculate_bosonic_variables(N):
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model for disordered systems of symmetry class D with N replicas.

    The formula for the number of variables is the dimension of the bosonic manifold:
    dim[SO(2N)/(SO(N)xSO(N))] + dim[Sp(2N,R)/(Sp(N,R)xSp(N,R))]
    """

    # Lie group dimension formulas
    # dim(SO(n)) = n(n-1)/2
    def dim_so(n):
        return n * (n - 1) // 2

    # dim(Sp(2n, R)) = n(2n+1)
    # The input to this function is the matrix size 'm' (e.g., m=2n)
    def dim_sp(m):
        if m % 2 != 0:
            raise ValueError("The dimension for Sp(m,R) must be even.")
        n = m // 2
        return n * (2 * n + 1)

    print(f"The calculation is for N = {N} replicas.")
    print("-" * 30)

    # --- SO part calculation ---
    dim_so_2N = dim_so(2 * N)
    dim_so_N = dim_so(N)
    dim_so_part = dim_so_2N - 2 * dim_so_N

    print("The number of variables from the SO sector is:")
    print(f"dim(SO({2*N})) - 2 * dim(SO({N}))")
    print(f"= {dim_so_2N} - 2 * {dim_so_N}")
    print(f"= {dim_so_2N} - {2 * dim_so_N}")
    print(f"= {dim_so_part}")
    print("-" * 30)

    # --- Sp part calculation ---
    dim_sp_2N = dim_sp(2 * N)
    dim_sp_N = dim_sp(N)
    dim_sp_part = dim_sp_2N - 2 * dim_sp_N

    print("The number of variables from the Sp sector is:")
    print(f"dim(Sp({2*N},R)) - 2 * dim(Sp({N},R))")
    print(f"= {dim_sp_2N} - 2 * {dim_sp_N}")
    print(f"= {dim_sp_2N} - {2 * dim_sp_N}")
    print(f"= {dim_sp_part}")
    print("-" * 30)

    # --- Total calculation ---
    total_variables = dim_so_part + dim_sp_part
    print("Total number of non-Grassman variables is the sum of the two parts:")
    print(f"Total = (dim SO part) + (dim Sp part)")
    print(f"Total = {dim_so_part} + {dim_sp_part} = {total_variables}")

if __name__ == '__main__':
    # Number of replicas for the given problem
    N_replicas = 2
    calculate_bosonic_variables(N_replicas)
