def dim_o(m):
    """Calculates the dimension of the orthogonal group O(m)."""
    return m * (m - 1) // 2

def dim_sp(m):
    """Calculates the dimension of the symplectic group Sp(2m)."""
    if m % 2 != 0:
        raise ValueError("Input to dim_sp must be an even number (2k).")
    k = m // 2
    return k * (m + 1)

def solve_class_d_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model with two replicas for symmetry class D.
    """
    # Number of replicas
    n = 2

    # Numerator G = OSp(2n|2n) = OSp(4|4)
    # Bosonic part is G_B = Sp(2n) x O(2n) = Sp(4) x O(4)
    dim_sp_2n = dim_sp(2 * n)
    dim_o_2n = dim_o(2 * n)
    dim_g_b = dim_sp_2n + dim_o_2n

    # Denominator K = O(n|n) x O(n|n) = O(2|2) x O(2|2)
    # Bosonic part is K_B = (O(n) x O(n)) x (O(n) x O(n))
    # which is 4 copies of O(n) = O(2)
    dim_o_n = dim_o(n)
    dim_k_b = 4 * dim_o_n

    # The number of variables is dim(G_B) - dim(K_B)
    result = dim_g_b - dim_k_b

    # Print the detailed calculation
    # The structure is (dim(Sp(2n)) + dim(O(2n))) - 4 * dim(O(n))
    print("The number of non-Grassman variables is the dimension of the bosonic manifold:")
    print(f"dim(Sp({2*n})) + dim(O({2*n})) - (4 * dim(O({n})))")
    print(f"= {dim_sp_2n} + {dim_o_2n} - (4 * {dim_o_n})")
    print(f"= {dim_g_b} - {dim_k_b}")
    print(f"= {result}")


solve_class_d_variables()
