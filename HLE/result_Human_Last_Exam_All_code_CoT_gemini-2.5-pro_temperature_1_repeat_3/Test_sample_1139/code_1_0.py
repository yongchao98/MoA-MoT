import math

def calculate_non_grassman_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model for disordered systems of symmetry class D with a given
    number of replicas.
    """
    # Number of replicas
    num_replicas = 2

    def dim_orthogonal_group(k):
        """Calculates the dimension of the orthogonal group O(k)."""
        if k < 1:
            return 0
        return k * (k - 1) // 2

    def dim_bosonic_supergroup(p, q):
        """
        Calculates the dimension of the bosonic part of the supergroup O(p|q).
        The bosonic subgroup of O(p|q) is O(p) x O(q).
        """
        return dim_orthogonal_group(p) + dim_orthogonal_group(q)

    # For symmetry class D, the target supermanifold is G/H.
    # G = O(2*N_R | 2*N_R)
    # H = O(N_R | N_R) x O(N_R | N_R)

    # Parameters for G = O(p_G | q_G)
    p_G = 2 * num_replicas
    q_G = 2 * num_replicas

    # Parameters for one component of H = O(p_H | q_H)
    p_H = num_replicas
    q_H = num_replicas

    # Calculate the bosonic dimension of the total group G
    dim_b_g = dim_bosonic_supergroup(p_G, q_G)

    # Calculate the bosonic dimension of the subgroup H
    # H has two identical components, so its dimension is twice the dimension of one component.
    dim_b_h = 2 * dim_bosonic_supergroup(p_H, q_H)

    # The number of non-Grassman variables is the difference in the bosonic dimensions.
    num_variables = dim_b_g - dim_b_h

    # Print the calculation steps
    print(f"The number of non-Grassman variables is the dimension of the bosonic part of the coset space G/H.")
    print(f"For symmetry class D with {num_replicas} replicas:")
    print(f"G = O({p_G}|{q_G}) and H = O({p_H}|{q_H}) x O({p_H}|{q_H}).")
    print("-" * 30)
    print(f"Dimension of the bosonic part of G, dim_B(G):")
    dim_o_p_g = dim_orthogonal_group(p_G)
    dim_o_q_g = dim_orthogonal_group(q_G)
    print(f"dim_B(O({p_G}|{q_G})) = dim(O({p_G})) + dim(O({q_G})) = {dim_o_p_g} + {dim_o_q_g} = {dim_b_g}")
    
    print(f"\nDimension of the bosonic part of H, dim_B(H):")
    dim_o_p_h = dim_orthogonal_group(p_H)
    dim_o_q_h = dim_orthogonal_group(q_H)
    dim_b_h_component = dim_bosonic_supergroup(p_H, q_H)
    print(f"dim_B(H) = 2 * (dim(O({p_H})) + dim(O({q_H}))) = 2 * ({dim_o_p_h} + {dim_o_q_h}) = {dim_b_h}")
    print("-" * 30)
    
    # Print the final equation with all numbers
    print("Final Calculation:")
    print(f"Number of variables = dim_B(G) - dim_B(H) = {dim_b_g} - {dim_b_h} = {num_variables}")

calculate_non_grassman_variables()