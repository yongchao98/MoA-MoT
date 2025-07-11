import math

def calculate_bosonic_variables():
    """
    Calculates the number of non-Grassmannian variables for the SUSY sigma-model
    of class D with a given number of replicas.
    """
    # Number of replicas for this problem
    n_replicas = 2
    
    # The symmetry class is D. The bosonic sector of the sigma-model
    # is parametrized by the symmetric space G/H, where for n replicas,
    # G = O(2*n) and H = O(n) x O(n).
    
    print("This problem asks for the number of non-Grassmannian (bosonic) variables")
    print(f"for the supersymmetric sigma-model for class D with {n_replicas} replicas.")
    print("This is equivalent to calculating the dimension of the target manifold, which is the symmetric space O(2n)/(O(n) x O(n)).")
    print(f"For n = {n_replicas}, the manifold is O({2*n_replicas}) / (O({n_replicas}) x O({n_replicas})).\n")
    
    def dim_orthogonal_group(k):
        """Calculates the dimension of the orthogonal group O(k)."""
        return k * (k - 1) // 2

    # Parameters for the groups G and H
    k_G = 2 * n_replicas
    k_H = n_replicas

    # Calculate dimension of the numerator group G = O(k_G)
    dim_G = dim_orthogonal_group(k_G)
    print(f"The dimension of the numerator group G = O({k_G}) is:")
    print(f"dim(O({k_G})) = {k_G} * ({k_G} - 1) / 2 = {dim_G}\n")

    # Calculate dimension of the denominator group H = O(k_H) x O(k_H)
    dim_H_single = dim_orthogonal_group(k_H)
    dim_H = dim_H_single + dim_H_single
    print(f"The dimension of the denominator group H = O({k_H}) x O({k_H}) is:")
    print(f"dim(O({k_H})) = {k_H} * ({k_H} - 1) / 2 = {dim_H_single}")
    print(f"dim(H) = dim(O({k_H})) + dim(O({k_H})) = {dim_H_single} + {dim_H_single} = {dim_H}\n")

    # The number of variables is the dimension of the manifold G/H
    num_variables = dim_G - dim_H
    
    print("The total number of non-Grassmannian variables is dim(G) - dim(H).")
    print("The final equation is:")
    print(f"{dim_G} - ({dim_H_single} + {dim_H_single}) = {num_variables}")

calculate_bosonic_variables()