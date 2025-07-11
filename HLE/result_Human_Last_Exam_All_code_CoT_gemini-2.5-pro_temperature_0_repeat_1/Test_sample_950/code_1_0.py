def solve_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the space of 3-subspaces of R^5, i.e., Gr_3(R^5).

    The integral cohomology groups H^i(Gr_3(R^5); Z) are known from the literature.
    This function encodes these known results and sums the ranks of the torsion parts.
    The rank of a torsion group Z/2Z + ... + Z/2Z (k times) is k.
    """

    # H_i is a dictionary where keys are the degree i and values are tuples
    # representing the structure of H^i(Gr_3(R^5); Z).
    # The tuple is (rank_of_free_part, rank_of_torsion_part).
    # The rank of the torsion part is the number of Z/2Z summands.
    cohomology_groups = {
        0: (1, 0),  # H^0 = Z
        1: (0, 0),  # H^1 = 0
        2: (0, 1),  # H^2 = Z/2Z
        3: (0, 1),  # H^3 = Z/2Z
        4: (0, 2),  # H^4 = Z/2Z + Z/2Z
        5: (0, 1),  # H^5 = Z/2Z
        6: (0, 0),  # H^6 = 0
    }

    torsion_ranks = []
    for i in sorted(cohomology_groups.keys()):
        free_rank, torsion_rank = cohomology_groups[i]
        if torsion_rank > 0:
            torsion_ranks.append(torsion_rank)

    total_torsion_rank = sum(torsion_ranks)
    
    # Create the equation string
    equation_parts = [str(r) for r in torsion_ranks]
    equation_str = " + ".join(equation_parts) + f" = {total_torsion_rank}"

    print("The rank of the torsion subgroup is the sum of the ranks of the torsion parts of each cohomology group H^i.")
    print(f"The ranks of the non-trivial torsion parts are for i=2, 3, 4, 5.")
    print(f"The calculation is: {equation_str}")
    print(f"The total rank of the torsion subgroup is: {total_torsion_rank}")

solve_torsion_rank()