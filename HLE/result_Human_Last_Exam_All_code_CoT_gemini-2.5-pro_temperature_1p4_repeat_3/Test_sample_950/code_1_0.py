import math

def solve_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology ring
    of the space of 3-subspaces of R^5, which is the Grassmannian Gr_3(R^5).

    The rank of the torsion subgroup of the cohomology ring H*(X; Z) is the sum
    of the ranks of the torsion parts of each cohomology group H^i(X; Z).
    All torsion in the integral cohomology of a real Grassmannian is 2-torsion.
    The rank of Tors(H^i) is its dimension as a Z/2Z vector space.

    By the Universal Coefficient Theorem, Tors(H^i(X; Z)) is isomorphic to Tors(H_{i-1}(X; Z)).
    So we need to sum the ranks of the torsion parts of the homology groups H_j(X; Z).

    The integral homology groups of Gr_3(R^5) are known:
    H_0 = Z                  -> Tors(H_0) = 0,             rank = 0
    H_1 = Z/2Z               -> Tors(H_1) = Z/2Z,          rank = 1
    H_2 = 0                  -> Tors(H_2) = 0,             rank = 0
    H_3 = Z + Z/2Z           -> Tors(H_3) = Z/2Z,          rank = 1
    H_4 = Z/2Z               -> Tors(H_4) = Z/2Z,          rank = 1
    H_5 = Z + (Z/2Z)^2       -> Tors(H_5) = (Z/2Z)^2,      rank = 2
    The manifold has dimension 6, so we sum the ranks for j=0 to 5.
    """

    # Ranks of the torsion subgroups of the homology groups H_j for j = 0 to 5.
    homology_torsion_ranks = [0, 1, 0, 1, 1, 2]

    # The total rank is the sum of these individual ranks.
    total_rank = sum(homology_torsion_ranks)

    print("The rank of the torsion subgroup of H^*(Gr_3(R^5); Z) is the sum of the ranks of Tors(H_j(Gr_3(R^5); Z)) for j=0..5.")
    print("The individual ranks are based on the known homology groups.")
    print("\nCalculation:")
    
    # Building the equation string with each number
    equation_str = " + ".join(map(str, homology_torsion_ranks))
    
    print(f"{equation_str} = {total_rank}")

    print(f"\nThe rank of the torsion subgroup of the integral cohomology ring is {total_rank}.")

solve_torsion_rank()