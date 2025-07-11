import sys

def solve():
    """
    This script calculates the total number of smooth coverings of PSL(2,p)
    for a prime p > 5.

    The problem asks for the total number of quasi-simple groups G such that
    G/Z(G) is isomorphic to the simple group S = PSL(2,p). The number of
    such non-isomorphic groups is given by the number of subgroups of the
    Schur multiplier of S.

    For S = PSL(2, p) with p > 5 prime, the Schur multiplier M(S) is the
    cyclic group of order 2, Z_2.

    A cyclic group of order 2 has exactly two subgroups:
    1. The trivial subgroup {e}. This corresponds to the universal covering group, SL(2,p).
    2. The group Z_2 itself. This corresponds to the trivial covering group, PSL(2,p).

    Therefore, there is 1 non-trivial covering and 1 trivial covering.
    The total number is their sum.
    """
    
    # We can set a sample prime p > 5, although the result is independent of it.
    p = 7

    # Number of trivial coverings (e.g., PSL(2,p) itself)
    num_trivial_covering = 1

    # Number of non-trivial coverings (e.g., SL(2,p))
    num_nontrivial_covering = 1
    
    # The total number of coverings is the sum.
    total_coverings = num_trivial_covering + num_nontrivial_covering
    
    # As requested, printing the final equation with each number.
    print(f"{num_trivial_covering} + {num_nontrivial_covering} = {total_coverings}")

solve()