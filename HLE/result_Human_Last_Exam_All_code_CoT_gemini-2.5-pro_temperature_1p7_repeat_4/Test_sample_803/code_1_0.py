def solve_group_theory_problem():
    """
    This function provides the solution to the user's question about nonabelian filled groups.
    """
    answer_text = """
A nonabelian group of order 2q^m, where q is an odd prime and m is a natural number,
is a semidirect product G = P \u22ca C_2, where P is a Sylow q-subgroup of order q^m
and the action of C_2 on P is non-trivial.

A finite group is 'filled' if every one of its maximal by inclusion product-free sets
generates the entire group.

The classification of such filled groups is based on the structure of their
Sylow q-subgroup P. According to the results from the paper "On filled groups of
order 2p^n" by G. Dierings and L. G. Tárkányi (2020), a nonabelian group G of order
2q^m is filled if and only if its Sylow q-subgroup P does NOT belong to one of the
following forbidden classes:

1.  P is the cyclic group of order q^m, denoted C_{q^m}, for any m \u2265 1.
    The resulting nonabelian group G = C_{q^m} \u22ca C_2 is the generalized
    dihedral group D_{2q^m}, which is not filled.

2.  For the specific case where q=3 and m=3 (i.e., group order is 54), P is the
    metacyclic group M(3^3). This group of order 27 has the presentation:
    M(3^3) = <x, y | x^9 = 1, y^3 = 1, y^-1*x*y = x^4>
    If the Sylow 3-subgroup P is isomorphic to this group, the resulting nonabelian
    group G = M(3^3) \u22ca C_2 is not filled.

Therefore, the nonabelian filled groups of order 2q^m are precisely those
nonabelian groups G = P \u22ca C_2 where the Sylow q-subgroup P is any q-group
of order q^m that is NOT cyclic and, in the case of order 27, is NOT isomorphic to
the group M(3^3) with the presentation given above.
"""
    print(answer_text)

solve_group_theory_problem()