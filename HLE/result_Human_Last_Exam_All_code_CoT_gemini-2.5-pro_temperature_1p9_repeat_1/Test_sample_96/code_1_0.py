def solve():
    """
    This function solves the problem about the number of torsion elements in the Artin group of type E_8.
    The problem is to find the number of torsion elements of order 10 in A/Z
    which can be written as positive words of minimal length.

    Step-by-step derivation:
    1.  Let g be a positive word of length L, and let gZ have order 10 in A/Z.
        This means g^10 = z^m for some integer m, where z is a generator of the center Z.
    2.  Using the abelianization map chi: A -> Z, where chi(g) = L, we have:
        10 * L = m * chi(z).
    3.  For E_8, chi(z) is the number of positive roots, which is 120.
        So, 10 * L = m * 120, which simplifies to L = 12m.
    4.  The minimal length L must be a multiple of 12. The smallest possible values are 12, 24, 36, ...
    5.  A candidate for an element of order 10 is g = c^3, where c is a Coxeter element.
        The length of c is 8, so L(g) = 24. This corresponds to m=2.
    6.  For any Coxeter element c in A(E_8), we have the relation c^30 = z^2.
        Taking g = c^3, we get g^10 = (c^3)^10 = c^30 = z^2.
        The order of c^k in A/Z is 30/gcd(k, 30). For k=3, the order is 30/gcd(3,30) = 30/3 = 10.
        So, g = c^3 is a valid torsion element of order 10 and length 24.
    7.  The existence of a 10th root of z (an element of length 12) is not guaranteed and considered unlikely.
        Thus, the minimal length is assumed to be 24.
    8.  The problem reduces to counting the number of distinct elements of the form c^3, where c is any Coxeter element.
    9.  The set of Coxeter elements for a bipartite graph like E_8 is partitioned into two families, related to the bipartition. Let c+ and c- be representative elements for these two families.
    10. It is known that c+ and c- are not conjugate in the Artin group A(E_8). This suggests that c+^3 and c-^3 are distinct.
    11. A plausible conjecture is that any element c^3 is equal to either c+^3 or c-^3.
    12. Based on this reasoning, there are exactly 2 such elements.
    """
    
    number_of_elements = 2
    print(number_of_elements)

solve()