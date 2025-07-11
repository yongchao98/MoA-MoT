from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.groups import PermutationGroup

def describe_group(name, G):
    """Prints basic properties of a given permutation group."""
    print(f"--- Group: {name} ---")
    print(f"Order: {G.order()}")
    print(f"Is Abelian: {G.is_abelian()}")
    print(f"Generators: {G.generators}")
    print("-" * (len(name) + 14))
    print()

def main():
    """
    This script demonstrates the construction of two non-isomorphic nonabelian
    filled groups of order 18, as examples of the general classification.
    A group of order 18 has the form 2*q^m with q=3, m=2.

    According to the theory, such groups are of the form Q x| C_2, where Q is a
    group of order 9 and the action of C_2 induces inversion on Q/Phi(Q).
    """

    # Example 1: The Dihedral Group D_18 = D(C_9)
    # This is a filled group where Q is the cyclic group of order 9.
    # It can be represented as a permutation group on 9 elements.
    # r is the rotation, s is the reflection.
    r_d18 = Permutation(0, 1, 2, 3, 4, 5, 6, 7, 8)
    s_d18 = Permutation(1, 8)(2, 7)(3, 6)(4, 5)
    D18 = PermutationGroup([r_d18, s_d18])
    describe_group("Dihedral Group D_18", D18)

    # Example 2: The Generalized Dihedral Group D(C_3 x C_3)
    # This is a filled group where Q is the elementary abelian group C_3 x C_3.
    # We can represent this group as a permutation group on the 9 elements
    # of Q = {(i,j) | i,j in {0,1,2}}.
    # Let's label the element (i,j) as 3*i + j.
    # (0,0)->0, (0,1)->1, (0,2)->2, (1,0)->3, ... , (2,2)->8
    
    # Generator 'a' corresponds to translation by (1,0) in Q.
    # (i,j) -> (i+1, j)
    a_perm = Permutation(0, 3, 6)(1, 4, 7)(2, 5, 8)
    
    # Generator 'b' corresponds to translation by (0,1) in Q.
    # (i,j) -> (i, j+1)
    b_perm = Permutation(0, 1, 2)(3, 4, 5)(6, 7, 8)
    
    # Generator 'x' corresponds to the inversion automorphism on Q.
    # x(i,j) = (-i, -j) mod 3.
    # x(0,0)=(0,0) -> 0->0
    # x(0,1)=(0,2) -> 1->2
    # x(0,2)=(0,1) -> 2->1
    # x(1,0)=(2,0) -> 3->6
    # x(1,1)=(2,2) -> 4->8
    # x(1,2)=(2,1) -> 5->7
    # x(2,0)=(1,0) -> 6->3
    # x(2,1)=(1,2) -> 7->5
    # x(2,2)=(1,1) -> 8->4
    x_perm = Permutation(1, 2)(3, 6)(4, 8)(5, 7)
    
    D_C3xC3 = PermutationGroup([a_perm, b_perm, x_perm])
    describe_group("Generalized Dihedral Group D(C_3 x C_3)", D_C3xC3)

if __name__ == '__main__':
    main()