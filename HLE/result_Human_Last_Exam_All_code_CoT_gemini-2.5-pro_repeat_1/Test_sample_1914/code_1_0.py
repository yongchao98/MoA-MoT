def count_categories():
    """
    Calculates the number of categories with 2 objects and 4 morphisms, up to isomorphism.

    The calculation is based on a case-by-case analysis of how the 4 morphisms
    are distributed among the four hom-sets of a 2-object category {A, B}.
    Two morphisms are always the identities id_A and id_B. We analyze the placement
    of the remaining two non-identity morphisms.
    """

    print("Calculating the number of categories with 2 objects and 4 morphisms.")
    print("This involves summing the number of non-isomorphic structures for each possible distribution of morphisms.\n")

    # Case 1: The two non-identity morphisms are both endomorphisms on object A.
    # The hom-set sizes are (3, 0, 0, 1). The set Hom(A, A) must form a monoid of order 3.
    # There are 5 known non-isomorphic monoids of order 3.
    # The case with sizes (1, 0, 0, 3) is isomorphic by swapping objects A and B, so it is not counted separately.
    case1_monoids_of_order_3 = 5
    print(f"Case 1: Morphisms form a monoid of order 3 on one object -> {case1_monoids_of_order_3} categories")

    # Case 2: Each object has one non-identity endomorphism.
    # The hom-set sizes are (2, 0, 0, 2). This corresponds to a product of two monoids of order 2.
    # There are 2 non-isomorphic monoids of order 2. Combining them gives 3 unique pairs
    # (up to isomorphism): (Z2, Z2), (Z2, idempotent), (idempotent, idempotent).
    case2_endo_pairs = 3
    print(f"Case 2: One non-identity endomorphism on each object -> {case2_endo_pairs} categories")

    # Case 3: Both non-identity morphisms are parallel arrows from A to B.
    # The hom-set sizes are (1, 2, 0, 1). No non-identity compositions are possible.
    # This defines a single, unique category structure.
    # The case (1, 0, 2, 1) is isomorphic.
    case3_parallel_arrows = 1
    print(f"Case 3: Two parallel arrows from one object to another -> {case3_parallel_arrows} category")

    # Case 4: One endomorphism on A and one arrow from A to B.
    # The hom-set sizes are (2, 1, 0, 1). This forces a composition rule.
    # There are 2 possible structures based on the monoid structure of Hom(A, A).
    case4_endo_and_arrow = 2
    print(f"Case 4: One endomorphism and one outgoing arrow -> {case4_endo_and_arrow} categories")

    # Case 5: One endomorphism on A and one arrow from B to A.
    # The hom-set sizes are (2, 0, 1, 1). This is the "opposite" of Case 4.
    # These structures are not isomorphic to those in Case 4 and also result in 2 categories.
    case5_opposite_of_case4 = 2
    print(f"Case 5: One endomorphism and one incoming arrow -> {case5_opposite_of_case4} categories")

    # Case 6: One arrow from A to B, and one from B to A.
    # The hom-set sizes are (1, 1, 1, 1). The composition rules force the objects A and B
    # to be isomorphic. This defines a single, unique category structure.
    case6_isomorphism = 1
    print(f"Case 6: An isomorphism between the two objects -> {case6_isomorphism} category")

    # Summing the counts from all non-isomorphic cases.
    total_categories = (case1_monoids_of_order_3 +
                        case2_endo_pairs +
                        case3_parallel_arrows +
                        case4_endo_and_arrow +
                        case5_opposite_of_case4 +
                        case6_isomorphism)

    print("\nFinal equation:")
    print(f"{case1_monoids_of_order_3} + {case2_endo_pairs} + {case3_parallel_arrows} + {case4_endo_and_arrow} + {case5_opposite_of_case4} + {case6_isomorphism} = {total_categories}")
    print(f"\nThere are {total_categories} categories with 2 objects and 4 morphisms, up to isomorphism.")
    return total_categories

if __name__ == '__main__':
    count_categories()