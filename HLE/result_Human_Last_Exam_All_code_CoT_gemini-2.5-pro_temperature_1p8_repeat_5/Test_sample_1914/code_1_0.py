def solve_category_count():
    """
    Calculates the number of non-isomorphic categories with 2 objects and 4 morphisms.
    """
    print("Step 1: Problem setup")
    print("Let the two objects be A and B.")
    print("Total morphisms = 4. Two must be identities: id_A and id_B.")
    print("This leaves 2 non-identity morphisms to distribute among Hom(A,A), Hom(B,B), Hom(A,B), Hom(B,A).\n")

    print("Step 2 & 3: Enumerate and group distributions of morphisms")
    # n_aa, n_bb, n_ab, n_ba
    # n_aa >= 1, n_bb >= 1. Total is 4.
    # Let x_aa = n_aa - 1, x_bb = n_bb - 1.
    # x_aa + x_ab + x_ba + x_bb = 2.
    # Solutions for (x_aa, x_bb, x_ab, x_ba) which sum to 2.
    solutions = []
    for i in range(3):
        for j in range(3-i):
            for k in range(3-i-j):
                l = 2 - i - j - k
                solutions.append((i,j,k,l))
    
    # Convert back to (n_aa, n_bb, n_ab, n_ba)
    distributions = [(s[0]+1, s[1]+1, s[2], s[3]) for s in solutions]
    
    # Group by isomorphism (swapping A and B)
    # The representative of each class will be the one that is lexicographically first.
    representatives = {}
    for d in distributions:
        # d = (n_aa, n_bb, n_ab, n_ba)
        # d_swapped = (n_bb, n_aa, n_ba, n_ab)
        d_swapped = (d[1], d[0], d[3], d[2])
        rep = min(d, d_swapped)
        if rep not in representatives:
            representatives[rep] = 0
            
    print("Found 6 unique distribution patterns (up to object isomorphism):")
    for r in sorted(representatives.keys()):
        print(f"  - {r}")
    print("\nStep 4: Analyze each case\n")

    total_categories = 0
    
    # Case 1: (3,1,0,0) - disjoint union of a 3-morphism monoid and a 1-morphism monoid.
    # Number of categories is the number of non-isomorphic monoids of order 3.
    num_monoids_3 = 7
    representatives[(3,1,0,0)] = num_monoids_3
    print(f"Case (3,1,0,0): Hom(A,A) has 3 morphisms, Hom(B,B) has 1. No interaction.")
    print("This corresponds to a disjoint union of two one-object categories (monoids).")
    print("The number of non-isomorphic monoids of order 3 is 7.")
    print(f"Number of categories = {num_monoids_3}\n")
    total_categories += num_monoids_3

    # Case 2: (1,1,2,0) - quiver with two parallel arrows.
    # All compositions are forced, leading to a unique structure.
    num_cat_1120 = 1
    representatives[(1,1,2,0)] = num_cat_1120
    print(f"Case (1,1,2,0): Hom(A,B) has 2 morphisms, others are just identities.")
    print("This is a simple quiver category. All compositions are fixed by identity laws.")
    print("There is no choice in the structure.")
    print(f"Number of categories = {num_cat_1120}\n")
    total_categories += num_cat_1120

    # Case 3: (2,2,0,0) - disjoint union of two 2-morphism monoids.
    # Monoids of order 2 are C2 (cyclic group) and S2 (semigroup {0,1}).
    # Pairs are (C2,C2), (S2,S2), (C2,S2). (S2,C2) is isomorphic to (C2,S2).
    num_cat_2200 = 3
    representatives[(2,2,0,0)] = num_cat_2200
    print(f"Case (2,2,0,0): Hom(A,A) has 2 morphisms, Hom(B,B) has 2.")
    print("A disjoint union of two monoids of order 2. There are 2 such monoids: C2 and S2.")
    print("Possible pairs (up to isomorphism): (C2, C2), (S2, S2), (C2, S2).")
    print(f"Number of categories = {num_cat_2200}\n")
    total_categories += num_cat_2200
    
    # Case 4: (2,1,1,0) - M=Hom(A,A) is a monoid of order 2, P=Hom(A,B) has 1 morphism g.
    # Action of M on P (composition g o f) is forced to be g.
    # The structure only depends on the choice for M (C2 or S2).
    num_cat_2110 = 2
    representatives[(2,1,1,0)] = num_cat_2110
    print(f"Case (2,1,1,0): Hom(A,A) has 2 morphs (f, id_A), Hom(A,B) has 1 (g).")
    print("Hom(A,A) can be monoid C2 or S2. The composition g o f must be g.")
    print("Both choices for Hom(A,A) lead to a valid, distinct category.")
    print(f"Number of categories = {num_cat_2110}\n")
    total_categories += num_cat_2110
    
    # Case 5: (2,1,0,1) - dual of the previous case.
    # Hom(A,A) is M_2 (2 choices), Hom(B,A) has 1 morphism h.
    # These are not isomorphic to the previous case's categories.
    num_cat_2101 = 2
    representatives[(2,1,0,1)] = num_cat_2101
    print(f"Case (2,1,0,1): Hom(A,A) has 2 morphs, Hom(B,A) has 1.")
    print("This is anti-isomorphic to the previous case, but not isomorphic.")
    print("Again, two choices for Hom(A,A) give two distinct categories.")
    print(f"Number of categories = {num_cat_2101}\n")
    total_categories += num_cat_2101
    
    # Case 6: (1,1,1,1) - f: A->B, g: B->A.
    # Composition must be g o f = id_A and f o g = id_B.
    # This defines the category of an isomorphism between A and B. Unique.
    num_cat_1111 = 1
    representatives[(1,1,1,1)] = num_cat_1111
    print(f"Case (1,1,1,1): Hom(A,B)={f}, Hom(B,A)={g}. Others are identities.")
    print("Compositions must be g o f = id_A and f o g = id_B.")
    print("This describes a single, rigid structure where A and B are isomorphic.")
    print(f"Number of categories = {num_cat_1111}\n")
    total_categories += num_cat_1111

    print("Step 5: Final Calculation")
    print("Total number of categories is the sum of counts from each case.")
    calculation_str = f"{num_monoids_3} + {num_cat_1120} + {num_cat_2200} + {num_cat_2110} + {num_cat_2101} + {num_cat_1111} = {total_categories}"
    print(calculation_str)
    
    return total_categories

final_answer = solve_category_count()
print(f"\nFinal Answer: The total number of categories with 2 objects and 4 morphisms, up to isomorphism, is {final_answer}.")
# This final print is for the user. The required output format is below.
# I'll just output the final number.
