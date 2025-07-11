def solve_homology_cobordism_problem():
    """
    Determines and explains the number of homology cobordism group elements
    representable by integral surgery on knots with at most four crossings.
    """
    print("This script calculates the number of unique homology cobordism group elements from surgeries on knots with at most 4 crossings.")
    print("-" * 70)

    # Step 1 & 2: Identify knots and valid surgeries
    print("Step 1: Identifying relevant knots and surgery types.")
    print("The knots with at most four crossings are the unknot (0_1), the trefoil knot (3_1), and the figure-eight knot (4_1).")
    print("To obtain a homology 3-sphere (an element of the group), surgery on a knot must use a coefficient of +1 or -1.")
    print("-" * 70)

    # Step 3: Analyze each knot
    print("Step 2: Analyzing the surgery results for each knot.")

    # Unknot (0_1)
    print("\n- For the unknot (0_1):")
    print("  Surgery with coefficients +1 or -1 on the unknot produces the standard 3-sphere, S³.")
    print("  The 3-sphere S³ is the identity element in the homology cobordism group.")
    identity_elements_found = 1
    print(f"  Unique elements found so far: {identity_elements_found}")

    # Figure-eight knot (4_1)
    print("\n- For the figure-eight knot (4_1):")
    print("  The figure-eight knot is a 'slice' knot.")
    print("  A theorem states that +1 or -1 surgery on any slice knot results in a homology sphere that is cobordant to S³.")
    print("  Therefore, surgery on the figure-eight knot also yields the identity element.")
    print(f"  This does not introduce a new element. Unique elements found so far: {identity_elements_found}")

    # Trefoil knot (3_1)
    print("\n- For the trefoil knot (3_1):")
    print("  The trefoil knot is not slice. Surgeries on it are known to produce non-trivial results.")
    print("  Specifically, ±1 surgery on the trefoil (and its mirror image) yields the Poincaré homology sphere.")
    print("  The Poincaré sphere represents a non-trivial element of order 2 in the homology cobordism group.")
    nontrivial_elements_found = 1
    print(f"  This introduces a new, non-trivial element.")
    
    total_distinct_elements = identity_elements_found + nontrivial_elements_found
    print("-" * 70)

    # Step 4: Count the unique elements
    print("Step 3: Final Count.")
    print("The distinct elements found are:")
    print("1. The identity element (from the unknot and figure-eight knot).")
    print("2. The Poincaré sphere class (from the trefoil knot).")
    
    print("\nFinal Calculation:")
    print(f"Number of identity elements: {identity_elements_found}")
    print(f"Number of non-trivial elements: {nontrivial_elements_found}")
    print(f"Total number of distinct elements = {identity_elements_found} + {nontrivial_elements_found} = {total_distinct_elements}")

solve_homology_cobordism_problem()