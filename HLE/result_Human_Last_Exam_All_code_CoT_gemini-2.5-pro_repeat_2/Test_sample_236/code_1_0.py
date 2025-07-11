def solve_homology_cobordism_elements():
    """
    Calculates the number of homology cobordism group elements from surgeries
    on knots with at most four crossings.
    """
    print("This program determines how many elements of the homology cobordism group can be")
    print("represented by an integral surgery on a knot with at most four crossings.")
    print("-" * 70)

    # We use a set to store the unique elements found.
    # We use symbolic names for the elements for clarity:
    # 'Id': The identity element (represented by S^3)
    # 'P': The element for the Poincare sphere Sigma(2,3,5)
    # 'P_inv': The inverse of P
    # 'S': The element for the manifold Sigma(2,3,7)
    # 'S_inv': The inverse of S
    unique_elements = set()
    knot_contributions = {}

    # Case 1: The Unknot (0_1)
    # It is achiral. ±1 surgery on the unknot yields S^3.
    print("1. Analyzing the Unknot (0_1):")
    unknot_elements = {'Id'}
    unique_elements.update(unknot_elements)
    knot_contributions['0_1'] = len(unknot_elements)
    print("   - Surgeries on the unknot produce the 3-sphere S^3.")
    print("   - S^3 represents the identity 'Id' in the homology cobordism group.")
    print(f"   - New unique elements found: {knot_contributions['0_1']} ({unknot_elements})")
    print("-" * 70)

    # Case 2: The Trefoil Knot (3_1)
    # It is chiral, so we consider both right-handed (+3_1) and left-handed (-3_1) versions.
    print("2. Analyzing the Trefoil Knot (3_1 and its mirror):")
    # Right-handed Trefoil
    # +1 surgery -> Poincare sphere Sigma(2,3,5), class 'P'
    # -1 surgery -> Sigma(2,3,7), class 'S'
    rh_trefoil_elements = {'P', 'S'}
    
    # Left-handed Trefoil
    # +1 surgery -> -Sigma(2,3,7), class 'S_inv'
    # -1 surgery -> -Sigma(2,3,5), class 'P_inv'
    lh_trefoil_elements = {'P_inv', 'S_inv'}
    
    trefoil_elements = rh_trefoil_elements.union(lh_trefoil_elements)
    newly_added = trefoil_elements - unique_elements
    unique_elements.update(trefoil_elements)
    knot_contributions['+/-3_1'] = len(newly_added)

    print("   - Right-handed trefoil surgeries yield the Poincare sphere ('P') and Sigma(2,3,7) ('S').")
    print("   - Left-handed trefoil surgeries yield their inverses ('P_inv' and 'S_inv').")
    print("   - These four elements ('P', 'P_inv', 'S', 'S_inv') are all distinct and non-trivial.")
    print(f"   - New unique elements found: {knot_contributions['+/-3_1']} ({trefoil_elements})")
    print("-" * 70)

    # Case 3: The Figure-Eight Knot (4_1)
    # It is amphichiral (non-chiral).
    # +1 surgery -> Dodecahedral Space, which is known to be homology cobordant to -Sigma(2,3,5) ('P_inv')
    # -1 surgery -> -Dodecahedral Space, which is homology cobordant to Sigma(2,3,5) ('P')
    print("3. Analyzing the Figure-Eight Knot (4_1):")
    fig8_elements = {'P_inv', 'P'}
    newly_added = fig8_elements - unique_elements
    unique_elements.update(fig8_elements)
    knot_contributions['4_1'] = len(newly_added)
    print("   - Surgeries on the figure-eight knot also produce elements represented by 'P' and 'P_inv'.")
    print("   - These elements have already been found from the trefoil knot.")
    print(f"   - New unique elements found: {knot_contributions['4_1']}")
    print("-" * 70)

    # Final Calculation
    print("4. Final Tally:")
    print("   The complete set of distinct elements is:", sorted(list(unique_elements)))
    
    total_count = len(unique_elements)
    
    c_01 = knot_contributions['0_1']
    c_31 = knot_contributions['+/-3_1']
    c_41 = knot_contributions['4_1']

    print("\nSumming the contributions of new unique elements from each knot type:")
    print(f"Unknot (0_1): {c_01}")
    print(f"Trefoil (+/-3_1): {c_31}")
    print(f"Figure-Eight (4_1): {c_41}")
    print(f"\nFinal Equation: {c_01} + {c_31} + {c_41} = {total_count}")
    print(f"There are {total_count} such elements.")

solve_homology_cobordism_elements()