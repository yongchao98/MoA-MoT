def solve_homology_cobordism_count():
    """
    Calculates the number of homology cobordism group elements
    representable by integral surgery on knots with at most four crossings.
    """
    print("This script calculates the number of distinct homology cobordism group elements")
    print("from integral surgeries on knots with at most four crossings.")
    print("-" * 70)

    # A set to store the unique homology cobordism elements found.
    # We use descriptive strings to represent the elements.
    # "Trivial" represents the identity element (from S^3).
    # Other strings represent known non-trivial homology spheres.
    elements_found = set()
    
    # --- Unknot (0_1) ---
    print("1. Analyzing the Unknot (0_1):")
    # +/-1 surgery on the unknot gives S^3, the trivial element.
    elements_from_unknot = {"Trivial"}
    new_elements_from_unknot = elements_from_unknot.difference(elements_found)
    elements_found.update(new_elements_from_unknot)
    print(f"   - Surgeries on the unknot contribute the trivial element.")
    print(f"   - New distinct elements found: {len(new_elements_from_unknot)}")
    print("-" * 70)

    # --- Trefoil Knot (3_1) ---
    print("2. Analyzing the Trefoil Knot (3_1):")
    print("   - The trefoil is chiral, so we consider both right-handed (RHT) and left-handed (LHT) versions.")
    # These surgeries yield four distinct non-trivial homology spheres.
    # RHT: +1 -> Sigma(2,3,7), -1 -> Poincare Sphere Sigma(2,3,5)
    # LHT: +1 -> -Sigma(2,3,5), -1 -> -Sigma(2,3,7)
    elements_from_trefoil = {
        "Poincare_Sphere", 
        "Mirror_Poincare_Sphere", 
        "Sigma(2,3,7)", 
        "Mirror_Sigma(2,3,7)"
    }
    new_elements_from_trefoil = elements_from_trefoil.difference(elements_found)
    elements_found.update(new_elements_from_trefoil)
    print(f"   - Surgeries on the trefoil and its mirror contribute 4 distinct non-trivial elements.")
    print(f"   - New distinct elements found: {len(new_elements_from_trefoil)}")
    print("-" * 70)

    # --- Figure-eight Knot (4_1) ---
    print("3. Analyzing the Figure-eight Knot (4_1):")
    # The figure-eight knot is slice. +/-1 surgery on a slice knot results
    # in a homology sphere that is homology cobordant to S^3 (the trivial element).
    elements_from_fig8 = {"Trivial"}
    new_elements_from_fig8 = elements_from_fig8.difference(elements_found)
    elements_found.update(new_elements_from_fig8)
    print(f"   - The figure-eight knot is slice, so its surgeries only contribute the trivial element.")
    print(f"   - New distinct elements found: {len(new_elements_from_fig8)}")
    print("-" * 70)

    # --- Final Calculation ---
    print("Summary of contributions:")
    count_unknot = len(new_elements_from_unknot)
    count_trefoil = len(new_elements_from_trefoil)
    count_fig8 = len(new_elements_from_fig8)
    total_elements = len(elements_found)

    print(f"The total number of distinct elements is the sum of new elements from each knot type:")
    print(f"{count_unknot} (from Unknot) + {count_trefoil} (from Trefoil) + {count_fig8} (from Figure-eight) = {total_elements}")
    
    print("\nThe final answer is:")
    print(total_elements)

solve_homology_cobordism_count()