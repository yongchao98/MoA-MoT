def solve_homology_cobordism_problem():
    """
    Calculates and explains the number of homology cobordism group elements
    represented by integral surgery on knots with at most four crossings.
    """

    print("To find the number of elements of the homology cobordism group from knots with at most four crossings, we proceed as follows:")
    print("\n1. Identify the knots: Unknot (0_1), Trefoil (3_1), and Figure-eight (4_1).")
    print("2. Consider integral surgery with framings +1 and -1, which produce homology spheres.")
    print("3. Identify the resulting distinct homology cobordism classes based on known theorems.")

    # These results are based on established theorems in knot theory.
    # The code serves to present the argument clearly.
    
    # Contribution from the Unknot (0_1)
    # Both +1 and -1 surgery on the unknot yield S^3, the trivial element.
    num_from_unknot = 1
    print(f"\n- Unknot (0_1): Surgeries result in the 3-sphere S^3.")
    print(f"  This represents the trivial element in the homology cobordism group.")
    print(f"  Number of distinct elements from Unknot = {num_from_unknot}")

    # Contribution from the Trefoil knot (3_1)
    # +1 surgery gives the Poincaré sphere (P), -1 surgery gives its mirror (-P).
    # These are two new, distinct elements.
    num_from_trefoil = 2
    print(f"\n- Trefoil (3_1): Surgeries result in the Poincaré sphere (P) and its mirror (-P).")
    print(f"  These are two distinct, non-trivial elements.")
    print(f"  Number of new distinct elements from Trefoil = {num_from_trefoil}")

    # Contribution from the Figure-eight knot (4_1)
    # Surgeries on 4_1 also produce P and -P, which we have already counted.
    num_new_from_fig8 = 0
    print(f"\n- Figure-eight (4_1): Surgeries result in the same two elements, P and -P.")
    print(f"  Therefore, there are no new elements from the figure-eight knot.")
    print(f"  Number of new distinct elements from Figure-eight = {num_new_from_fig8}")

    # Final Calculation
    total_elements = num_from_unknot + num_from_trefoil + num_new_from_fig8

    print("\n" + "="*40)
    print("Final Calculation:")
    print("The total number of distinct elements is the sum of the unique contributions from each knot type.")
    print(f"Total elements = (from Unknot) + (new from Trefoil) + (new from Figure-eight)")
    print(f"Total elements = {num_from_unknot} + {num_from_trefoil} + {num_new_from_fig8} = {total_elements}")
    print("="*40)

solve_homology_cobordism_problem()