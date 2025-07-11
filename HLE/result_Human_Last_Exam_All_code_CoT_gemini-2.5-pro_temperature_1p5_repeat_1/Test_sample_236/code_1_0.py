def solve_homology_cobordism_count():
    """
    Calculates the number of homology cobordism group elements from surgery
    on knots with at most four crossings.
    """

    # Step 1: Define knots with at most 4 crossings.
    # Note: 3_1* is the mirror of 3_1. 0_1 and 4_1 are amphichiral.
    knots = {
        '0_1': 'Unknot',
        '3_1': 'Right-handed Trefoil',
        '3_1*': 'Left-handed Trefoil',
        '4_1': 'Figure-eight Knot'
    }
    print("This program calculates the number of distinct elements of the homology cobordism group")
    print("represented by integral surgery on knots with at most four crossings.\n")

    # Step 2 & 3: Identify resulting manifolds from +/-1 surgery.
    # Let E = S^3 (trivial element)
    # Let P = Poincaré homology sphere class, [Σ(2,3,5)]
    # Let T = [Σ(2,3,7)]
    # Let -T = [-Σ(2,3,7)]
    surgery_map = {
        ('0_1', 1): ('S^3', 'Trivial element'),
        ('0_1', -1): ('S^3', 'Trivial element'),
        ('3_1', 1): ('Σ(2,3,7)', 'Σ(2,3,7) class'),
        ('3_1', -1): ('Σ(2,3,5)', 'Poincaré sphere class'),
        ('3_1*', 1): ('-Σ(2,3,5)', 'Poincaré sphere class'), # Class has order 2
        ('3_1*', -1): ('-Σ(2,3,7)', 'Inverse of Σ(2,3,7) class'),
        ('4_1', 1): ('Σ(2,3,5)', 'Poincaré sphere class'),
        ('4_1', -1): ('-Σ(2,3,5)', 'Poincaré sphere class'), # Class has order 2
    }

    print("--- Surgeries and Resulting Manifolds ---")
    distinct_classes = set()
    for (knot, coeff), (manifold, class_desc) in surgery_map.items():
        print(f"Surgery on {knots[knot]:<22} with coefficient {coeff:>2} -> {manifold:<8} (represents: {class_desc})")
        distinct_classes.add(class_desc)

    # Step 4: Count the distinct classes.
    print("\n--- Counting Distinct Classes ---")
    print("The distinct classes obtained are:")
    class_list = sorted(list(distinct_classes))
    for i, class_name in enumerate(class_list):
        print(f"{i+1}. {class_name}")

    count = len(distinct_classes)

    print("\n--- Final Calculation ---")
    contributions = {
        "Trivial element": "Unknot",
        "Poincaré sphere class": "Trefoil and Figure-eight knots",
        "Σ(2,3,7) class": "Right-handed Trefoil",
        "Inverse of Σ(2,3,7) class": "Left-handed Trefoil"
    }

    print("The final count is the sum of these unique classes:")
    sum_str = " + ".join(["1" for _ in class_list])
    print(f"Total = {sum_str} = {count}")
    print("\nIn summary:")
    for class_name, source in contributions.items():
        if class_name in class_list:
            print(f"- 1 element from {source} ({class_name})")

solve_homology_cobordism_count()