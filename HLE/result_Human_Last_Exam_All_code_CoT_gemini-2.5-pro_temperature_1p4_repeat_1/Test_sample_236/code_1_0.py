import collections

def solve_knot_surgery_classes():
    """
    Determines how many elements of the homology cobordism group can be represented
    by an integral surgery on a knot with at most four crossings.
    """

    # Step 1: Define the knots and their relevant properties.
    # Knots considered: Unknot (0_1), Trefoil (3_1), Figure-eight (4_1).
    # We must also consider the mirror of the non-amphichiral trefoil.
    knots = {
        "0_1": {"name": "Unknot"},
        "3_1": {"name": "Right-handed Trefoil"},
        "-3_1": {"name": "Left-handed Trefoil"},
        "4_1": {"name": "Figure-eight knot"},
    }

    # Step 2 & 3: Map surgeries on these knots to their resulting 3-manifolds and
    # identify their homology cobordism class.
    # The class is given a unique identifier based on established topological results.
    # 'P' denotes the Poincaré homology sphere Σ(2,3,5). Its class has order 2.
    # 'S7' denotes the Brieskorn sphere Σ(2,3,7). Its class has infinite order.
    #
    # Manifold properties:
    # [S^3]: The trivial class (identity).
    # [P]: Order 2, so [P] = [-P]. Non-trivial.
    # [S7]: Infinite order, so [S7] != [-S7]. Non-trivial.
    # All non-trivial classes listed here ([P], [S7], [-S7]) are distinct.
    surgery_results = {
        # (Knot ID, Surgery Coeff) -> (Resulting Manifold, Class Identifier)
        ("0_1", 1): ("S^3", "Identity"),
        ("0_1", -1): ("S^3", "Identity"),
        ("3_1", 1): ("Σ(2,3,7)", "Sigma(2,3,7)"),
        ("3_1", -1): ("Poincaré Sphere", "Poincare"),
        ("-3_1", 1): ("-Poincaré Sphere", "Poincare"), # class of -P is P
        ("-3_1", -1): ("-Σ(2,3,7)", "-Sigma(2,3,7)"),
        ("4_1", 1): ("Poincaré Sphere", "Poincare"),
        ("4_1", -1): ("-Poincaré Sphere", "Poincare"), # class of -P is P
    }

    print("Analyzing homology cobordism classes from knot surgeries:")
    print("=" * 60)

    represented_classes = set()
    surgeries_by_knot = collections.defaultdict(list)

    # Populate the surgeries for each knot
    for (knot_id, n), (manifold, class_id) in surgery_results.items():
        surgeries_by_knot[knot_id].append(
            (n, manifold, class_id)
        )
        represented_classes.add(class_id)

    # Print the analysis for each knot
    for knot_id in ["0_1", "3_1", "-3_1", "4_1"]:
        print(f"Knot: {knots[knot_id]['name']} ({knot_id})")
        for n, manifold, class_id in surgeries_by_knot[knot_id]:
            print(f"  -> {n: >+2}-surgery yields {manifold:<20} -> class '{class_id}'")
        print("-" * 60)

    # Step 4 & 5: Count the distinct classes found.
    count = len(represented_classes)
    class_list = sorted(list(represented_classes))
    
    print("The set of distinct homology cobordism classes found is:")
    print(class_list)
    
    # Create the final equation string
    equation_str = " + ".join(["1" for _ in class_list])

    print("\nFinal count:")
    print("Each unique class contributes 1 to the total count.")
    print(f"The calculation is: {equation_str} = {count}")

solve_knot_surgery_classes()
<<<4>>>