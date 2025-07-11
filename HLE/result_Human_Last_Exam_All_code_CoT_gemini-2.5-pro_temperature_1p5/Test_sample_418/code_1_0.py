def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space.
    """

    # The problem asks for the number of equivalence classes on a space X,
    # which is the disjoint union of five specific topological spaces.
    # The equivalence relation is defined by the action of the group of
    # auto-homeomorphisms of X.

    # Step 1: Define the subspaces and their properties.
    # Each dictionary entry represents a subspace and its number of equivalence classes.
    # A space is 'homogeneous' if any point can be mapped to any other point
    # via an auto-homeomorphism. Homogeneous spaces have exactly one equivalence class.
    spaces = {
        "The torus": {
            "explanation": "Homogeneous (points are equivalent via translations)",
            "classes": 1
        },
        "The sphere": {
            "explanation": "Homogeneous (points are equivalent via rotations)",
            "classes": 1
        },
        "The real line": {
            "explanation": "Homogeneous (points are equivalent via translations)",
            "classes": 1
        },
        "A three point discrete space": {
            "explanation": "Homogeneous (any permutation of points is a homeomorphism)",
            "classes": 1
        },
        "A five point discrete space": {
            "explanation": "Homogeneous (any permutation of points is a homeomorphism)",
            "classes": 1
        }
    }

    # Step 2: Explain the reasoning.
    # The five subspaces are topologically distinct (due to properties like
    # compactness, connectedness, fundamental group, and cardinality).
    # Therefore, any auto-homeomorphism of the total space X must map each
    # subspace to itself.
    # This means the total number of equivalence classes is the sum of the
    # number of classes in each individual subspace.
    print("The total number of equivalence classes is the sum of the classes from each subspace.")
    print("Number of classes per subspace:")

    total_classes = 0
    class_counts = []

    for name, properties in spaces.items():
        num_classes = properties["classes"]
        print(f"- {name}: {num_classes} class(es).")
        total_classes += num_classes
        class_counts.append(str(num_classes))

    # Step 3: Print the final calculation as an equation.
    equation = " + ".join(class_counts)
    print(f"\nFinal calculation: {equation} = {total_classes}")


if __name__ == "__main__":
    solve_equivalence_classes()