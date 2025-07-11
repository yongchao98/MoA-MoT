def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space X.

    The space X is a disjoint union of five components. An auto-homeomorphism of X
    must map each component to itself, as no two components are homeomorphic.
    Therefore, the total number of equivalence classes is the sum of the number of
    classes in each component.

    For each component, we determine if it is a "homogeneous space," meaning its
    group of auto-homeomorphisms acts transitively on it. If so, there is only
    one equivalence class.
    """

    # 1. Torus: Homogeneous space. For any x, y on the torus, there is a
    # translation (a homeomorphism) mapping x to y.
    torus_classes = 1

    # 2. Sphere: Homogeneous space. For any x, y on the sphere, there is a
    # rotation (a homeomorphism) mapping x to y.
    sphere_classes = 1

    # 3. Real line: Homogeneous space. For any x, y on the real line, there
    # is a translation (a homeomorphism) mapping x to y.
    real_line_classes = 1

    # 4. Three-point discrete space: Homogeneous space. Any permutation of the
    # points is a homeomorphism. The permutation group acts transitively.
    discrete_3_classes = 1

    # 5. Five-point discrete space: Homogeneous space. Similar to the 3-point
    # space, any permutation is a homeomorphism, and the group acts transitively.
    discrete_5_classes = 1

    class_counts = [
        torus_classes,
        sphere_classes,
        real_line_classes,
        discrete_3_classes,
        discrete_5_classes,
    ]

    total_classes = sum(class_counts)

    # Building the equation string
    equation_str = " + ".join(map(str, class_counts)) + f" = {total_classes}"

    print(f"Number of equivalence classes in the Torus: {torus_classes}")
    print(f"Number of equivalence classes in the Sphere: {sphere_classes}")
    print(f"Number of equivalence classes in the Real Line: {real_line_classes}")
    print(f"Number of equivalence classes in the 3-point discrete space: {discrete_3_classes}")
    print(f"Number of equivalence classes in the 5-point discrete space: {discrete_5_classes}")
    print("\nThe total number of equivalence classes is the sum:")
    print(equation_str)

solve_equivalence_classes()