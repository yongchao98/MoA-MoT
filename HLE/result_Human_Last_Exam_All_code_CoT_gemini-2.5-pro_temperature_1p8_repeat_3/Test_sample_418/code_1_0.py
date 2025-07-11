def solve_equivalence_classes():
    """
    Calculates the total number of equivalence classes for the given topological space.

    The space is a disjoint union of five components. We determined that any
    auto-homeomorphism of the total space must map each component to itself.
    Therefore, the total number of equivalence classes is the sum of the number
    of equivalence classes on each component space.

    For each of the five spaces, the group of auto-homeomorphisms acts transitively
    on its points, meaning each space is homogeneous and has only one equivalence class.
    """

    # Number of equivalence classes on the torus
    torus_classes = 1
    # Number of equivalence classes on the sphere
    sphere_classes = 1
    # Number of equivalence classes on the real line
    real_line_classes = 1
    # Number of equivalence classes on a three-point discrete space
    three_point_classes = 1
    # Number of equivalence classes on a five-point discrete space
    five_point_classes = 1

    # The total number of equivalence classes is the sum.
    total_classes = torus_classes + sphere_classes + real_line_classes + three_point_classes + five_point_classes

    print("The total number of equivalence classes is the sum of the classes on each component space.")
    print("Number of classes for each component:")
    print(f"- Torus: {torus_classes}")
    print(f"- Sphere: {sphere_classes}")
    print(f"- Real Line: {real_line_classes}")
    print(f"- Three-point discrete space: {three_point_classes}")
    print(f"- Five-point discrete space: {five_point_classes}")
    print("\nThe final equation is:")
    print(f"{torus_classes} + {sphere_classes} + {real_line_classes} + {three_point_classes} + {five_point_classes} = {total_classes}")

solve_equivalence_classes()