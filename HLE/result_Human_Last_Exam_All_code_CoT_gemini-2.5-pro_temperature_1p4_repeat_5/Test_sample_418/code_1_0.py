def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space.

    The space is a disjoint union of five components. The number of equivalence classes
    for the total space is the sum of the number of equivalence classes for each component,
    as no homeomorphism can map points between these topologically distinct components.

    Each of the component spaces (Torus, Sphere, Real Line, discrete spaces) is
    homogeneous, meaning all its points belong to a single equivalence class.
    """

    # Number of equivalence classes for the Torus
    num_classes_torus = 1

    # Number of equivalence classes for the Sphere
    num_classes_sphere = 1

    # Number of equivalence classes for the Real Line
    num_classes_real_line = 1

    # Number of equivalence classes for a three-point discrete space
    num_classes_discrete_3 = 1

    # Number of equivalence classes for a five-point discrete space
    num_classes_discrete_5 = 1

    # Total number of equivalence classes is the sum of the classes from each component
    total_classes = (num_classes_torus +
                     num_classes_sphere +
                     num_classes_real_line +
                     num_classes_discrete_3 +
                     num_classes_discrete_5)

    # Print the equation as requested
    print(f"{num_classes_torus} (Torus) + "
          f"{num_classes_sphere} (Sphere) + "
          f"{num_classes_real_line} (Real Line) + "
          f"{num_classes_discrete_3} (3-pt Discrete) + "
          f"{num_classes_discrete_5} (5-pt Discrete) = {total_classes}")

solve_equivalence_classes()