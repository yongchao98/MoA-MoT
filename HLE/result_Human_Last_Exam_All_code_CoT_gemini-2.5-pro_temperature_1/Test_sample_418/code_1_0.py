def solve_topology_equivalence_classes():
    """
    Calculates the total number of equivalence classes for a disjoint union of topological spaces.

    The equivalence relation is x ~ y if there is an auto-homeomorphism sending x to y.
    For a disjoint union, auto-homeomorphisms map components to themselves.
    So, we sum the number of equivalence classes for each component.
    """

    # Each of the given spaces is homogeneous, meaning all points within that space
    # form a single equivalence class.
    classes_torus = 1
    classes_sphere = 1
    classes_real_line = 1
    classes_three_point_discrete = 1
    classes_five_point_discrete = 1

    # The total number of equivalence classes is the sum of the classes from each component.
    total_classes = (classes_torus +
                     classes_sphere +
                     classes_real_line +
                     classes_three_point_discrete +
                     classes_five_point_discrete)

    # Print the breakdown of the calculation.
    print(f"Number of equivalence classes in the Torus: {classes_torus}")
    print(f"Number of equivalence classes in the Sphere: {classes_sphere}")
    print(f"Number of equivalence classes in the Real Line: {classes_real_line}")
    print(f"Number of equivalence classes in the Three Point Discrete Space: {classes_three_point_discrete}")
    print(f"Number of equivalence classes in the Five Point Discrete Space: {classes_five_point_discrete}")
    print("\nTotal number of equivalence classes is the sum:")
    print(f"{classes_torus} + {classes_sphere} + {classes_real_line} + {classes_three_point_discrete} + {classes_five_point_discrete} = {total_classes}")

solve_topology_equivalence_classes()