def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space X.

    The space X is the disjoint union of:
    - The torus
    - The sphere
    - The real line
    - A three-point discrete space
    - A five-point discrete space

    The equivalence relation is x ~ y if there is an auto-homeomorphism of X sending x to y.
    """

    # Step 1: Analyze the components.
    # The five components are the torus, the sphere, the real line, a 3-point discrete space,
    # and a 5-point discrete space.
    # No two of these spaces are homeomorphic to each other. For example:
    # - Torus and Sphere have different fundamental groups.
    # - The real line is non-compact, while the others are compact (except the line).
    # - The discrete spaces are not connected, while the others are.
    # - The 3-point and 5-point discrete spaces have different numbers of points.
    #
    # Since no two components are homeomorphic, any auto-homeomorphism of the total space X
    # must map each component to itself. This means that two points in different components
    # can never be in the same equivalence class.
    #
    # Therefore, the total number of equivalence classes is the sum of the number of
    # equivalence classes within each component.

    # Step 2: Determine the number of equivalence classes in each component.
    # A space is homogeneous if for any two points x, y, there is an auto-homeomorphism
    # mapping x to y. A homogeneous space has only one equivalence class.

    # The torus is a homogeneous space (e.g., via translations).
    torus_classes = 1

    # The sphere is a homogeneous space (e.g., via rotations).
    sphere_classes = 1

    # The real line is a homogeneous space (e.g., via translations).
    real_line_classes = 1

    # A discrete space with n points is homogeneous because any permutation of the points
    # is a homeomorphism.
    three_point_discrete_classes = 1
    five_point_discrete_classes = 1

    # Step 3: Sum the number of classes from each component.
    total_classes = (torus_classes +
                     sphere_classes +
                     real_line_classes +
                     three_point_discrete_classes +
                     five_point_discrete_classes)

    # Step 4: Print the results, including the final equation.
    print("The total number of equivalence classes is the sum of the classes in each component.")
    print(f"Number of classes in the torus: {torus_classes}")
    print(f"Number of classes in the sphere: {sphere_classes}")
    print(f"Number of classes in the real line: {real_line_classes}")
    print(f"Number of classes in the three-point discrete space: {three_point_discrete_classes}")
    print(f"Number of classes in the five-point discrete space: {five_point_discrete_classes}")
    
    # The final equation as requested.
    print("\nThe final calculation is:")
    print(f"{torus_classes} + {sphere_classes} + {real_line_classes} + {three_point_discrete_classes} + {five_point_discrete_classes} = {total_classes}")

solve_equivalence_classes()