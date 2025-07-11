def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space.

    The space X is a disjoint union of five components. The equivalence relation
    x ~ y holds if there is an auto-homeomorphism of X mapping x to y.

    The five components are topologically distinct, so any auto-homeomorphism of X
    maps each component to itself. This means the total number of equivalence
    classes is the sum of the number of classes in each component.

    Each of the five components is a homogeneous space, meaning all its points
    belong to a single equivalence class.
    """

    # Number of classes for the Torus (homogeneous)
    num_classes_torus = 1

    # Number of classes for the Sphere (homogeneous)
    num_classes_sphere = 1

    # Number of classes for the Real line (homogeneous)
    num_classes_real_line = 1

    # Number of classes for a 3-point discrete space (homogeneous)
    num_classes_discrete_3 = 1

    # Number of classes for a 5-point discrete space (homogeneous)
    num_classes_discrete_5 = 1

    # Total number of classes is the sum
    total_classes = (num_classes_torus +
                     num_classes_sphere +
                     num_classes_real_line +
                     num_classes_discrete_3 +
                     num_classes_discrete_5)

    print("Each of the five component spaces is homogeneous, meaning all points within")
    print("that component belong to a single equivalence class.")
    print("The number of classes for each component is therefore 1.")
    print("\nThe total number of equivalence classes is the sum of the classes from each component.")
    print("\nFinal Calculation:")
    print(f"{num_classes_torus} (torus) + "
          f"{num_classes_sphere} (sphere) + "
          f"{num_classes_real_line} (real line) + "
          f"{num_classes_discrete_3} (3-pt discrete) + "
          f"{num_classes_discrete_5} (5-pt discrete) = "
          f"{total_classes}")

solve_equivalence_classes()