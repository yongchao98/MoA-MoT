def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes on the given topological space X.

    The space X is the disjoint union of:
    - The torus
    - The sphere
    - The real line
    - A three-point discrete space
    - A five-point discrete space

    The equivalence relation x ~ y means there is an auto-homeomorphism of X sending x to y.
    """

    # An auto-homeomorphism of X must map each topologically distinct component to itself.
    # Therefore, we can find the number of equivalence classes on each component space
    # and sum them up.

    # 1. The Torus (T^2)
    # The torus is a homogeneous space. For any two points x, y on the torus, there is
    # a translation (which is a homeomorphism) that maps x to y.
    # This means all points on the torus are in the same equivalence class.
    classes_torus = 1
    print(f"Number of equivalence classes on the torus: {classes_torus}")

    # 2. The Sphere (S^2)
    # The sphere is also a homogeneous space. For any two points x, y on the sphere, there
    # is a rotation (a homeomorphism) that maps x to y.
    # All points on the sphere are in one equivalence class.
    classes_sphere = 1
    print(f"Number of equivalence classes on the sphere: {classes_sphere}")

    # 3. The Real Line (R)
    # The real line is a homogeneous space. The translation f(z) = z + (y-x) is a
    # homeomorphism mapping point x to point y.
    # All points on the real line are in one equivalence class.
    classes_real_line = 1
    print(f"Number of equivalence classes on the real line: {classes_real_line}")

    # 4. A three-point discrete space
    # In a discrete space, any permutation of the points is a homeomorphism. The group
    # of permutations acts transitively, so any point can be mapped to any other.
    # All three points are in one equivalence class.
    classes_discrete_3 = 1
    print(f"Number of equivalence classes on the 3-point discrete space: {classes_discrete_3}")

    # 5. A five-point discrete space
    # Similar to the three-point space, the space is homogeneous under its group
    # of auto-homeomorphisms (permutations).
    # All five points are in one equivalence class.
    classes_discrete_5 = 1
    print(f"Number of equivalence classes on the 5-point discrete space: {classes_discrete_5}")

    # Total number of classes is the sum of the classes from each component.
    total_classes = (classes_torus + classes_sphere + classes_real_line +
                     classes_discrete_3 + classes_discrete_5)

    print("\nTotal number of equivalence classes is the sum from each component:")
    print(f"{classes_torus} + {classes_sphere} + {classes_real_line} + {classes_discrete_3} + {classes_discrete_5} = {total_classes}")

solve_equivalence_classes()