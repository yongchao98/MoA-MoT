def solve_equivalence_classes():
    """
    Calculates the number of equivalence classes for the given topological space X.

    The space X is the disjoint union of:
    - The torus
    - The sphere
    - The real line
    - A three-point discrete space
    - A five-point discrete space

    Two points x, y are equivalent if there is an auto-homeomorphism of X sending x to y.
    """

    # Step 1: Justification for decomposing the problem.
    # The five spaces are the connected components of X. We check if they are
    # topologically distinct.
    # - Torus vs Sphere: Different fundamental groups.
    # - Torus/Sphere vs Real Line: Compactness (Torus/Sphere are compact, Real Line is not).
    # - Torus/Sphere/Real Line vs Discrete Spaces: Connectivity (the first three are
    #   connected, discrete spaces with more than one point are not).
    # - 3-point discrete vs 5-point discrete: Cardinality (different number of points).
    # Since all five components are topologically distinct, any auto-homeomorphism of X
    # must map each component to itself. This means equivalence classes are contained
    # entirely within a single component. Thus, we can sum the number of equivalence
    # classes from each component.

    # Step 2: Determine the number of equivalence classes for each component.
    # This is equivalent to finding the number of orbits under the action of the
    # auto-homeomorphism group of each component. A space with only one orbit is
    # called homogeneous.

    # The torus is a homogeneous space. Any point can be mapped to any other
    # point via a translation, which is a homeomorphism.
    classes_torus = 1

    # The sphere is a homogeneous space. Any point can be mapped to any other
    # point via a rotation, which is a homeomorphism.
    classes_sphere = 1

    # The real line is a homogeneous space. Any point can be mapped to any other
    # point via a translation f(z) = z + c, which is a homeomorphism.
    classes_real_line = 1

    # A discrete space with n points has the symmetric group S_n as its
    # auto-homeomorphism group. The action of S_n on the n points is transitive.
    # Therefore, all points are in the same orbit.
    classes_discrete_3_point = 1
    classes_discrete_5_point = 1

    # Step 3: Calculate the total number of equivalence classes.
    component_classes = [
        classes_torus,
        classes_sphere,
        classes_real_line,
        classes_discrete_3_point,
        classes_discrete_5_point
    ]
    
    total_classes = sum(component_classes)

    # Output the result as an equation.
    class_strings = [str(c) for c in component_classes]
    equation = " + ".join(class_strings)
    
    print("The number of equivalence classes for the components are:")
    print(f"Torus: {classes_torus}")
    print(f"Sphere: {classes_sphere}")
    print(f"Real line: {classes_real_line}")
    print(f"Three-point discrete space: {classes_discrete_3_point}")
    print(f"Five-point discrete space: {classes_discrete_5_point}")
    print("\nThe total number of equivalence classes is the sum of the classes from each component:")
    print(f"{equation} = {total_classes}")

solve_equivalence_classes()