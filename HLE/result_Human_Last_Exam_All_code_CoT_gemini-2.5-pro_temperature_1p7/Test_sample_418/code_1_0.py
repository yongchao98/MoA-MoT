def solve_equivalence_classes():
    """
    Calculates the total number of equivalence classes for the given topological space.

    The space X is a disjoint union of 5 topologically distinct components.
    An auto-homeomorphism of X must map each component to itself.
    Therefore, the total number of equivalence classes is the sum of the
    number of equivalence classes within each component.

    Each of the 5 component spaces is homogeneous, meaning that for any two
    points x, y in the component, there exists an auto-homeomorphism f
    of that component sending x to y. A homogeneous space has exactly one
    equivalence class.
    """

    # Number of classes for each component space
    # 1. The Torus
    torus_classes = 1
    # 2. The Sphere
    sphere_classes = 1
    # 3. The Real Line
    real_line_classes = 1
    # 4. A three-point discrete space
    three_point_discrete_classes = 1
    # 5. A five-point discrete space
    five_point_discrete_classes = 1

    class_counts = [
        torus_classes,
        sphere_classes,
        real_line_classes,
        three_point_discrete_classes,
        five_point_discrete_classes,
    ]

    total_classes = sum(class_counts)

    # Building the equation string to display the calculation
    # "1 + 1 + 1 + 1 + 1 = 5"
    equation_str = " + ".join(map(str, class_counts))
    final_equation = f"{equation_str} = {total_classes}"

    print("The number of equivalence classes for each topologically distinct component is 1, as each component is a homogeneous space.")
    print("The total number of equivalence classes is the sum of these counts.")
    print("The calculation is as follows:")
    print(final_equation)

solve_equivalence_classes()