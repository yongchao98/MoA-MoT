def solve_topology_problem():
    """
    Calculates the number of equivalence classes for the given topological space.

    The space X is a disjoint union of five components. An auto-homeomorphism
    of X must map each component to a homeomorphic component. We first determine
    that none of the given components are homeomorphic to each other, so any
    auto-homeomorphism maps each component to itself.

    This means the total number of equivalence classes is the sum of the
    number of equivalence classes within each component.

    We then analyze each component:
    1. Torus: Homogeneous space. All points are equivalent. 1 class.
    2. Sphere: Homogeneous space. All points are equivalent. 1 class.
    3. Real line: Homogeneous space. All points are equivalent. 1 class.
    4. 3-point discrete space: Homogeneous, as any permutation is a homeomorphism. 1 class.
    5. 5-point discrete space: Homogeneous, for the same reason. 1 class.

    The total is the sum of the classes from each component.
    """

    # Number of equivalence classes for each component space
    classes_per_component = {
        "The torus": 1,
        "The sphere": 1,
        "The real line": 1,
        "A three point discrete space": 1,
        "A five point discrete space": 1,
    }

    # Extract the numbers for the final equation
    numbers = list(classes_per_component.values())

    # Calculate the total number of equivalence classes
    total_classes = sum(numbers)

    # Format the numbers into a string representing the sum
    equation_str = " + ".join(map(str, numbers))

    # Print the final equation and the result
    print(f"The number of equivalence classes is the sum of the classes in each non-homeomorphic component.")
    print(f"Equation: {equation_str} = {total_classes}")

solve_topology_problem()