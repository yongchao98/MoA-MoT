import sys

def solve_topology_problem():
    """
    Calculates the number of equivalence classes for a disjoint union of topological spaces.

    The equivalence relation x ~ y holds if there is an auto-homeomorphism
    of the space sending x to y.
    """

    # An auto-homeomorphism must preserve connected components. Therefore, we can find the
    # number of equivalence classes for each component space and sum them.

    # 1. The torus is a homogeneous space. Any point can be mapped to any other
    # by a combination of rotations/translations.
    torus_classes = 1

    # 2. The sphere is also homogeneous. Any point can be mapped to any other via rotation.
    sphere_classes = 1

    # 3. The real line is homogeneous. A translation h(z) = z + a is a homeomorphism.
    real_line_classes = 1

    # 4. In a discrete space, any permutation of points is a homeomorphism.
    # Thus, a discrete space is always homogeneous.
    three_point_discrete_classes = 1

    # 5. Same logic applies to the five-point discrete space.
    five_point_discrete_classes = 1

    # The total number of equivalence classes is the sum of the classes from each component.
    all_classes = [
        torus_classes,
        sphere_classes,
        real_line_classes,
        three_point_discrete_classes,
        five_point_discrete_classes
    ]
    total_classes = sum(all_classes)

    print("Each of the five component spaces is homogeneous, meaning all points within that")
    print("component belong to a single equivalence class. The total number of classes is the")
    print("sum of the number of classes from each component.\n")
    print(f"Number of classes in the torus: {torus_classes}")
    print(f"Number of classes in the sphere: {sphere_classes}")
    print(f"Number of classes in the real line: {real_line_classes}")
    print(f"Number of classes in the 3-point discrete space: {three_point_discrete_classes}")
    print(f"Number of classes in the 5-point discrete space: {five_point_discrete_classes}")
    print("-" * 20)
    
    # Constructing the equation string as requested
    equation_str = " + ".join(map(str, all_classes))
    print(f"Total number of equivalence classes = {equation_str} = {total_classes}")

if __name__ == "__main__":
    solve_topology_problem()
