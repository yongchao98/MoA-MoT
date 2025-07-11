def calculate_equivalence_classes():
    """
    This script calculates the number of equivalence classes for a given topological space X.

    The equivalence relation is defined as x ~ y if there is an auto-homeomorphism
    of the space sending x to y.

    The space X is the disjoint union of:
    - The torus
    - The sphere
    - The real line
    - A three-point discrete space
    - A five-point discrete space
    """

    # Step 1: Determine the number of equivalence classes for each component space.
    # A space where any point can be mapped to any other point by an auto-homeomorphism
    # is called homogeneous. Such a space has exactly 1 equivalence class.

    # The Torus, Sphere, and Real Line are all homogeneous spaces.
    num_classes_torus = 1
    num_classes_sphere = 1
    num_classes_real_line = 1

    # Any discrete space is homogeneous because any permutation of its points
    # is a homeomorphism.
    num_classes_3_point_discrete = 1
    num_classes_5_point_discrete = 1

    # Step 2: Sum the classes from each component.
    # Since the five component spaces are topologically distinct, any auto-homeomorphism
    # of the total space must map each component space to itself.
    # Therefore, the total number of equivalence classes is the sum of the classes
    # from each component.
    
    classes_per_space = [
        num_classes_torus,
        num_classes_sphere,
        num_classes_real_line,
        num_classes_3_point_discrete,
        num_classes_5_point_discrete
    ]
    
    total_classes = sum(classes_per_space)

    # Step 3: Print the explanation and the final result.
    print("The total number of equivalence classes is the sum of the classes from each of the 5 component spaces.")
    print("Each component space is homogeneous, containing only one equivalence class.")
    print("\nThe calculation is:")

    # Create and print the equation string.
    equation = " + ".join(map(str, classes_per_space)) + f" = {total_classes}"
    print(equation)


if __name__ == "__main__":
    calculate_equivalence_classes()