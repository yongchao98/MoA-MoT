def count_homeomorphism_classes():
    """
    This function lists and counts the known homeomorphism classes of
    homogeneous planar continua based on the classification theorem.
    """

    # According to the mathematical classification theorem, these are the only classes.
    classes = {
        "The point": 1,
        "The circle": 1,
        "The pseudo-arc": 1,
        "The circle of pseudo-arcs": 1
    }

    print("The distinct homeomorphism classes of homogeneous planar continua are:")
    for class_name in classes.keys():
        print(f"- {class_name}")

    # The total number is the sum of these distinct classes.
    # The user asked to output each number in the final equation.
    # We represent this by summing 1 for each class.
    class_counts = list(classes.values())
    total_count = sum(class_counts)

    # Building the equation string: 1 + 1 + 1 + 1 = 4
    equation_str = " + ".join(map(str, class_counts))
    
    print("\nThe total number is found by summing the count of each class:")
    print(f"{equation_str} = {total_count}")

# Execute the function to print the result.
count_homeomorphism_classes()