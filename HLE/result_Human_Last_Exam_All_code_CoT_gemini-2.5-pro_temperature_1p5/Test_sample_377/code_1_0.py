def count_conjugacy_classes_of_extraspecial_group():
    """
    Calculates the number of conjugacy classes of the extraspecial group S = 3^(1+2).
    """
    p = 3
    n = 1

    # Order of the group S = p^(2n+1)
    order_S = p**(2 * n + 1)
    # Order of the center Z(S) = p
    order_Z_S = p
    
    # The center Z(S) contains p elements. Each central element forms a conjugacy class of size 1.
    num_central_classes = order_Z_S
    print(f"The order of the extraspecial group S is {order_S}.")
    print(f"The center Z(S) has {order_Z_S} elements.")
    print(f"Each central element forms its own conjugacy class.")
    print(f"Number of central classes: {num_central_classes}")
    
    # The remaining elements are non-central.
    num_non_central_elements = order_S - order_Z_S
    print(f"Number of non-central elements: {num_non_central_elements}")
    
    # For any non-central element g, the size of its centralizer C_S(g) is p^(2n).
    size_centralizer_non_central = p**(2 * n)
    
    # The size of any non-central conjugacy class is |S| / |C_S(g)| = p.
    size_non_central_class = order_S / size_centralizer_non_central
    print(f"The size of any non-central conjugacy class is {int(size_non_central_class)}.")
    
    # The number of non-central conjugacy classes is the number of non-central elements
    # divided by the size of each non-central class.
    num_non_central_classes = num_non_central_elements / size_non_central_class
    print(f"Number of non-central classes: {int(num_non_central_classes)}")
    
    # Total number of conjugacy classes is the sum of central and non-central classes.
    total_classes = num_central_classes + num_non_central_classes
    print(f"\nTotal number of conjugacy classes in S is the sum of central and non-central classes.")
    print(f"Total number of classes = {num_central_classes} + {int(num_non_central_classes)}")
    
    # According to the block theory reasoning, this is the number of blocks of kG.
    print(f"\nThe number of blocks of kG is equal to the number of conjugacy classes of S.")
    print(f"Number of blocks of kG = {int(total_classes)}")
    
    return int(total_classes)

# Execute the function
count_conjugacy_classes_of_extraspecial_group()
