import math

def calculate_blocks():
    """
    Calculates the number of blocks of kG by finding the number of conjugacy classes of S.

    The number of blocks of kG is equal to the number of blocks of kS,
    which is equal to the number of conjugacy classes of S.
    """
    # S is the extraspecial group of order 27, S = 3^(1+2)_+
    order_S = 27

    # The center Z(S) of an extraspecial p-group p^(1+2n) has order p.
    # For S, p=3.
    order_Z = 3

    # Each element of the center forms its own conjugacy class.
    num_central_classes = order_Z
    print(f"The group S has order {order_S}.")
    print(f"The center Z(S) has order {order_Z}.")
    print(f"This gives {num_central_classes} conjugacy classes of size 1.")

    # Calculate the number of non-central elements.
    num_non_central_elements = order_S - order_Z

    # For any non-central element g, the size of its centralizer C_S(g) is 9.
    order_centralizer_non_central = 9
    
    # The size of a conjugacy class for a non-central element is |S| / |C_S(g)|.
    size_non_central_class = order_S // order_centralizer_non_central
    
    # The number of non-central classes is the number of non-central elements
    # divided by the size of each non-central class.
    num_non_central_classes = num_non_central_elements // size_non_central_class
    print(f"The number of non-central elements is {num_non_central_elements}.")
    print(f"Each non-central element belongs to a conjugacy class of size {size_non_central_class}.")
    print(f"This gives {num_non_central_classes} non-central conjugacy classes.")

    # Total number of classes is the sum of central and non-central classes.
    total_classes = num_central_classes + num_non_central_classes
    
    print("\nThe total number of conjugacy classes of S is the final answer.")
    print(f"Total classes = {num_central_classes} + {num_non_central_classes} = {total_classes}")

calculate_blocks()