def count_blocks():
    """
    This function calculates the number of blocks for the group algebra kG.
    Based on the theoretical argument, this number is equal to the number of
    conjugacy classes of the Sylow 3-subgroup S. S is an extraspecial
    p-group of order p^(2n+1) with p=3 and n=1.
    """
    p = 3
    n = 1
    
    # The center Z(S) gives p conjugacy classes of size 1.
    num_central_classes = p
    print(f"The group S has a center of order {p}. This gives {p} central conjugacy classes.")
    
    # The remaining elements are non-central.
    total_elements = p**(2*n + 1)
    num_non_central_elements = total_elements - p
    
    # Each non-central class has size p.
    size_non_central_class = p
    num_non_central_classes = num_non_central_elements // size_non_central_class
    print(f"The number of non-central elements is {total_elements} - {p} = {num_non_central_elements}.")
    print(f"These elements form classes of size {p}, so there are {num_non_central_elements} / {p} = {num_non_central_classes} non-central classes.")

    # Total number of classes.
    total_classes = num_central_classes + num_non_central_classes
    print(f"The total number of conjugacy classes in S is the sum:")
    print(f"{num_central_classes} + {num_non_central_classes} = {total_classes}")
    
    print("\nSince the number of blocks of kG is equal to the number of conjugacy classes of S,")
    print(f"the number of blocks is {total_classes}.")

count_blocks()