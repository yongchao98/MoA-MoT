def solve():
    """
    Calculates the largest possible dimension for the quotient ring R/I.
    """
    # The order of the group G is given.
    group_order = 10000

    # The dimension of the coinvariant algebra R/I is the order of the group that
    # acts effectively on the polynomial ring. This effective group is the image of
    # the group homomorphism phi: G -> GL(10, C) that defines the action.
    # The dimension is given by |G| / |ker(phi)|.

    # To find the largest possible dimension, we need to find the smallest possible
    # size of the kernel of the action phi.
    # The smallest a subgroup can be is the trivial group, of size 1.
    # This corresponds to a faithful action.
    min_kernel_size = 1

    # This maximum is achievable if a group of order 10000 can have a faithful
    # 10-dimensional representation. The group (Z/10Z)^4 is an example
    # of a group of order 10000 with a faithful representation in GL(4,C),
    # and thus in GL(10,C).

    # Therefore, the largest possible dimension is the order of G.
    largest_dimension = group_order / min_kernel_size

    print("The largest possible dimension for R/I is calculated as follows:")
    print("Dimension = |G| / |ker(phi)|")
    print("To maximize the dimension, we minimize the kernel size.")
    print(f"The calculation for the largest dimension is:")
    print(f"{group_order} / {min_kernel_size} = {int(largest_dimension)}")

solve()