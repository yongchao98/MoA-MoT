def calculate_rank():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # Orders of the fundamental groups of X1, X2, and X3
    n1 = 5
    n2 = 8
    n3 = 2

    # The rank of the kernel K is given by the formula for the rank of the
    # Cartesian subgroup of a free product of finite abelian groups:
    # rank = 1 + (n1*n2*n3) - (n2*n3 + n1*n3 + n1*n2)

    # Calculate the product of all orders
    prod_all = n1 * n2 * n3

    # Calculate the sum of partial products
    sum_partial_prods = (n2 * n3) + (n1 * n3) + (n1 * n2)

    # Calculate the rank
    rank = 1 + prod_all - sum_partial_prods

    # Print the equation with all the numbers
    print(f"Rank = 1 + ({n1} * {n2} * {n3}) - ({n2} * {n3} + {n1} * {n3} + {n1} * {n2})")
    print(f"Rank = 1 + {prod_all} - ({n2*n3} + {n1*n3} + {n1*n2})")
    print(f"Rank = 1 + {prod_all} - {sum_partial_prods}")
    print(f"Rank = {1 + prod_all - sum_partial_prods}")

calculate_rank()