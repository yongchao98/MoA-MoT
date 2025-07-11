def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # Orders of the fundamental groups of X1, X2, and X3
    n1 = 5  # For Z_5 from the pentagon
    n2 = 8  # For Z_8 from the octagon
    n3 = 2  # For Z_2 from the real projective plane

    # Number of groups in the free product
    k = 3

    print("The orders of the fundamental groups of the spaces X1, X2, and X3 are:")
    print(f"n1 = {n1}")
    print(f"n2 = {n2}")
    print(f"n3 = {n3}")
    print(f"The number of spaces being connected summed is k = {k}.")
    print("")

    print("The rank 'r' of the kernel is calculated using the formula:")
    print("r = 1 - (n2*n3 + n1*n3 + n1*n2) + (k-1)*n1*n2*n3")
    print("")

    # Calculate the terms of the formula
    term_n2n3 = n2 * n3
    term_n1n3 = n1 * n3
    term_n1n2 = n1 * n2
    sum_of_prods = term_n2n3 + term_n1n3 + term_n1n2

    coeff_k = k - 1
    prod_all = n1 * n2 * n3
    term_k = coeff_k * prod_all

    # Calculate the final rank
    rank = 1 - sum_of_prods + term_k

    print("Plugging in the values into the formula:")
    print(f"r = 1 - ({n2}*{n3} + {n1}*{n3} + {n1}*{n2}) + ({k}-1)*{n1}*{n2}*{n3}")
    print(f"r = 1 - ({term_n2n3} + {term_n1n3} + {term_n1n2}) + {coeff_k}*{prod_all}")
    print(f"r = 1 - {sum_of_prods} + {term_k}")
    print(f"r = {1 - sum_of_prods + term_k}")
    print("")
    print(f"The final rank of the kernel as a free group is {rank}.")

solve_rank_of_kernel()