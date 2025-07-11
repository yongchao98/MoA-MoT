def solve_cap_set_bound():
    """
    This function explains and identifies the best-known lower bound for the size
    of a cap set in dimension 8.
    """
    dimension = 8

    # A simple lower bound is constructed from vectors with only {0, 1} coordinates.
    # The size of this set is 2^n.
    simple_bound = 2 ** dimension
    print(f"A basic construction for a cap set in dimension {dimension} gives a lower bound of 2^{dimension} = {simple_bound}.")

    # Better bounds are found using the product construction with known smaller values.
    # The size of the largest cap set in dimension 6, r_3(6), is 112.
    # The size of the largest cap set in dimension 2, r_3(2), is 4.
    r3_6 = 112
    r3_2 = 4
    product_bound = r3_6 * r3_2
    print(f"The product construction using dimensions 6 and 2 gives a better bound:")
    print(f"r_3({dimension}) >= r_3(6) * r_3(2) = {r3_6} * {r3_2} = {product_bound}")

    # The best-known lower bound comes from a specific construction by Yves Edel (2004).
    best_known_bound = 496
    print(f"\nHowever, the best-known lower bound, established by a specific construction, is larger.")
    print(f"The record for dimension {dimension} is {best_known_bound}.")

solve_cap_set_bound()