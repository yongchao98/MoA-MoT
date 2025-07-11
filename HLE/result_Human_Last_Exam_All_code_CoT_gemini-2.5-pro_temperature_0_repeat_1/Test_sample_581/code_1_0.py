def find_cap_set_lower_bound():
    """
    This function calculates a lower bound for the size of a cap set in dimension 8
    using the product construction and known values for smaller dimensions.
    It then states the best-known lower bound from mathematical literature.
    """
    # Known exact sizes of cap sets for smaller dimensions
    # r3_n represents r_3(n), the size of the largest cap set in dimension n.
    r3_2 = 4
    r3_6 = 112

    # The product construction gives a lower bound: r_3(n+m) >= r_3(n) * r_3(m)
    # We can find a lower bound for dimension 8 by using n=6 and m=2.
    dimension_n = 6
    dimension_m = 2
    total_dimension = dimension_n + dimension_m

    # Calculate the lower bound
    lower_bound = r3_6 * r3_2

    print("A simple method to find a lower bound for a cap set's size is the product construction.")
    print(f"The formula is: r_3(n+m) >= r_3(n) * r_3(m).")
    print(f"Using known values for dimension {dimension_n} (size={r3_6}) and dimension {dimension_m} (size={r3_2}):")
    print(f"r_3({total_dimension}) >= r_3({dimension_n}) * r_3({dimension_m})")
    print(f"So, a lower bound is {r3_6} * {r3_2} = {lower_bound}.")
    print("\nWhile 448 is a valid lower bound, a better one has been found using more advanced methods.")

    # State the best-known lower bound from literature
    best_known_lower_bound = 496
    print(f"The best known lower bound for the size of a cap set in dimension 8 is {best_known_lower_bound}.")

find_cap_set_lower_bound()