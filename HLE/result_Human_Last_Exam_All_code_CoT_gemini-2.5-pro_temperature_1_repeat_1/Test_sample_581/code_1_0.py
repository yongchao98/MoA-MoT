def solve_cap_set_bound():
    """
    This function explains and calculates a lower bound for the size of a cap set
    in dimension 8, and then provides the best-known current bound.
    """
    print("To find a lower bound for the size of a cap set in dimension 8, we can use the product construction.")
    print("This method combines known cap sets from smaller dimensions.")
    print("Let r_3(n) be the maximum size of a cap set in dimension n.")
    print("The product construction tells us that r_3(n1 + n2) >= r_3(n1) * r_3(n2).\n")

    # Known exact sizes for some smaller dimensions
    # r_3(6) = 112 was found by computer search by Davis and Maclagan (2003)
    # r_3(2) = 4 is a simple, well-known result.
    r3_6 = 112
    r3_2 = 4

    # We want a bound for dimension 8, which can be written as 6 + 2.
    n1 = 6
    n2 = 2
    n_target = n1 + n2

    # Calculate the lower bound from this construction
    constructed_bound = r3_6 * r3_2

    print(f"Using the known sizes for dimension n={n1} (size={r3_6}) and n={n2} (size={r3_2}), we can construct a lower bound for dimension n={n_target}:")
    # Output the final equation with each number, as requested.
    print(f"The equation is: {r3_6} * {r3_2} = {constructed_bound}")
    print(f"This construction proves that a cap set of at least size {constructed_bound} exists in dimension 8.\n")
    
    print("However, better constructions have been found.")
    best_known_bound = 496
    print(f"The best-known lower bound for the size of a cap set in dimension 8 is currently {best_known_bound}.")
    print("This result was established by C. Elsholtz and P.P. Palincza in 2020.")

solve_cap_set_bound()