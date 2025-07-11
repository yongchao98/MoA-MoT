def solve_cap_set_bound():
    """
    Explains and calculates lower bounds for the cap set problem in dimension 8.
    """
    # The cap set problem asks for the size of the largest subset of (Z/3)^n
    # without a 3-term arithmetic progression. This size is denoted r_3(n).
    # We are looking for the best-known lower bound for n=8.

    # Exact values are known for small n:
    r3_values = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112
    }

    # A good lower bound can be found using the product construction:
    # r_3(n+m) >= r_3(n) * r_3(m).
    # We can calculate a bound for r_3(8) using r_3(6) and r_3(2).
    n1 = 6
    n2 = 2
    r3_n1 = r3_values[n1]
    r3_n2 = r3_values[n2]
    product_bound = r3_n1 * r3_n2

    print("A method for finding a lower bound for r_3(8) is the product construction.")
    print(f"Using known values r_3({n1}) = {r3_n1} and r_3({n2}) = {r3_n2}:")
    print(f"A lower bound is r_3({n1}) * r_3({n2}) = {r3_n1} * {r3_n2} = {product_bound}")
    print("-" * 40)

    # However, this is not the best-known bound. More advanced, specific
    # constructions have produced better results. The current record comes
    # from mathematical literature, not a simple calculation.
    best_known_lower_bound = 496
    print("The best known lower bound for the size of a cap set in dimension 8 is:")
    print(best_known_lower_bound)

solve_cap_set_bound()