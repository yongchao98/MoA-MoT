def solve_cap_set_bound():
    """
    Calculates some simple lower bounds for the cap set problem in dimension 8
    and then states the best-known lower bound.
    """
    # The size of the largest cap set in dimension n, r_3(n), is known for n <= 6.
    r3_known = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112
    }

    print("The cap set problem seeks the maximum number of points in GF(3)^n with no three on a line.")
    print("Lower bounds can be found by constructing large cap sets.")
    print("A simple method is the product construction: r_3(a+b) >= r_3(a) * r_3(b).\n")
    print("Let's calculate some lower bounds for dimension 8 using this method:")

    # Calculate r_3(8) >= r_3(4) * r_3(4)
    dim1, dim2 = 4, 4
    val1, val2 = r3_known[dim1], r3_known[dim2]
    bound1 = val1 * val2
    print(f"Using r_3({dim1}) * r_3({dim2}): {val1} * {val2} = {bound1}")

    # Calculate r_3(8) >= r_3(5) * r_3(3)
    dim1, dim2 = 5, 3
    val1, val2 = r3_known[dim1], r3_known[dim2]
    bound2 = val1 * val2
    print(f"Using r_3({dim1}) * r_3({dim2}): {val1} * {val2} = {bound2}")

    # Calculate r_3(8) >= r_3(6) * r_3(2)
    dim1, dim2 = 6, 2
    val1, val2 = r3_known[dim1], r3_known[dim2]
    bound3 = val1 * val2
    print(f"Using r_3({dim1}) * r_3({dim2}): {val1} * {val2} = {bound3}")

    print("\nWhile these product constructions provide good lower bounds, they are not the best known.")
    print("More advanced constructions yield a better result.")

    # The best known lower bound for r_3(8)
    best_known_lower_bound = 496
    print(f"\nThe best known lower bound for the size of a cap set in dimension 8 is {best_known_lower_bound}.")

solve_cap_set_bound()