def solve_cap_set_bound():
    """
    Calculates a lower bound for the cap set size in dimension 8 using the product method,
    and then states the best-known lower bound from literature.
    """
    # Known optimal sizes of cap sets in lower dimensions
    r3 = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112,
    }

    # Calculate lower bounds for n=8 using product construction
    # We check all partitions of 8 into n1 + n2 where n1, n2 are in our dict.
    n = 8
    best_bound = 0
    best_pair = (0, 0)

    for n1 in range(1, n // 2 + 1):
        n2 = n - n1
        if n1 in r3 and n2 in r3:
            bound = r3[n1] * r3[n2]
            if bound > best_bound:
                best_bound = bound
                best_pair = (n1, n2)

    print("Investigating lower bounds for r_3(8) using product construction:")
    n1, n2 = best_pair
    val1, val2 = r3[n1], r3[n2]
    print(f"The best bound from this method comes from dimensions {n1} and {n2}:")
    print(f"r_3({n1}) * r_3({n2}) = {val1} * {val2} = {best_bound}")
    print("\nThis provides a lower bound of 448.")

    # State the best-known lower bound from literature
    best_known_lower_bound = 496
    print(f"However, a better lower bound has been found using more complex constructions.")
    print(f"The best-known lower bound for the size of a cap set in dimension 8 is {best_known_lower_bound}.")

solve_cap_set_bound()