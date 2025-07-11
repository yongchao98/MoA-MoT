def solve_cap_set_bound():
    """
    Calculates a lower bound for the size of a cap set in dimension 8
    using the product construction and known values, then states the
    best-known lower bound from the literature.
    """
    print("The size of the largest cap set in the vector space (Z/3Z)^n is denoted by r_3(n).")
    print("The exact values of r_3(n) are only known for n <= 6.")
    print("We can find a lower bound for r_3(8) using the product construction rule: r_3(n+m) >= r_3(n) * r_3(m).\n")

    # Known exact values for r_3(n) for n=1 to 6
    r3_known = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112,
    }

    print("Known values:")
    for n_val, size in r3_known.items():
        print(f"r_3({n_val}) = {size}")

    n_target = 8
    best_product_bound = 0
    print(f"\nCalculating lower bounds for r_3({n_target}) using partitions of {n_target}:")

    # Iterate through partitions of n_target (e.g., 8 = 2+6, 3+5, 4+4)
    for i in range(1, n_target // 2 + 1):
        j = n_target - i
        if i in r3_known and j in r3_known:
            bound = r3_known[i] * r3_known[j]
            # Output each number in the final equation as requested
            print(f"For partition {i}+{j}: r_3({i}) * r_3({j}) = {r3_known[i]} * {r3_known[j]} = {bound}")
            if bound > best_product_bound:
                best_product_bound = bound

    print(f"\nThe best lower bound from this simple product construction is {best_product_bound}.")

    # State the actual best-known lower bound from research
    best_known_bound = 496
    print("\nHowever, more advanced constructions have produced a better result.")
    print(f"The best known lower bound for the size of a cap set in dimension 8 is {best_known_bound}.")

solve_cap_set_bound()