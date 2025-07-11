def calculate_product_bound():
    """
    Calculates lower bounds for cap set sizes using the product construction.
    """
    # Best known lower bounds for r_3(n) for small n.
    # For n=1..6, these are the exact values.
    r3_bounds = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112
    }

    print("Calculating product construction bound for r_3(7):")
    n = 7
    max_bound_n7 = 0
    # Iterate through partitions of 7 (a, b) where a+b=7
    for a in range(1, (n // 2) + 1):
        b = n - a
        if a in r3_bounds and b in r3_bounds:
            bound = r3_bounds[a] * r3_bounds[b]
            print(f"r_3({a}) * r_3({b}) = {r3_bounds[a]} * {r3_bounds[b]} = {bound}")
            if bound > max_bound_n7:
                max_bound_n7 = bound
    
    r3_bounds[n] = max_bound_n7
    print(f"\nBest product construction bound for r_3(7) is {max_bound_n7}. Adding this to our known bounds.")
    
    print("\n-------------------------------------\n")
    
    print("Calculating product construction bound for r_3(8):")
    n = 8
    max_bound_n8 = 0
    # Iterate through partitions of 8 (a, b) where a+b=8
    for a in range(1, (n // 2) + 1):
        b = n - a
        if a in r3_bounds and b in r3_bounds:
            bound = r3_bounds[a] * r3_bounds[b]
            print(f"r_3({a}) * r_3({b}) = {r3_bounds[a]} * {r3_bounds[b]} = {bound}")
            if bound > max_bound_n8:
                max_bound_n8 = bound

    print(f"\nThe best lower bound for r_3(8) from the product construction is {max_bound_n8}.")
    print("\nHowever, the best *known* lower bound from any construction method must be considered.")
    
    best_known_bound = 496
    print(f"The current record lower bound for r_3(8) is {best_known_bound}.")

calculate_product_bound()