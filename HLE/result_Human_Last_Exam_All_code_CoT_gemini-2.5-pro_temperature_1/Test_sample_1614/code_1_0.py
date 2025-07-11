def calculate_tiling_ways():
    """
    Calculates the number of ways to tile a 2x4 board using 2x1, 2x2, and 2x4 tiles.
    """
    # T is a dictionary to store the values of T_n.
    T = {}

    # Base case T_0: One way to tile a 2x0 board (with no tiles).
    T[0] = 1

    # Base case T_1: One way to tile a 2x1 board (one 2x1 tile).
    T[1] = 1

    # Base case T_2: Two ways to tile a 2x2 board (two 2x1 tiles or one 2x2 tile).
    T[2] = 2

    # Calculate T_3 using the recurrence T_n = T_{n-1} + T_{n-2}.
    # The 2x4 tile is too large for a 2x3 board.
    T[3] = T[2] + T[1]

    # Calculate T_4 using the full recurrence relation: T_n = T_{n-1} + T_{n-2} + T_{n-4}.
    T[4] = T[3] + T[2] + T[0]

    # Print the step-by-step calculation leading to T_4.
    print(f"To find T_4, we first need T_0, T_1, T_2, and T_3.")
    print(f"T_0 = {T[0]}")
    print(f"T_1 = {T[1]}")
    print(f"T_2 = {T[2]}")
    print(f"T_3 = T_2 + T_1 = {T[2]} + {T[1]} = {T[3]}")
    print(f"Finally, we can calculate T_4:")
    # The final equation with numbers, as requested.
    print(f"T_4 = T_3 + T_2 + T_0 = {T[3]} + {T[2]} + {T[0]} = {T[4]}")


calculate_tiling_ways()