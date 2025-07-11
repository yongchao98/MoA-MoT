def calculate_t4():
    """
    Calculates the number of ways to tile a 2x4 board using 2x1, 2x2, and 2x4 tiles.
    """
    # T_n represents the number of ways to tile a 2xn board.
    # We will use a dictionary to store the values of T_n.
    T = {}

    # Base cases determined by enumeration
    T[0] = 1  # For a 2x0 board, there is one way (the empty tiling).
    T[1] = 1  # For a 2x1 board, one 2x1 tile.
    T[2] = 2  # For a 2x2 board, two 2x1 tiles or one 2x2 tile.
    T[3] = 3  # For a 2x3 board, (1,1,1), (1,2), (2,1).

    # The recurrence relation for n >= 4 is T(n) = T(n-1) + T(n-2) + T(n-4).
    # We apply this for n=4.
    T[4] = T[3] + T[2] + T[0]

    print("To calculate T(4), we use the recurrence relation T(n) = T(n-1) + T(n-2) + T(n-4).")
    print("First, we need the values for T(3), T(2), and T(0).")
    print(f"T(0) = {T[0]}")
    print(f"T(1) = {T[1]}")
    print(f"T(2) = {T[2]}")
    print(f"T(3) = {T[3]}")
    print("\nNow, we can calculate T(4):")
    print(f"T(4) = T(3) + T(2) + T(0)")
    print(f"T(4) = {T[3]} + {T[2]} + {T[0]} = {T[4]}")

calculate_t4()