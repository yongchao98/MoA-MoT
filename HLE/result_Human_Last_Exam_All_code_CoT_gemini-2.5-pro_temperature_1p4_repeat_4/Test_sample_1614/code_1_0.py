def solve_tiling():
    """
    Calculates T_4, the number of ways to tile a 2x4 board
    with 2x1, 2x2, and 2x4 tiles.
    """
    # T_n represents the number of ways to tile a 2xn board.
    # We establish a recurrence relation: T_n = T_{n-1} + T_{n-2} + T_{n-4}

    # Base cases determined by enumeration:
    # T_0: 1 way (empty tiling)
    # T_1: 1 way (one 2x1 tile)
    # T_2: 2 ways (two 2x1 tiles, or one 2x2 tile)
    # T_3: 3 ways (three 2x1 tiles; or a 2x2 and a 2x1; or a 2x1 and a 2x2)
    T = {
        0: 1,
        1: 1,
        2: 2,
        3: 3
    }

    # Calculate T_4 using the recurrence relation:
    # T_4 = T_3 + T_2 + T_0
    t3 = T[3]
    t2 = T[2]
    t0 = T[0]
    t4 = t3 + t2 + t0

    # Print the final calculation step by step
    print(f"To find T_4, we use the recurrence relation T_n = T_(n-1) + T_(n-2) + T_(n-4).")
    print(f"For n=4, this gives T_4 = T_3 + T_2 + T_0.")
    print(f"Using the pre-calculated base cases:")
    print(f"T_4 = {t3} + {t2} + {t0} = {t4}")

solve_tiling()