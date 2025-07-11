def calculate_tiling_ways():
    """
    Calculates T_4, the number of ways to tile a 2x4 board with
    2x1, 2x2, and 2x4 tiles, by using a recurrence relation.
    """
    # T_n represents the number of ways to tile a 2xn board.
    # The recurrence relation is T_n = T_{n-1} + T_{n-2} + T_{n-4}.

    # We first define the base cases.
    # T[n] will store the value of T_n.
    T = {}

    # T_0: One way to tile a 2x0 board (the empty tiling).
    T[0] = 1

    # T_1: One way to tile a 2x1 board (one 2x1 tile).
    T[1] = 1
    
    # T_2: Two ways for a 2x2 board (two 2x1 tiles, or one 2x2 tile).
    # This can also be calculated as T_2 = T_1 + T_0
    T[2] = T[1] + T[0]

    # T_3: Ways to tile a 2x3 board. We can place a 2x1 tile (leaving a 2x2 board)
    # or a 2x2 tile (leaving a 2x1 board). So T_3 = T_2 + T_1.
    T[3] = T[2] + T[1]

    # Now we calculate T_4 using the full recurrence: T_4 = T_3 + T_2 + T_0
    t3 = T[3]
    t2 = T[2]
    t0 = T[0]
    t4 = t3 + t2 + t0

    # Print the explanation and the final equation.
    print(f"To calculate T_4, we use the recurrence relation T_n = T_(n-1) + T_(n-2) + T_(n-4).")
    print(f"For n=4, this is T_4 = T_3 + T_2 + T_0.")
    print(f"First, we find the base cases:")
    print(f"T_0 = {t0}")
    print(f"T_1 = {T[1]}")
    print(f"T_2 = T_1 + T_0 = {T[1]} + {t0} = {t2}")
    print(f"T_3 = T_2 + T_1 = {t2} + {T[1]} = {t3}")
    print(f"Plugging these values into the equation for T_4:")
    print(f"T_4 = {t3} + {t2} + {t0} = {t4}")

calculate_tiling_ways()
<<<6>>>