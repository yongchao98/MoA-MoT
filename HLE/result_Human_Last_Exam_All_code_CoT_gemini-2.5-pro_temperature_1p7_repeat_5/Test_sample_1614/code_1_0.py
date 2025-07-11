def solve_tiling_problem():
    """
    This function calculates T_4, the number of ways to tile a 2x4 board
    using 2x1, 2x2, and 2x4 tiles, and prints the final calculation.
    """
    # A dictionary will be used to store the values of T_n for our dynamic programming approach.
    T = {}

    # Based on our analysis, we have the following recurrence relation:
    # T_n = T_{n-1} + T_{n-2} + T_{n-4}

    # We first establish the base cases by manual calculation.
    # T_0: There is one way to tile a 2x0 board (the empty tiling).
    T[0] = 1
    
    # T_1: A 2x1 board can only be tiled with one 2x1 tile.
    T[1] = 1
    
    # T_2: A 2x2 board can be tiled with two 2x1 tiles or one 2x2 tile.
    T[2] = 2
    
    # T_3: A 2x3 board can be tiled in T_2 + T_1 = 2 + 1 = 3 ways.
    T[3] = T[2] + T[1]

    # Now, we use the recurrence relation to calculate T_4.
    # T_4 = T_3 + T_2 + T_0
    T_4_value = T[3] + T[2] + T[0]

    # Finally, print the equation with each number used in the calculation, as requested.
    print(f"T_4 = {T[3]} + {T[2]} + {T[0]} = {T_4_value}")

solve_tiling_problem()