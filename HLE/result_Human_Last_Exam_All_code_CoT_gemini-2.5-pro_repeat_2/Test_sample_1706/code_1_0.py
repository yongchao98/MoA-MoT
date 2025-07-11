def solve_go_problem():
    """
    This function determines the correct move for Black to capture all White stones
    and prints the coordinates of that move.
    """

    # Based on the analysis, the vital point to capture the white group is (2, 4).
    # This move splits the white stones and creates an immediate atari,
    # allowing Black to maintain the initiative and capture all stones.
    best_move_row = 2
    best_move_col = 4

    # The problem asks to output each number in the final move.
    print("The best first move for Black is to play at the coordinate below.")
    print(f"Row: {best_move_row}")
    print(f"Column: {best_move_col}")

solve_go_problem()