def solve_go_problem():
    """
    This function determines the correct move for Black on the given Go board.
    The goal is to capture the two white stones and ensure the black group's survival.

    The two white stones are at (0, 1) and (1, 1).
    Their liberties are at (0, 0) and (1, 0).

    1. If Black plays at (0, 0), it starts a ladder that fails because White can escape to the right.
    2. If Black plays at (1, 0), it creates a "net" (geta). White's only escape attempt is to play at (0, 0),
       but this results in the white stones being captured immediately.

    Therefore, the only correct move is (1, 0).
    """
    # The correct row and column for Black's move
    correct_row = 1
    correct_col = 0

    # The problem asks to output the numbers in the final "equation" or answer.
    # We will format the output as (row, col).
    print(f"({correct_row}, {correct_col})")

solve_go_problem()