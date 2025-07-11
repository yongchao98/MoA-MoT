def solve_go_problem():
    """
    This function analyzes the Go board situation and determines the optimal move for Black.

    The board state is:
    1 EWBEEEBWE
    2 EWBEEWBWE
    3 EWWBEWBWE
    4 EEWBEEBWE
    5 EEWWBBBWE
    6 EEEEWWWWE
    7 EEEEEEEEE
    8 EEEEEEEEE

    It is Black's turn.

    Analysis:
    1. Survival: The Black stones at (3, 4) and (4, 4) are weak and under attack.
       To survive, they must be connected to a stronger group. Placing a stone at (4, 5)
       connects the stone at (4, 4) to the strong Black chain starting at (5, 5).
       This single move solidifies the entire central Black group, ensuring its survival.

    2. Capture: This move also creates a decisive attack. The two White stones at (3, 3)
       and (4, 3), which were previously attacking Black, are now cut off and their
       position becomes untenable. The move at (4, 5) ensures they can be captured.

    Therefore, the vital point for both survival and capture is the same.
    The coordinates are 1-indexed from the top-left corner.
    """
    # The solution is the point that ensures survival for the weak black group
    # and simultaneously sets up a capture of the two attacking white stones.
    # This vital point is at row 4, column 5.
    best_move_row = 4
    best_move_col = 5

    # The problem asks for the answer in the format (row, col)
    # and to order them lexicographically if there are multiple.
    # In this classic tsumego, there is only one optimal point.
    answer = f"({best_move_row}, {best_move_col})"
    print(answer)

solve_go_problem()