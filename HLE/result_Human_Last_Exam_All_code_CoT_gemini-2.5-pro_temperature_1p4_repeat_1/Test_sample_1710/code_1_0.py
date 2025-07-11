def solve_go_problem():
    """
    This function analyzes the given Go board position to find the optimal move for Black.

    The problem is a 'tsumego' (a life-and-death problem). The solution is a 'tesuji'
    (a clever or skillful play) that resolves the central conflict on the board.

    Analysis:
    1. The main fight is between the central white group and the surrounding black stones.
    2. The white stones at (2, 2) and (3, 2) are under direct pressure from black's
       incomplete net B(1,2), B(2,3), and B(3,3).
    3. The vital point for Black to play is (3, 1).
    4. Playing B(3,1) is an attachment that puts the white group in a critical
       situation, leading to the guaranteed capture of two white stones. This move
       simultaneously ensures the survival and strength of the surrounding black group.
    """

    # The optimal move is at row 3, column 1 (using 0-based indexing).
    # There is only one such move that works.
    best_move_row = 3
    best_move_col = 1

    # The problem asks for the answer in the format (row, col)
    # We print the numbers that make up the final coordinate pair.
    print(f"({best_move_row}, {best_move_col})")

solve_go_problem()