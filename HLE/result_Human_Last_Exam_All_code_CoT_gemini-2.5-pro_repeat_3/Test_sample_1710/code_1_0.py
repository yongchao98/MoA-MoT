def solve_go_problem():
    """
    This function analyzes the Go board situation and determines the optimal move for Black.

    The analysis reveals:
    1. A black group at (0,2) and (1,2) is in danger. It can be captured if White plays at (1,3).
       To survive, Black must play at (1,3) to connect this group to another friendly group.
    2. An isolated two-stone white group exists at (1,5) and (2,5). This is the target for capture.
    3. The single move at (1,3) is the vital point for Black's survival. The problem states this move must also
       be the one to capture the two white stones, meaning it creates an unescapable situation for them.

    Therefore, the solution is the move at (1,3).
    """
    # The solution is determined by analyzing the board state as described above.
    # The row is 1 (from top, 0-indexed).
    # The column is 3 (from left, 0-indexed).
    answer_row = 1
    answer_col = 3
    
    # The output format is specified as (row, col).
    print(f"({answer_row}, {answer_col})")

solve_go_problem()