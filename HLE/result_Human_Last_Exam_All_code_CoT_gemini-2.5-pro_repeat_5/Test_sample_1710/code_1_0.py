def solve_go_problem():
    """
    This function determines the optimal move for Black on the given Go board.

    The problem is a life-and-death situation. A black group on the right is in danger with only two liberties: (0, 5) and (3, 5).
    Black must play at one of these points to survive.
    The goal is to also capture the two white stones at (1, 5) and (2, 5).

    - If Black plays at (0, 5): The black group is saved, but White can then play at (3, 5) to connect their two stones to safety. The capture fails.
    - If Black plays at (3, 5): The black group is saved by connecting to another black group, making it very strong. This move also cuts the two white stones off from any help, leading to their eventual capture.

    Therefore, the unique correct move is (3, 5).
    """
    # The solution is found by logical deduction based on the rules of Go.
    # The optimal move is at row 3, column 5.
    row = 3
    col = 5

    # The final answer is formatted as (row, col)
    # The problem asks to output each number in the final equation.
    print(f"({row}, {col})")

solve_go_problem()