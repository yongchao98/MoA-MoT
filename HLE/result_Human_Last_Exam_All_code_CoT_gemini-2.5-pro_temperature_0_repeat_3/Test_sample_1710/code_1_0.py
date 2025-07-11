def solve_go_puzzle():
    """
    This function determines the correct move for Black on the given Go board.

    The problem is a capturing race on the right side of the board.
    - A black group {B(0,7), B(1,7), B(2,7), B(3,7)} has one liberty at (3,8).
    - A white group {W(0,8), W(1,8), W(2,8)} also has one liberty at (3,8).

    It is Black's turn. By playing at the shared liberty (3,8), Black captures
    the white group and saves its own group. This is the only move that
    satisfies the prompt's conditions of survival and capture.
    """
    # The solution is the coordinate of the critical liberty.
    # Row is 3, Column is 8 (using 0-based indexing).
    solution_coordinate = (3, 8)

    # Print the solution in the specified format: (row, col)
    print(f"({solution_coordinate[0]}, {solution_coordinate[1]})")

solve_go_puzzle()