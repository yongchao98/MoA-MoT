def solve_go_problem():
    """
    This function determines the optimal move(s) for Black on the given Go board.

    The analysis identifies two vital points for Black that both ensure the survival of
    Black's central group and lead to the capture of a two-stone white group.
    This is a "miai" situation, where Black has two equally good options.

    1. Move at (2, 4): This move threatens the white stones at (1,5) and (2,5)
       while strengthening the black stone at (2,3).
    2. Move at (3, 4): This move solidifies the black cutting stones at (3,3) and (4,4),
       which launches a decisive attack to capture the white stones at (4,2) and (4,3).

    Since both moves are valid answers, and the request is to provide both if they exist,
    we list them in lexicographically increasing order.
    """
    # The coordinates are 0-indexed (row, col).
    answer1_row = 2
    answer1_col = 4
    answer2_row = 3
    answer2_col = 4

    # Format the answers as specified.
    # The output format is (row, col), (row, col)
    print(f"({answer1_row}, {answer1_col}), ({answer2_row}, {answer2_col})")

solve_go_problem()