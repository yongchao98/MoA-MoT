def solve_go_problem():
    """
    This function determines the correct move for Black on the given Go board.

    The logic is as follows:
    1. The two white stones to be captured are at (1, 5) and (2, 5) (0-indexed).
    2. Their liberties are at (0, 5) and (1, 4). Black must play at one of these points.
    3. Playing at B(0, 5) allows White to play at W(1, 4) and create a living group with 4 liberties. This fails.
    4. Playing at B(1, 4) forces White to respond at W(0, 5). The resulting white group can be systematically captured through a correct sequence of moves by Black.
    5. Therefore, the unique winning move for Black is at (1, 4).
    6. We convert this 0-indexed coordinate to a 1-indexed format for the final answer.
    """
    # The correct move, identified through analysis (0-indexed)
    correct_move_row_0_indexed = 1
    correct_move_col_0_indexed = 4

    # Convert to 1-indexed for the output format
    final_row = correct_move_row_0_indexed + 1
    final_col = correct_move_col_0_indexed + 1

    # Print the final answer in the format (row, col)
    # The problem statement mentions "output each number in the final equation",
    # so we format the string carefully.
    print(f"({final_row}, {final_col})")

solve_go_problem()