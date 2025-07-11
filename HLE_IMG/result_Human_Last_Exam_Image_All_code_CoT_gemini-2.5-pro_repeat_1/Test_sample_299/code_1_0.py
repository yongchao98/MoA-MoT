def solve_chess_puzzle():
    """
    This function prints the solution to the mate-in-2 chess puzzle.
    White's first move, Qa4-d1+, forces checkmate on the second move
    regardless of Black's reply. This function outlines the two possible
    lines of play.
    """
    first_move_white = "Qa4-d1+"

    print("The key move for White is 1. " + first_move_white)
    print("This forces checkmate in two different ways depending on Black's response:")
    print("")

    # Variation 1: If Black King moves to d3
    black_response_1 = "Kd3"
    second_move_white_1 = "Qd1-b3#"
    print("Variation 1:")
    print("1. " + first_move_white + " " + black_response_1)
    print("2. " + second_move_white_1)
    print("")

    # Variation 2: If Black King captures the Rook on f3
    black_response_2 = "Kxf3"
    second_move_white_2 = "Ne7-g6#"
    print("Variation 2:")
    print("1. " + first_move_white + " " + black_response_2)
    print("2. " + second_move_white_2)

solve_chess_puzzle()