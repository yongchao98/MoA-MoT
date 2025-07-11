def solve_chess_puzzle():
    """
    This function determines the pieces needed for the specified checkmate
    and prints them in the required format.
    """

    # Based on the analysis, the minimal set of additional pieces is:
    # - White: 1 Bishop
    # - Black: 1 Pawn
    # The lists are already sorted by piece value (ascending).
    white_pieces = ["White Bishop"]
    black_pieces = ["Black Pawn"]

    # Combine the lists, with White's pieces first.
    all_pieces = white_pieces + black_pieces

    # Format the output as a comma-separated string.
    result_string = ", ".join(all_pieces)

    print(result_string)

solve_chess_puzzle()