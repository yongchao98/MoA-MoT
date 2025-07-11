def solve_chess_puzzle():
    """
    This function determines the minimum pieces for a Diagonal Corridor Mate
    with the White King on a1 and the Black King on h8, and prints the result.
    """

    # The checkmate pattern requires a White piece checking along the a1-h8 diagonal,
    # and Black pieces blocking the King's escape squares (g8 and h7).

    # To minimize piece count and then value:
    # 1. White's checking piece: A Bishop (value 3) is cheaper than a Queen (value 9).
    # 2. Black's blocking pieces: Pawns (value 1) are the cheapest pieces to block the corridor.

    white_pieces = ["White Bishop"]
    black_pieces = ["Black Pawn", "Black Pawn"]

    # The problem asks for the pieces to be listed as a comma-separated string.
    # White pieces first, then Black pieces. Each sub-list is sorted by piece value.
    # Our lists are already sorted correctly.
    final_piece_list = white_pieces + black_pieces

    # The final answer is the string representation of this list.
    answer = ", ".join(final_piece_list)

    print(answer)

solve_chess_puzzle()
<<<White Bishop, Black Pawn, Black Pawn>>>