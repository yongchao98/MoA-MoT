def solve_chess_puzzle():
    """
    This function constructs the solution to the Diagonal Corridor Mate puzzle
    based on the specified constraints.
    """

    # Piece values for calculation (not strictly needed for final output formatting but good for clarity)
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9
    }

    # The minimal set of additional pieces determined by the logic.
    # White pieces are sorted by value (Bishop and Knight are equal).
    # The prompt requires Black to have a material advantage.
    # White's value: Bishop (3) + Knight (3) = 6
    # Black's value: Queen (9) = 9
    # Since 9 > 6, this configuration is valid.
    white_pieces = ["White Bishop", "White Knight"]
    black_pieces = ["Black Queen"]

    # Combine the lists to form the final comma-separated string
    all_pieces = white_pieces + black_pieces
    result_string = ", ".join(all_pieces)

    print(result_string)

solve_chess_puzzle()