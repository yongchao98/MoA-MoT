def solve_diagonal_corridor_mate():
    """
    Determines the minimum pieces for a specific Diagonal Corridor Mate.

    The logic is as follows:
    1.  The Black King is on h8, checked along the a1-h8 diagonal.
    2.  The mating piece must be a White Bishop or Queen. A Bishop has the minimum value.
    3.  The King's escape squares are g8 and h7. These must be blocked by Black pieces.
    4.  To block h7, a Black Pawn is the cheapest option.
    5.  To block g8, a Black Pawn is illegal. A Black Knight is the cheapest piece that
        does not block the diagonal check from the White Bishop.
    6.  The final list is constructed and sorted according to the rules.
    """

    # Piece values for sorting (not strictly needed for this hardcoded solution, but good practice)
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9
    }

    # Minimum pieces determined by the logic above
    white_pieces = ["White Bishop"]
    black_pieces = ["Black Pawn", "Black Knight"]

    # The lists are already sorted by value:
    # White: Bishop (3)
    # Black: Pawn (1), Knight (3)

    # Combine the lists and format as a comma-separated string
    final_piece_list = white_pieces + black_pieces
    result_string = ", ".join(final_piece_list)

    print(result_string)

solve_diagonal_corridor_mate()
<<<White Bishop, Black Pawn, Black Knight>>>