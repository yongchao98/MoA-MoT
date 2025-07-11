def solve_diagonal_corridor_mate():
    """
    This function determines the minimum pieces for a Diagonal Corridor Mate
    with the White King on a1 and the Black King on h8, and prints the result.
    """

    # To deliver the mate, White needs a piece checking the h8 King along the long diagonal.
    # A Bishop is the lowest-value piece that can do this.
    white_pieces = ["White Bishop"]

    # The Black King on h8 has three escape squares: g8, h7, and g7.
    # These must be blocked. The most efficient way is with Black's own lowest-value pieces (Pawns).
    black_pieces = ["Black Pawn", "Black Pawn", "Black Pawn"]

    # Combine the lists, with White's pieces first.
    all_pieces = white_pieces + black_pieces

    # Create the final comma-separated string.
    result_string = ", ".join(all_pieces)

    print(result_string)

solve_diagonal_corridor_mate()