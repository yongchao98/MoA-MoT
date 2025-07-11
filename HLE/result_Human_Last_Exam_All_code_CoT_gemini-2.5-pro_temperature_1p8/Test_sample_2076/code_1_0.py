def solve_diagonal_corridor_mate():
    """
    This function determines and prints the pieces for the specified chess puzzle.

    The logic for finding the pieces is as follows:
    1.  The problem requires checkmating the Black King on h8 along the a1-h8 diagonal.
    2.  The most efficient way to both deliver check and control the King's escape squares (h7 and g8) is with a White Queen on g7.
    3.  The Queen on g7 is attacked by the King on h8, so it must be protected.
    4.  The piece with the lowest value that can protect the Queen is a White Pawn on h6.
    5.  This two-piece solution (White Queen, White Pawn) is the minimum number of pieces possible. It is also the lowest-value two-piece solution.
    6.  The final list of pieces is then formatted according to the rules: sorted by value, with White's pieces listed first.
    """
    
    # The minimal set of pieces required for the mate.
    white_pieces = ["Pawn", "Queen"]
    black_pieces = []

    # A dictionary of piece values to enable sorting.
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9
    }

    # Sort each list of pieces by their value in ascending order.
    white_pieces.sort(key=lambda piece: piece_values[piece])
    black_pieces.sort(key=lambda piece: piece_values[piece])
    
    # Prepare the final list for printing, adding the color prefix.
    final_list = ["White " + piece for piece in white_pieces] + ["Black " + piece for piece in black_pieces]

    # Join the list into a comma-separated string for the final output.
    result = ", ".join(final_list)
    
    print(result)

solve_diagonal_corridor_mate()
<<<White Pawn, White Queen>>>