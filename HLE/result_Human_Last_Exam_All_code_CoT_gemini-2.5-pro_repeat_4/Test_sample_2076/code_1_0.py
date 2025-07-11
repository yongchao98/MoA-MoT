def solve_diagonal_corridor_mate():
    """
    This function determines the minimum pieces for a Diagonal Corridor Mate
    and prints them in the specified format.
    """

    # Pieces are determined by the logic described above.
    # White pieces are listed first, sorted by value (Pawn=1, Bishop=3).
    white_pieces = ["White Pawn", "White Bishop"]

    # Black pieces are listed next, sorted by value (Pawn=1).
    black_pieces = ["Black Pawn", "Black Pawn"]

    # Combine the lists to form the final set of additional pieces.
    final_pieces = white_pieces + black_pieces

    # Format the list as a comma-separated string for the output.
    # This represents the "final equation" of pieces needed.
    final_answer = ", ".join(final_pieces)

    print(final_answer)

solve_diagonal_corridor_mate()