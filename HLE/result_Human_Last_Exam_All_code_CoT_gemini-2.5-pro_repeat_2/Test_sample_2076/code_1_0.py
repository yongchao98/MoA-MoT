def solve_diagonal_corridor_mate():
    """
    Determines the minimum pieces for a Diagonal Corridor Mate with Kings on a1 and h8.

    This function calculates the optimal set of additional pieces to create the checkmate,
    adhering to the principles of minimum piece count and minimum piece value.
    The result is then printed in the specified comma-separated format.
    """

    # White's attacking piece: A Bishop is the lowest-value piece that can deliver mate along a long diagonal.
    white_pieces = ["White Bishop"]

    # Black's blocking pieces: The Black King on h8 has three escape squares (g8, h7, g7).
    # These must be blocked by its own pieces. Pawns are the lowest-value pieces to do this.
    black_pieces = ["Black Pawn", "Black Pawn", "Black Pawn"]

    # Combine the lists as per the output requirements (White first, then Black).
    final_piece_list = white_pieces + black_pieces

    # Format the list into a comma-separated string for the final output.
    # The problem statement implies an equation, which we can represent as the list of pieces.
    # We will print each piece name as requested.
    output_string = ", ".join(final_piece_list)

    print(output_string)

solve_diagonal_corridor_mate()