def solve_diagonal_corridor_mate():
    """
    This function constructs the solution for the Diagonal Corridor Mate puzzle
    based on the specified constraints and prints the result.

    The puzzle requires finding the minimum number of additional pieces for a
    Diagonal Corridor Mate with the White King on a1 and the Black King on h8,
    with White delivering the mate.

    1.  The minimum number of pieces is determined to be 2. A White Queen on e5
        delivers check and controls two escape squares (g7, g8). A Black Pawn on h7
        blocks the final escape square.
    2.  White Pieces: A single White Queen.
    3.  Black Pieces: A single Black Pawn.
    4.  The output must be a comma-separated list, with White's pieces listed first,
        followed by Black's pieces.
    """
    white_pieces = ["White Queen"]
    black_pieces = ["Black Pawn"]

    # Combine the lists as per the output requirements (White first, then Black)
    final_piece_list = white_pieces + black_pieces

    # Format the list into a comma-separated string for printing
    output_string = ", ".join(final_piece_list)

    print(output_string)

solve_diagonal_corridor_mate()