def solve_diagonal_corridor_mate():
    """
    This function determines the minimal pieces for a specific Diagonal Corridor Mate
    and prints the result in the required format.
    """

    # Based on the analysis, the minimal set of additional pieces is:
    # White: A Bishop to deliver check and a Pawn to protect it.
    # Black: Two Pawns to form the corridor trapping their own king.

    # List the pieces for each side.
    # Piece values: Pawn=1, Bishop=3. The lists are already sorted by value.
    white_pieces = ["Pawn", "Bishop"]
    black_pieces = ["Pawn", "Pawn"]

    # Format the piece names with their color.
    formatted_pieces = []
    for piece in white_pieces:
        formatted_pieces.append(f"White {piece}")
    for piece in black_pieces:
        formatted_pieces.append(f"Black {piece}")

    # Join the list into a single comma-separated string.
    final_answer = ", ".join(formatted_pieces)

    print(final_answer)

solve_diagonal_corridor_mate()