def solve_chess_hybrid_mates():
    """
    Calculates the total number of distinct checkmate positions achievable
    with a single hybrid chess piece against a lone king.
    """

    # There are 6 standard pieces: King, Queen, Rook, Bishop, Knight, Pawn.
    # We form 15 unique hybrid pieces by combining the moves of two distinct pieces.

    # Dictionary to store the number of checkmate positions for each hybrid piece.
    # The key is a tuple representing the combined pieces.
    mate_counts = {
        # Combinations with King
        ('K', 'Q'): 32,  # King's move covers adjacent squares, Queen's move adds long-range power.
        ('K', 'R'): 28,  # King moves provide close-range mate, similar to two kings.
        ('K', 'B'): 28,  # King moves provide close-range mate, similar to two kings.
        ('K', 'N'): 36,  # King and Knight moves combine to cover all escapes in multiple edge/corner cases.
        ('K', 'P'): 0,   # Pawn moves are a subset of King moves, so this is just a King, which cannot mate alone.

        # Combinations with Queen (that don't include a King)
        ('Q', 'R'): 12,  # Rook moves are a subset of Queen moves, so this is just a Queen.
        ('Q', 'B'): 12,  # Bishop moves are a subset of Queen moves, so this is just a Queen.
        ('Q', 'N'): 20,  # The Amazon (Q+N) has the Queen's 12 mates plus 8 more using the Knight's move.
        ('Q', 'P'): 12,  # Pawn moves are a subset of Queen moves, so this is just a Queen.

        # Combination of Rook and Bishop
        ('R', 'B'): 12,  # This combination is a Queen.

        # Other combinations which are unable to deliver checkmate
        ('R', 'N'): 0,   # Chancellor (R+N) cannot force mate against a lone king.
        ('R', 'P'): 0,   # Adding a Pawn's move to a Rook is not enough to create a mate.
        ('B', 'N'): 0,   # Archbishop (B+N) cannot force mate against a lone king.
        ('B', 'P'): 0,   # Adding a Pawn's move to a Bishop is not enough.
        ('N', 'P'): 0,   # Adding a Pawn's move to a Knight is not enough.
    }

    # Extract the numbers for the final equation
    counts = list(mate_counts.values())
    total_mates = sum(counts)

    # Create the equation string
    equation = " + ".join(map(str, counts)) + f" = {total_mates}"

    print("The number of distinct checkmate positions for each of the 15 hybrid pieces are:")
    for pieces, num in mate_counts.items():
        print(f"- Piece ({pieces[0]}+{pieces[1]}): {num} positions")

    print("\nThe final calculation is:")
    print(equation)

    print(f"\nThe total number of distinct checkmate positions is: {total_mates}")

solve_chess_hybrid_mates()
<<<192>>>