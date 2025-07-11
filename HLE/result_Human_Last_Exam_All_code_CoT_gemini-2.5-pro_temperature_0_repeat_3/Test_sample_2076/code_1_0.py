def solve_diagonal_corridor_mate():
    """
    Constructs the piece list for a minimal Diagonal Corridor Mate.

    The problem asks for a minimal position for a Diagonal Corridor Mate
    with the White King on a1 and the Black King on h8, where White delivers the mate.

    1.  Attacker: A White Bishop is the lowest-value piece (3) that can attack
        h8 along the long diagonal.
    2.  Blockers: The Black King's escape squares (g8, h7) must be blocked by
        Black's own pieces (g7 is covered by the White Bishop).
        - A Black Pawn (1) is the cheapest piece to block h7.
        - A Black Knight (3) is the cheapest piece (besides a pawn, which is illegal)
          to block g8.
    3.  The final list of additional pieces is sorted by color (White then Black)
        and then by value.
    """

    # Define the additional pieces required for the mate
    white_pieces = [
        ("White Bishop", 3)
    ]
    black_pieces = [
        ("Black Pawn", 1),
        ("Black Knight", 3)
    ]

    # Sort each list by piece value (ascending)
    white_pieces.sort(key=lambda x: x[1])
    black_pieces.sort(key=lambda x: x[1])

    # Combine the lists, with White's pieces first
    final_piece_names = [p[0] for p in white_pieces] + [p[0] for p in black_pieces]

    # Format the result as a comma-separated string
    result_string = ", ".join(final_piece_names)

    print(result_string)

solve_diagonal_corridor_mate()