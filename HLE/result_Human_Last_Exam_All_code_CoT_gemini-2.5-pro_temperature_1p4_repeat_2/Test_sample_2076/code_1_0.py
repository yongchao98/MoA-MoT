def solve_diagonal_corridor_mate():
    """
    This function determines and prints the minimal pieces required for a specific
    Diagonal Corridor Mate chess position.

    The problem asks for a minimal position where:
    - White King is on a1, Black King is on h8.
    - White delivers a "Diagonal Corridor Mate".
    - Minimality is first by number of pieces, then by total piece value.

    The analysis leads to the following conclusion:
    1.  The check must be delivered along a diagonal. The cheapest piece for this is a White Bishop.
        To make the check unblockable, we can place it at g7.
    2.  The King's escape squares (g8, h7) must be blocked by its own pieces to fit the
        "corridor" definition. The cheapest pieces for this are two Black Pawns on g8 and h7.
    3.  This results in 3 additional pieces: 1 White Bishop and 2 Black Pawns.
    4.  The piece values are: Bishop (3), Pawn (1), Pawn (1). The total value is 5, which is minimal.

    The final list is formatted as requested: White pieces first, then Black, with each
    group sorted by piece value.
    """

    # Define the pieces for the minimal solution
    white_pieces = ["White Bishop"]
    black_pieces = ["Black Pawn", "Black Pawn"]

    # The lists are already sorted by piece value (Bishop=3, Pawn=1).
    # Combine the lists for the final output.
    final_piece_list = white_pieces + black_pieces

    # Print the result as a comma-separated string.
    print(", ".join(final_piece_list))

solve_diagonal_corridor_mate()