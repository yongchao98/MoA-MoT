def solve_diagonal_corridor_mate():
    """
    Determines the minimum pieces for a specific Diagonal Corridor Mate.

    The problem asks for the minimum number of additional pieces to create
    a Diagonal Corridor Mate with the White King on a1 and the Black King on h8.
    White is delivering the mate.

    1.  The Attacking Piece: The check must be along the a1-h8 diagonal.
        A White Bishop is the lowest-value piece that can do this.
    2.  The Blocking Pieces: The Black King's escape squares are h7 and g8.
        In this specific mate, they must be blocked by Black's own pieces.
        Two Black Pawns are the lowest-value pieces for this job.

    This results in a 3-piece solution: 1 White Bishop and 2 Black Pawns.
    This is the minimal solution both in number of pieces and total piece value.
    """

    # List the required pieces for each side.
    # These are already sorted by value (as there's only one type per side).
    white_pieces = ["White Bishop"]
    black_pieces = ["Black Pawn", "Black Pawn"]

    # Combine the lists as per the output specification (White first, then Black).
    final_piece_list = white_pieces + black_pieces

    # Format the list into a comma-separated string.
    # The instruction "output each number in the final equation" seems to be a
    # leftover from a different problem template, as there is no equation here.
    # The primary instruction and example specify a string list of piece names.
    # For example: White Bishop (value: 3), Black Pawn (value: 1), Black Pawn (value: 1).
    # The final output string will just be the names.
    final_answer = ", ".join(final_piece_list)

    print(final_answer)

solve_diagonal_corridor_mate()