def solve_chess_puzzle():
    """
    This function determines the minimal pieces for a specific Diagonal Corridor Mate.

    The problem asks to construct a checkmate position where:
    - White King is on a1, Black King is on h8.
    - It's a "Diagonal Corridor Mate" delivered by White.
    - The solution must use the minimum number of additional pieces.
    - If multiple solutions have the same piece count, the one with the lowest
      total piece value is chosen.

    Our analysis leads to the following conclusion:
    1. A White Bishop on g7 checks the Black King on h8. This piece also
       attacks g8 and occupies g7, neutralizing two of the three escape squares.
       A Bishop is chosen over a Queen for its lower piece value.
    2. A Black Pawn on h7 blocks the final escape square for the Black King.
       A Pawn is chosen as it's the lowest-value piece that can perform this block.

    This two-piece solution (White Bishop, Black Pawn) is minimal in both
    piece count and total value.
    """
    # List of additional pieces for the solution.
    # The lists are pre-sorted by piece value (Pawn < Knight/Bishop < Rook < Queen).
    white_pieces = ["White Bishop"]
    black_pieces = ["Black Pawn"]

    # Combine the piece lists into a single list for output.
    final_piece_list = white_pieces + black_pieces

    # Format the list into the required comma-separated string.
    output = ", ".join(final_piece_list)

    print(output)

solve_chess_puzzle()