def solve_chess_puzzle():
    """
    This function determines the minimum pieces for a Diagonal Corridor Mate
    with the Black King on h8 and the White King on a1, and formats the answer.
    """

    # The minimal set of pieces for the pattern is one White Bishop and two Black Pawns.
    # White pieces are listed first, then Black pieces.
    # Each sub-list is sorted by piece value (Bishop=3, Pawn=1).
    
    white_pieces = ["White Bishop"]
    black_pieces = ["Black Pawn", "Black Pawn"]
    
    # Combine the lists to form the final answer
    all_pieces = white_pieces + black_pieces
    
    # Format the result as a comma-separated string
    final_answer = ", ".join(all_pieces)
    
    print(final_answer)

solve_chess_puzzle()