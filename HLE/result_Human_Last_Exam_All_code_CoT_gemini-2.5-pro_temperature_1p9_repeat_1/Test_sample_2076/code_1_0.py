def solve_chess_puzzle():
    """
    This function determines the minimum pieces for a Diagonal Corridor Mate
    with the White King on a1 and the Black King on h8, with White to mate.
    """
    
    # Pieces are sorted by value: Pawn (1), Knight (3), Bishop (3), Rook (5), Queen (9)
    # The minimum set of white pieces found are a Pawn and a Bishop.
    white_pieces = ["White Pawn", "White Bishop"]
    
    # The minimum set of black pieces found are two Pawns.
    black_pieces = ["Black Pawn", "Black Pawn"]
    
    # Combine the lists as per the problem description.
    all_pieces = white_pieces + black_pieces
    
    # Format the final output string.
    result_string = ", ".join(all_pieces)
    
    print(result_string)

solve_chess_puzzle()