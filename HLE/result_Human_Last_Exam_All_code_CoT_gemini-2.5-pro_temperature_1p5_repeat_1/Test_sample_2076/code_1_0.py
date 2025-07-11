def solve_chess_puzzle():
    """
    Constructs and prints the list of pieces for the Diagonal Corridor Mate problem.
    """
    # The solution requires one White piece and two Black pieces.
    
    # White pieces list, sorted by value.
    # A Bishop is required to deliver the check along the diagonal.
    # It has a value of 3.
    white_pieces = ["White Bishop"]

    # Black pieces list, sorted by value.
    # Two Pawns are needed to block the king's escape squares (g7, h7).
    # Pawns are the lowest value piece (1).
    black_pieces = ["Black Pawn", "Black Pawn"]

    # Combine the lists as specified: White pieces first, then Black.
    all_pieces = white_pieces + black_pieces

    # Create the final comma-separated string for the output.
    result_string = ", ".join(all_pieces)
    
    print(result_string)

solve_chess_puzzle()
<<<White Bishop, Black Pawn, Black Pawn>>>