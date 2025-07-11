def solve_diagonal_corridor_mate():
    """
    This function determines the minimum pieces for a Diagonal Corridor Mate
    with the given constraints and prints the result.
    """

    # To satisfy the minimum piece value, White's checking piece is a Bishop.
    white_pieces = ["White Bishop"]

    # To satisfy the minimum piece count and value, Black's blocking pieces are Pawns.
    # The King on h8 needs its escape squares g8 and h7 blocked by its own pieces.
    # The escape square g7 is covered by the White Bishop's attack.
    black_pieces = ["Black Pawn", "Black Pawn"]

    # Sort each list by piece value (though not strictly necessary with these pieces)
    # Pawn=1, Knight=3, Bishop=3, Rook=5, Queen=9
    # The lists are already effectively sorted.
    
    # Combine the lists as required
    final_pieces = white_pieces + black_pieces
    
    # Format the final output string
    result_string = ", ".join(final_pieces)
    
    print(result_string)

solve_diagonal_corridor_mate()