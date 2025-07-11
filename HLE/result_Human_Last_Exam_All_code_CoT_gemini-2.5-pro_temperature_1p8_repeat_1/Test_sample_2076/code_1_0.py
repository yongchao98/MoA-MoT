def solve_chess_puzzle():
    """
    This function determines and prints the pieces for the specified chess puzzle.

    The problem asks for the minimum number of additional pieces to create a
    Diagonal Corridor Mate with the Black King on h8 and the White King on a1.
    If multiple solutions have the same number of pieces, the one with the
    minimum total piece value is chosen.

    The optimal solution found is:
    - White: Bishop, Knight
    - Black: Pawn

    This setup involves a White Bishop on the a1-h8 diagonal (e.g., f6) to deliver
    the check. A Black Pawn on h7 provides the "corridor" element by blocking an
    escape square. A White Knight (e.g., on e6) efficiently controls the remaining
    two escape squares (g7 and g8). This uses a minimal 3 pieces with a total
    value of 3 (Bishop) + 3 (Knight) + 1 (Pawn) = 7.
    """
    
    # Piece values for sorting: Pawn=1, Knight=3, Bishop=3, Rook=5, Queen=9
    # The chosen pieces for the minimal solution.
    white_pieces = [
        ("White Bishop", 3), 
        ("White Knight", 3)
    ]
    
    black_pieces = [
        ("Black Pawn", 1)
    ]
    
    # Sort the pieces based on their value.
    # For pieces with the same value (Bishop and Knight), alphabetical order is a stable secondary sort.
    white_pieces.sort(key=lambda x: (x[1], x[0]))
    black_pieces.sort(key=lambda x: (x[1], x[0]))
    
    # Extract just the names for the final list
    sorted_white_piece_names = [p[0] for p in white_pieces]
    sorted_black_piece_names = [p[0] for p in black_pieces]
    
    # Combine the lists, with White's pieces first
    final_piece_list = sorted_white_piece_names + sorted_black_piece_names
    
    # Print the result as a comma-separated string
    print(", ".join(final_piece_list))

solve_chess_puzzle()
<<<White Bishop, White Knight, Black Pawn>>>