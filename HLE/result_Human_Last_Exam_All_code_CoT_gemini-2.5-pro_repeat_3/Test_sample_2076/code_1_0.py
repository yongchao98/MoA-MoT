def solve_chess_puzzle():
    """
    Determines the minimum pieces for a Diagonal Corridor Mate under given constraints.

    This function follows the logic outlined:
    1. The mating piece for White must be a Bishop to minimize value.
    2. The Black King's escape square (g7) must be blocked by a Black piece.
    3. To ensure Black has a material advantage (as stipulated), the Black piece's
       value must be greater than the White Bishop's (3).
    4. The lowest-value Black piece greater than 3 is a Rook (5).
    5. This results in the minimal set: one White Bishop and one Black Rook.
    """
    
    # Define piece values for clarity in logic, though not strictly needed for the final string
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9
    }

    # Step 1: Determine the minimal white checking piece
    white_piece = "White Bishop"
    white_value = piece_values["Bishop"]

    # Step 2: Determine the minimal black blocking piece that gives Black material advantage
    black_piece = "Black Rook"
    black_value = piece_values["Rook"]

    # Step 3: Assemble the list of pieces as per the output format
    # White pieces list, sorted by value (only one piece)
    white_pieces_list = [white_piece]
    
    # Black pieces list, sorted by value (only one piece)
    black_pieces_list = [black_piece]

    # Combine the lists
    final_pieces = white_pieces_list + black_pieces_list
    
    # Format as a comma-separated string
    result_string = ", ".join(final_pieces)
    
    print(result_string)

solve_chess_puzzle()
<<<White Bishop, Black Rook>>>