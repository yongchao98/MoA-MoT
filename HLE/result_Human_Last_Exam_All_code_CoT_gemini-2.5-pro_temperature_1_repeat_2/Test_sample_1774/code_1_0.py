def solve_chess_mate_in_4():
    """
    This function prints the step-by-step solution to the provided chess puzzle,
    which is a forced checkmate in 4 moves for White.
    """
    
    moves_to_mate = 4
    
    # The full algebraic notation for the winning sequence
    # Each tuple contains the move number, White's move, and Black's move.
    sequence = [
        (1, "Ng5", "hxg5"),
        (2, "Qh7+", "Kf8"),
        (3, "Nxf7", "Kxf7"),
        (4, "Bg6#", "")
    ]

    print(f"White has a forced checkmate in {moves_to_mate} moves.")
    print("Here is the winning sequence:\n")

    # Print the moves in the sequence
    for move in sequence:
        move_number = move[0]
        white_move = move[1]
        black_move = move[2]
        
        # The final move has no response from Black
        if black_move:
            print(f"Move {move_number}: White plays {white_move}, Black must play {black_move}")
        else:
            print(f"Move {move_number}: White plays {white_move}, delivering checkmate.")

    print("\n--- Final Equation (Standard Notation) ---")
    for move in sequence:
        move_number = move[0]
        white_move = move[1]
        black_move = move[2]
        print(f"{move_number}. {white_move} {black_move}")

solve_chess_mate_in_4()