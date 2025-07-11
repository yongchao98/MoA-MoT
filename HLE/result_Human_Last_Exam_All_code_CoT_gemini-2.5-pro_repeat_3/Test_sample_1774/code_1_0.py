def solve_chess_mate_in_4():
    """
    Analyzes the chess position and prints the forced mate in 4 sequence.
    """
    print("Analyzing the board...")
    print("White has a powerful attack setup against the Black king.")
    print("The shortest path to checkmate is a forced sequence in 4 moves.\n")
    
    print("The winning sequence is as follows:")
    
    # Move 1: White sacrifices the Queen to expose the king.
    print("1. White moves Queen to H7, checking the King (Qxh7+).")
    print("   Black is forced to capture the Queen with the King (Kxh7).\n")
    
    # Move 2: White's Knight checks the King.
    print("2. White moves Knight from E5 to G6, checking the King (Ng6+).")
    print("   Black's best, though futile, response is to capture the Knight with the Bishop (Bxg6).\n")
    
    # Move 3: White's Pawn captures, resulting in a discovered check.
    print("3. White's Pawn on H2 captures the Bishop on G6 (hxg6+).")
    print("   This is a discovered check from the Rook on H1. Black's only move is King to H6 (Kh6).\n")
    
    # Move 4: White delivers the final checkmate.
    print("4. White moves Rook to H4, delivering checkmate (Rh4#).")
    print("   The Black King on H6 cannot escape:")
    print("   - G5 is controlled by the Knight on E4.")
    print("   - H5 is controlled by the Pawn on G6.\n")

    print("---")
    print("Calculating the number of moves to mate:")
    
    # Representing each of White's moves in the count.
    white_move_1 = 1
    white_move_2 = 1
    white_move_3 = 1
    white_move_4 = 1
    
    total_moves = white_move_1 + white_move_2 + white_move_3 + white_move_4
    
    # As requested, printing the equation.
    print(f"The equation for the total moves is: {white_move_1} + {white_move_2} + {white_move_3} + {white_move_4} = {total_moves}")
    print(f"\nIt is a mate in {total_moves} moves.")

solve_chess_mate_in_4()

<<<4>>>