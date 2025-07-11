def solve_chess_puzzle():
    """
    This function explains and prints the solution to the chess puzzle.
    """
    moves_to_mate = 2
    
    print("The plan is to use a queen sacrifice to expose the black king,")
    print("followed by a knight move to deliver checkmate.")
    print("-" * 40)
    
    # Move 1
    white_move_number_1 = 1
    black_move_number_1 = 1
    print(f"Move {white_move_number_1} (White): Queen takes pawn on H7, check (Qxh7+).")
    print(f"Move {black_move_number_1} (Black): King takes Queen on H7 (Kxh7).")
    
    # Move 2
    white_move_number_2 = 2
    print(f"Move {white_move_number_2} (White): Knight moves to G6, checkmate (Ng6#).")
    
    print("-" * 40)
    print("This is checkmate because the Black King has no legal moves.")
    
    print(f"\nFinal result: It takes White {moves_to_mate} moves to deliver checkmate.")

solve_chess_puzzle()
<<<2>>>