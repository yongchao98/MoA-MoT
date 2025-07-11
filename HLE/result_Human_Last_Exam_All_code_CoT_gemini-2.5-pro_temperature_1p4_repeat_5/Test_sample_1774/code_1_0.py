def solve_chess_puzzle():
    """
    This function outlines the step-by-step solution to the chess puzzle
    and prints the number of moves to checkmate.
    """
    print("Analyzing the board, we can find a forced checkmate sequence for White.")
    print("The plan is to sacrifice the Queen to expose the Black King, then use the two Knights to deliver the final blow.")
    print("-" * 40)

    # Move 1: Queen Sacrifice
    move_1 = 1
    print(f"Move {move_1}: White plays Queen to H7, check (Qxh7+).")
    print("The Black King has no choice but to capture the White Queen (Kxh7).")
    print("The Black King is now exposed on H7.")
    print("-" * 40)

    # Move 2: Knight Check
    move_2 = 2
    print(f"Move {move_2}: White plays Knight from E5 to G6, check (Ng6+).")
    print("The Black King is forced to move back to its starting square (Kg8).")
    print("-" * 40)

    # Move 3: Checkmate
    move_3 = 3
    print(f"Move {move_3}: White plays the other Knight from E4 to G5, checkmate (Ng5#).")
    print("This is checkmate because:")
    print("  - The King on G8 is attacked by the Knight on G5.")
    print("  - The escape square H7 is covered by the White Bishop on D3.")
    print("  - The escape square H8 is covered by the White Knight on G6.")
    print("  - The check cannot be blocked or the checking piece captured.")
    print("-" * 40)
    
    total_moves = 3
    print(f"The total number of moves for White to checkmate the Black King is {total_moves}.")

solve_chess_puzzle()