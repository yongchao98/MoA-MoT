def solve_chess_puzzle():
    """
    This function explains the solution to the given chess puzzle.
    It prints the step-by-step sequence for a forced checkmate.
    """
    print("White has a forced checkmate in 4 moves.")
    print("The winning sequence, a pattern known as Anastasia's Mate, starts with a Queen sacrifice.")
    print("\nHere is the main line of play:")
    
    # Move 1
    print("1. White: Qxh7+  |  Black: Kxh7")
    print("   White sacrifices the Queen to expose the Black King. Black is forced to capture.")

    # Move 2
    print("2. White: Ng6+   |  Black: Kh8")
    print("   The Knight checks the King and is protected by the h7-pawn. The King must move.")
    print("   (Note: If Black played 2...Kh6, White would play 3.Rh1+ Kg5 4.h4#, also a mate in 4).")
    
    # Move 3
    print("3. White: Rh1+   |  Black: Kg8")
    print("   The Rook on the now-open h-file checks the King, forcing it back to the back rank.")
    
    # Move 4
    print("4. White: Rh8#")
    print("   This is checkmate. The King is attacked by the Rook, which is protected by the Knight on g6.")
    print("   All escape squares are blocked or attacked.")

solve_chess_puzzle()