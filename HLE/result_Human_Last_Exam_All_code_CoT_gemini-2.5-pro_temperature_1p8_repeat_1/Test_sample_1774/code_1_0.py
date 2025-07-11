def solve_chess_puzzle():
    """
    This function explains the step-by-step solution to the chess puzzle.
    It details the forced checkmate sequence for White.
    """
    print("The puzzle is a mate in 3 moves.")
    print("Here is the sequence of moves that forces checkmate against the best defense:\n")

    # Move 1
    print("--- Move 1 ---")
    print("White's Move: Knight from E5 to F7 (Nf7+)")
    print("This puts the Black King on G8 in check. Black must respond.")
    print("Black's best defense is to capture the knight, prolonging the game. Black plays King from G8 captures F7 (Kxf7).\n")

    # Move 2
    print("--- Move 2 ---")
    print("White's Move: Knight from E4 to G5 (Ng5+)")
    print("This is a second check. The Knight on G5 now covers the King's escape squares F7 and H7.")
    print("Black's only legal move is to move the King back to G8. Black plays King from F7 to G8 (Kg8).\n")

    # Move 3
    print("--- Move 3 ---")
    print("White's Move: Queen from H5 to H7 (Qh7#)")
    print("This is checkmate. The Black King on G8 is attacked by the Queen and cannot escape:")
    print("- It cannot move to F7 or H7, which are attacked by the Queen.")
    print("- It cannot move to F8, which is blocked by its own Rook.")
    print("- It cannot capture the Queen on H7, as the square is also defended by the Knight on G5.")

if __name__ == "__main__":
    solve_chess_puzzle()
<<<3>>>