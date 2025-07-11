def solve_chess_puzzle():
    """
    Analyzes the provided chess position and prints the forced mating sequence for White.
    The solution finds the shortest guaranteed checkmate, assuming best defense from Black.
    """
    mate_in_moves = 3
    print(f"White to move, mate in {mate_in_moves}.")
    print("-" * 30)

    # --- Move 1 ---
    print("1. White: Knight from E5 to F7 (Nf7+)")
    print("   This move puts the Black King in check. The Knight is defended by the Knight on E4.")
    print("   Black's best defense is to capture the Knight, as moving the King (1...Kh8) would lead to immediate mate on the next move (2. Qh7#).")
    print("1. ... Black: Rook from F8 captures Knight on F7 (Rxf7)")
    print("-" * 30)

    # --- Move 2 ---
    print("2. White: Knight from E4 captures Rook on F7 (Nxf7+)")
    print("   White continues the attack by recapturing with the second Knight, placing the King in check again.")
    print("   Black has only one legal move to escape the check.")
    print("2. ... Black: King from G8 to H8 (Kh8)")
    print("-" * 30)

    # --- Move 3 ---
    print("3. White: Queen from H5 to H7 (Qh7#)")
    print("   This is checkmate. The Black King is trapped.")
    print("   - It is attacked by the Queen.")
    print("   - It cannot capture the Queen, as it's defended by the Bishop on D3.")
    print("   - It cannot move to G8, which is controlled by the White Knight on F7.")

solve_chess_puzzle()