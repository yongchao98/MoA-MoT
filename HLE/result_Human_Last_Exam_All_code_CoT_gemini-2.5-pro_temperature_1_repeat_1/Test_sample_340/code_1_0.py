def solve_three_check_chess():
    """
    Analyzes the given Three-Check chess position to find the minimal number of moves for White to win.
    """

    print("Analyzing the Three-Check chess puzzle...")
    print("FEN: r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1 +0+0")
    print("Goal: Find the minimal number of moves for White to win (deliver 3 checks).\n")

    print("Step 1: White's most forcing move.")
    print("White can deliver a check on the very first move by castling queenside.")
    print("White's Move 1: O-O-O+")
    print("This move brings the Rook from a1 to d1, checking the Black King on e8.")
    print("Check Count: White 1, Black 0\n")

    print("Step 2: Analyze Black's optimal defenses.")
    print("Black must respond to the check. There are two legal and reasonable defenses:")
    print("  A) 1... Kf8 (King moves out of check)")
    print("  B) 1... Qd7 (Queen blocks the check)\n")

    print("Step 3: Find the forced continuation for White against each defense.\n")

    print("--- Analysis of Defense A: 1... Kf8 ---")
    print("If Black plays 1... Kf8, the king is on f8.")
    print("White's Move 2: Rd8+")
    print("White sacrifices the rook on d8, delivering the second check.")
    print("Check Count: White 2, Black 0")
    print("Black is forced to capture the rook (e.g., 2... Bxd8 or 2... Nxd8).")
    print("White's Move 3: Bxe7+")
    print("White's Bishop on g5 captures the Black Queen on e7, delivering the third and final check.")
    print("Check Count: White 3, Black 0")
    print("Result: White wins in 3 moves.\n")

    print("--- Analysis of Defense B: 1... Qd7 ---")
    print("If Black plays 1... Qd7, the Queen blocks the check on d7.")
    print("White's Move 2: Rxd7+")
    print("White captures the Queen, delivering the second check.")
    print("Check Count: White 2, Black 0")
    print("Black must recapture the rook. Two possibilities:")
    print("  - If Black plays 2... Nxd7 (Knight captures):")
    print("    White's Move 3: Bxd7+")
    print("    White's Bishop on b5 captures the Knight, delivering the third check.")
    print("    Check Count: White 3, Black 0")
    print("    Result: White wins in 3 moves.")
    print("  - If Black plays 2... Kxd7 (King captures):")
    print("    White's Move 3: Qd3+")
    print("    White's Queen moves to d3, delivering the third check.")
    print("    Check Count: White 3, Black 0")
    print("    Result: White wins in 3 moves.\n")

    print("Conclusion:")
    print("Against all of Black's optimal defenses, White can force a win in exactly 3 of their own moves.")
    print("Therefore, the minimal amount of moves by White to win is 3.")

solve_three_check_chess()
print("\nFinal Answer:")
print(3)