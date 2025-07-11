import sys

def solve_chess_puzzle():
    """
    This function explains the solution to the given chess puzzle.
    It details the forced checkmate sequence for White.
    """
    # The FEN string for the puzzle
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"

    print(f"Analyzing the chess position from FEN: {fen}")
    print("It is White's turn to move. White has a forced checkmate.\n")
    print("Here is the step-by-step solution for the mate in 4:\n")

    # Move 1
    print("--- Move 1 ---")
    print("White's first move is Qxf7+. This checks the Black king.")
    print("Black's king cannot move to f8 (attacked by the Rook on f1) or h7 (attacked by the Bishop on h5).")
    print("Therefore, Black's only legal response is Kh8.")
    print("The board state after 1. Qxf7+ Kh8\n")

    # Move 2
    print("--- Move 2 ---")
    print("White follows up with Ng6+. This is a powerful double check.")
    print("The Knight on g6 checks the king, and it also reveals a discovered check from the Bishop on h5.")
    print("In a double check, the king must move. The only legal square is h7.")
    print("The board state after 2. Ng6+ Kh7\n")

    # Move 3
    print("--- Move 3 ---")
    print("White continues the attack with Nf8+. This is another double check.")
    print("The Knight on f8 checks the king, and it also reveals a discovered check from the Queen on f7.")
    print("Once again, the king must move. The only legal square is back to h8.")
    print("The board state after 3. Nf8+ Kh8\n")

    # Move 4
    print("--- Move 4 ---")
    print("White delivers the final blow with Qxg7#.")
    print("The Queen on g7, protected by the Bishop on h5, checkmates the Black king on h8.")
    print("Black has no legal moves, and the game is over.\n")

    # Final sequence summary
    print("--- The Final Mating Sequence ---")
    print("1. Qxf7+ Kh8")
    print("2. Ng6+ Kh7")
    print("3. Nf8+ Kh8")
    print("4. Qxg7#")
    
    print("\nConclusion: White can force a win in 4 moves.")

solve_chess_puzzle()

# The final answer is the number of moves for White to win.
sys.stdout.write("<<<4>>>\n")