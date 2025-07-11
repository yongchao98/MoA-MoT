def solve_chess_puzzle():
    """
    Analyzes the chess position and explains the best move for White.
    """
    best_move = "Nh6+"
    explanation = [
        "The best move for White is Nh6+.",
        "This move initiates a forced checkmate sequence.",
        "Let's break down the line:",
        "1. Nh6+ : White's knight delivers a check. The Black king cannot move to f8 (due to the Queen on d5) or h8 (due to the knight itself).",
        "1. ... gxh6 : This is Black's only legal move, forcing the capture of the knight and opening the h-file.",
        "2. Qh5+ : The White queen moves to h5, checking the king. The only legal move for Black is to move the king to g7.",
        "2. ... Kg7 : Black is forced to move the king.",
        "3. Qg4+ : The queen checks again from g4, forcing the king back to the h-file.",
        "3. ... Kh8 : The only legal move for the king.",
        "4. Qg8# : White delivers checkmate. The king is trapped and cannot escape the attack from the queen on g8.",
        f"Therefore, {best_move} is the winning move."
    ]

    for line in explanation:
        print(line)

solve_chess_puzzle()