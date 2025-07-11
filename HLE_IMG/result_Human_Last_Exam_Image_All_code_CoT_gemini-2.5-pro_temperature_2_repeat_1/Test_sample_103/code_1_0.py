def solve_chess_puzzle():
    """
    Analyzes the chess position and determines the best move for White.

    The position analysis reveals that White has a powerful attack.
    The move Nh6+ is the most forcing and decisive option.

    1. Nh6+ ...
       - If Black plays ...gxh6, White follows up with 2. Qh5, creating
         a crushing and winning attack on the exposed king.
       - If Black plays ...Kh8, White follows up with 2. Qd7, threatening
         unstoppable checkmate on two squares (g7 and h7).

    This move is superior to other options like Nxd8+, which also wins but is less
    forcing, and Qxe4, which could lead to a material disadvantage for White.
    """
    best_move_analysis = {
        "move": "Nh6+",
        "evaluation": "Decisive advantage for White.",
        "reason": "Forces either a checkmating attack or a massive material gain."
    }

    print("The best move for White is I. Nh6+")
    print(f"Evaluation: {best_move_analysis['evaluation']}")
    print(f"Reason: {best_move_analysis['reason']}")

solve_chess_puzzle()
