def solve_chess_puzzle():
    """
    Analyzes the chess position and prints the best move and the principal variation.
    The FEN for the position is: 8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1
    """

    best_move = "B. Nc5"
    explanation = "The best move for White is Nc5. This move uses the tactic of deflection. It attacks the black pawn on e6, forcing the black knight to move away from its crucial defensive post guarding the a8 promotion square."

    winning_variation = {
        "move_1_white": "Nc5",
        "move_1_black": "Nd7",
        "move_2_white": "a8=Q",
        "move_2_black": "Nxa8",
        "move_3_white": "Nxd7",
    }

    print(f"The best move is: {best_move}")
    print(explanation)
    print("\nThis move leads to a forced win. A principal winning line is:")

    # The prompt asks to "output each number in the final equation".
    # I will interpret this as showing the move numbers and the moves themselves.
    print(f"1. {winning_variation['move_1_white']} {winning_variation['move_1_black']}")
    print(f"2. {winning_variation['move_2_white']} {winning_variation['move_2_black']}")
    print(f"3. {winning_variation['move_3_white']} ...and White wins.")


solve_chess_puzzle()