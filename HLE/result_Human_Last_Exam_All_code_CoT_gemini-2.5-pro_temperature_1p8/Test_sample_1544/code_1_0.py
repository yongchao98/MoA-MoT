def solve_chess_puzzle():
    """
    Analyzes the chess position from FEN 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
    and prints the best move for White with a clear explanation.
    """
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    best_move = "Qxc5+"
    explanation = (
        "The optimal move for White is to force a Queen exchange, "
        "transitioning to a winning King and Pawn endgame."
    )
    winning_line = [
        ("1.", "Qxc5+", "bxc5+"),
        ("2.", "Kf1", "Kc7"),
        ("3.", "Ke2", "Kd7"),
        ("4.", "Kd3", "...")
    ]

    print(f"Analyzing the chess position: {fen}")
    print("\nWhite's best move is a capture that simplifies the position.")
    print(f"\nThe best move is: {best_move}")
    print(f"\nRationale: {explanation}")
    print("\nThe game follows a forced sequence leading to a winning endgame for White:")
    
    # "Outputting each number in the final equation" here means showing the move numbers.
    for move_num, white_move, black_move in winning_line:
        if white_move and black_move:
            print(f"{move_num} {white_move} {black_move}")
        elif white_move:
             print(f"{move_num} {white_move}")
    
    print("\nAfter these moves, the Black King is tied to defending against the d6-pawn,")
    print("while the White King is free to attack and win Black's other pawns.")

solve_chess_puzzle()
print("\n<<<Qxc5+>>>")