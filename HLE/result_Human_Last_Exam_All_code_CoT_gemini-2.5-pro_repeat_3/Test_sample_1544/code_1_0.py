def solve_chess_puzzle():
    """
    This function determines the best move for White in the given chess position.
    The FEN for the position is: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
    """
    # Based on the analysis, the winning move is Qf8+.
    # 1. Qf8+ forces Black's king to c7 (Kc7).
    # 2. Qe7! creates a dual threat on the black queen and a mate threat on c7.
    # 3. Black's best defense involves a series of checks which White can parry.
    # 4. White has a subtle winning resource (6. Qe5+) that avoids a queen trade
    #    and leads to a winning endgame.
    best_move = "Qf8+"
    print(best_move)

solve_chess_puzzle()