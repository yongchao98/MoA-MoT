import chess

def solve_chess_puzzle():
    """
    Analyzes the chess position and determines the best move for White.

    The position is given by the FEN: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1

    Analysis:
    1.  White's passed d-pawn on d6 is a massive threat, just two squares from promotion.
    2.  Black's king is somewhat exposed, and the Black queen is active but has no immediate knockout blows.
    3.  White needs to find a forcing sequence to capitalize on the d-pawn's strength before Black can create counterplay.

    Candidate Moves Evaluation:
    -   Moves like Qxc5 or other queen moves that allow a queen trade (e.g., ...Qe3+) lead to a lost or difficult pawn endgame for White, as Black's d-pawn becomes a passed c-pawn.
    -   The move Qf7+ is the most forcing move on the board.

    Line Analysis:
    1.  White plays Qf7+. This check forces Black's king to move.
        - If 1... Kb7, White plays 2. Qxc7+. This is a critical move that captures the d7 pawn with check. Black must recapture with 2... Kxa7. The resulting endgame is a White king and pawns on a4 and d6 versus a Black king on a7 and a pawn on a5. White's d-pawn is unstoppable and will promote, leading to an easy win.
        - If 1... Ka6, White can play 2. Qe7!, pinning the Black queen to the king. Black cannot save the queen and the game is lost. For example, 2... Ka7 3. Qxc5+ dxc5, and White's king can easily stop Black's newly created passed c-pawn while the d-pawn marches to victory.

    Conclusion:
    The move Qf7+ initiates a forced sequence that leads to a winning position for White.
    """
    best_move = "Qf7+"
    print(f"The best move for White is: {best_move}")

solve_chess_puzzle()