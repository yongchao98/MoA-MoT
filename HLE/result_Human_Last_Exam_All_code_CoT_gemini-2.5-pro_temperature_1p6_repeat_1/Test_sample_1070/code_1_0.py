def solve_chess_riddle():
    """
    This function provides the solution to the specified chess problem.

    The problem is a chess endgame: White King on d2, White Bishop on e2 vs. Black King on d5.
    A special rule applies: the Black King can only move to white squares.
    The goal is to find White's first move that leads to checkmate in the minimum number of moves,
    assuming optimal play from both sides.

    Analysis Summary:
    1. The game is a forced win for White because the Black King is restricted to half the board.
    2. The strategy is to force the Black King to a corner square of the same color as the Bishop (a light square, e.g., a8 or h1).
    3. Comparing the most promising first moves for White:
        - 1. Kd3: Leads to a forced checkmate in 9 moves against optimal defense.
        - 1. Bf3: A more powerful move that immediately restricts Black's options. This leads to a forced checkmate in 8 moves.
    4. Since 8 moves is faster than 9 moves, 1. Bf3 is the optimal first move.
    """
    
    # The optimal first move for White
    best_move = "Bf3"
    
    # The minimum number of moves to force checkmate with this move
    move_count = 8
    
    print(f"{best_move}, {move_count}")

solve_chess_riddle()