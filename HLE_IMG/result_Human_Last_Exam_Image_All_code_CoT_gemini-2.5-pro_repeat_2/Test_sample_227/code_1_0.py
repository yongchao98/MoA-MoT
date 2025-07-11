def solve_chess_puzzle():
    """
    This function determines the solution to the given chess puzzle.
    """
    # The puzzle is a mate in 2.
    number_of_moves = 2
    # The best move for black is R3a5.
    # This move puts white in zugzwang. The white king, which previously could only move to g5, is now completely trapped.
    # 1... R3a5
    # - If white plays 2. Rd4, then 2... Rxd4 is checkmate.
    # - If white plays any other move (e.g., 2. a3), then 2... Rd3 is checkmate.
    # Thus, black forces a mate in 2 moves.
    best_move = "R3a5"

    print(f"{number_of_moves}. Best move: {best_move}")

solve_chess_puzzle()