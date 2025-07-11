def solve_chess_puzzle():
    """
    This function provides the solution to the specified chess endgame puzzle.

    The analysis concludes that the best first move for White is one that maximally
    restricts the Black King's movement, forcing it into a mating net.
    The move 1. Bd3 cuts off two of the Black King's four initial escape squares.
    Assuming optimal play from both sides (where Black tries to delay mate as
    long as possible), this move leads to a forced checkmate in 10 moves.
    """
    # The optimal first move for White.
    best_move = "Bd3"

    # The number of moves to achieve checkmate with optimal play from both sides.
    moves_to_mate = 10

    # Print the result in the specified "move, number" format.
    print(f"{best_move}, {moves_to_mate}")

solve_chess_puzzle()