def solve_chess_endgame_puzzle():
    """
    This script provides the solution to the specified chess endgame problem.

    Problem State:
    - White: King on d2, Bishop on e2.
    - Black: King on d5.
    - Special Rule: Black's king is cursed to only move on white squares.
    - Goal: Find White's best first move to achieve the fastest checkmate.

    Analysis Summary:
    The optimal move for White is 1. Kd3. This move restricts the Black King's movement by controlling the c4 and e4 squares. Black is forced to move to either c6 or e6.

    - If Black plays 1... Ke6, White can force a checkmate in 3 more moves, for a total of 4 moves (e.g., 1. Kd3 Ke6 2. Kd4 Kf7 3. Ke5 Ke8 4. Bf3#). Black can get mated in 3 moves if they play sub-optimally (2... Kf5), but they will choose the line that lasts 4 moves.
    - If Black plays 1... Kc6, a symmetrical line of play also leads to a forced checkmate in 4 moves (e.g., 1. Kd3 Kc6 2. Kc4 Kb7 3. Kb5 Ka8 4. Ba6#).

    Since both of Black's responses lead to a mate in 4 moves, this is the guaranteed outcome. Other White moves, like 1. Bc4+, allow Black to survive longer.
    """
    # The optimal first move for white.
    best_move = "Kd3"

    # The minimum number of moves to a guaranteed checkmate.
    number_of_moves = 4

    # Print the answer in the required "move, number" format.
    # The final equation is the combination of the move and the number of turns.
    print(f"{best_move}, {number_of_moves}")

solve_chess_endgame_puzzle()