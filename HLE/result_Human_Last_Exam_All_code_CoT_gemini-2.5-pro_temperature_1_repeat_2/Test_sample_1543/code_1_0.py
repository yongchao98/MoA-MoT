def solve():
    """
    This problem is a "mate in N" puzzle in Capablanca chess.
    The FEN is: 9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1

    The solution is a mate in 2 moves for White.

    1. White's key move: Qd2
    This move creates an unstoppable threat of Qj2#.

    2. Black's response:
    Black has no move to stop the coming checkmate. Any legal move by Black
    will be followed by White's mating move.

    3. White's mating move: Qj2#
    The Queen on j2 checks the King on j8.
    - The King cannot move to i8 (attacked by Queen).
    - The King cannot move to j7 (attacked by Queen).
    - The King cannot move to i7 (occupied by its own bishop).
    - The check cannot be blocked.
    - The Queen on j2 cannot be captured.

    This is a checkmate. The minimal number of moves for White to win is 2.
    """
    winning_moves_count = 2
    print(f"The minimal amount of moves by White to win is: {winning_moves_count}")

solve()