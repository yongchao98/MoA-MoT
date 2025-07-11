def solve():
    """
    This function explains the solution to the Capablanca chess puzzle.
    The puzzle asks for the minimum number of moves for White to force a checkmate.

    The solution is found through logical analysis of forced move sequences.

    1. White starts with the move Qd8+. This checks the Black King.
       Equation: 1. Qd8+

    2. Black has two main ways to respond to this check to prolong the game:
       a) ... Cf8 (Chancellor moves f7 -> f8)
       b) ... Ch8 (Chancellor moves f7 -> h8 as a knight)

    3. If Black plays 1... Cf8, White can deliver checkmate on the next move with 2. Ag4#.
       This is a mate in 2 moves for White.
       Equation: 1. Qd8+ Cf8 2. Ag4#

    4. However, if Black plays the optimal defense 1... Ch8, the game continues.
       White's next move is 2. Qe7+. This again checks the King.
       Black's only response is to block with 2... Cf7.
       Equation: 1. Qd8+ Ch8 2. Qe7+ Cf7

    5. Finally, White delivers checkmate with 3. Qxf7#.
       Equation: 1. Qd8+ Ch8 2. Qe7+ Cf7 3. Qxf7#

    Because Black can choose the line that extends the game to 3 moves for White,
    the minimal number of moves for White to *force* a win is 3.
    """
    
    # Final Answer
    # The longest path to a forced mate determines the number of moves.
    # The path is 1. Qd8+ Ch8 2. Qe7+ Cf7 3. Qxf7#
    # This involves 3 moves by White.
    
    number_of_moves = 3
    print(f"The minimal amount of moves by White to win is: {number_of_moves}")

solve()