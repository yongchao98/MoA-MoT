def solve_capablanca_puzzle():
    """
    This function explains the solution to the Capablanca chess puzzle
    and prints the minimum number of moves for White to win.
    """
    
    fen_string = "9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1"
    
    print("Analyzing the Capablanca chess puzzle:")
    print(f"FEN: {fen_string}")
    print("Board: 10x8")
    print("White to move and find the fastest win.")
    print("-" * 30)

    explanation = """
The position is a mate in 2 moves for White. Here is one decisive line:

1. White moves the Archbishop from h2 to c7 (Ac7+).
   - This is a Bishop move that puts the Black King on j8 in check.
   - The Archbishop on c7 also attacks the escape square i8 with its Knight-move ability.
   - The Black King cannot move. Black's only legal response is to block the check.

2. Black is forced to play Chancellor from f7 to d8 (...Cd8).
   - This is the only move that blocks the Bishop check from the Archbishop.

3. White plays Queen takes Chancellor on d8 (Qxd8#).
   - This is checkmate. The Black King is checked by the Queen.
   - The King still cannot move because the Archbishop on c7 covers the i8 square.
   - The Queen cannot be captured, and the check cannot be blocked.
"""
    
    print(explanation)
    
    minimal_moves = 2
    
    print("Thus, the minimal amount of moves by White to force a win is:")
    print(minimal_moves)

solve_capablanca_puzzle()