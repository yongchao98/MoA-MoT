def solve_chess_puzzle():
    """
    Analyzes the chess position and determines the best move for white.
    """
    position_description = """
White's position is dangerous due to the Black Queen on h4 and Bishop on e4, which threaten the White King.
However, White has a winning move that neutralizes the threat and secures a material advantage.
Let's analyze the best option: Ne5+
"""

    move_analysis = """
1. White plays Ne5+. The Knight moves from f7 to e5, delivering a check to the Black King on g8.
   - Knight at f7 -> e5
   - King at g8 is in check.

2. Black must respond to the check. The only safe square for the king is h8.
   - Black plays Kh8.

3. Now, White's Knight on e5 performs a critical defensive role: it protects the g2 pawn.
   This prevents Black's main counter-attacking idea of ...Qxg2+.

4. With the king-side secured, White can now play Qxe4.
   - White Queen at d5 captures Black Bishop at e4.

5. Result: White wins a bishop and completely stops Black's attack. This leads to a decisive advantage for White.
Other moves like Qxe4 directly are good but allow Black counterplay with ...Qxg2+. Ne5+ is the most precise and winning move.
"""

    print("Step-by-step analysis of the best move:")
    print(position_description)
    print("The best move is Ne5+.")
    print(move_analysis)


solve_chess_puzzle()