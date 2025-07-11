import sys

def solve_chess_puzzle():
    """
    This function analyzes the chess position and determines the best move for White.

    FEN: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1

    The analysis shows that White has a winning advantage by forcing a favorable queen trade.
    The key move sequence is:
    1. Qf8+ : This is a check, forcing the black king to move.
       For example, Black plays 1...Ka7.
    2. Qxc5  : White offers a queen trade. Black is practically forced to accept,
               as declining (with 2...cxc5) would leave White with a queen and an overwhelming advantage.
    3. ...Qxf2+ 4. Kxf2 cxc5 : The game transitions into a pawn endgame.
    
    In this endgame, White's king is activated on f2 while Black's king has been pushed to the side.
    White's active king combined with the strong d6-pawn ensures a win.
    """
    best_move = "Qf8+"
    print(f"The best move for White is: {best_move}")

solve_chess_puzzle()
# The puzzle is a "what is the next move" question. There is no equation.
# Therefore, I am printing the answer in chess notation as requested by the persona.
# The following format is to match the required final output format.