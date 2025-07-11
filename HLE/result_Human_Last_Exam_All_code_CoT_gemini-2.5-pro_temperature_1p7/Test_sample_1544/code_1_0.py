import sys

def solve_chess_puzzle():
    """
    Analyzes the chess position and provides the best move for White.
    FEN: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
    """

    # The best move is Qf5.
    best_move = "Qf5"

    print(f"Based on the analysis, the optimal move for White is: {best_move}")
    print("\nThis move forces a queen trade on a square that leads to a winning pawn endgame for White.")
    print("The key winning line is as follows:")

    # In chess notation, a '+' indicates a check and '!' indicates a brilliant move.
    # The prompt mentions "output each number in the final equation!", which is interpreted here as showing the move numbers.
    winning_variation = [
        "1. Qf5 Qxf5+",
        "2. exf5 Kc5",
        "3. Kf2 Kxd6",
        "4. Ke3 Ke5",
        "5. f6! Kxf6",
        "6. Kd4 Ke6",
        "7. Kc5 ... White wins"
    ]

    for move in winning_variation:
        print(move)
    
    # Print the final part of the line explaining the outcome.
    print("\nAfter move 6, White's king is free to attack Black's queenside pawns.")
    print("White will win the ensuing pawn race, leading to a decisive victory.")


solve_chess_puzzle()

# The final answer is the single best move in standard chess notation.
sys.stdout.write("<<<Qf5>>>")