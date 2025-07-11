import sys

def solve_chess_puzzle():
    """
    This function prints the step-by-step solution to the chess puzzle.
    The correct sequence leads to a "smothered mate".
    """
    print("The correct sequence of moves that forces checkmate is:")
    
    # Move 1 for Black
    move_number_1 = "1."
    black_move_1 = "... Qg1+"
    print(f"{move_number_1} {black_move_1}")
    print("Explanation: Black sacrifices the Queen by moving it to g1. This is a check, and White's only legal response is to capture the Queen with the Rook on f1.")

    # Move 2 for White and Black
    move_number_2 = "2."
    white_move_2 = "Rxg1"
    black_move_2 = "Nf2#"
    print(f"{move_number_2} {white_move_2} {black_move_2}")
    print("Explanation: White's Rook captures the Queen. Black's Knight then moves to f2, delivering a checkmate. The White King on h1 is trapped by its own Rook and pawn, a pattern known as a smothered mate.")

solve_chess_puzzle()