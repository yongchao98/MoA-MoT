def solve_chess_puzzle():
    """
    This function returns the solution to the chess puzzle.
    The puzzle is a mate in 3 for black.
    The winning line is:
    1... Ra4+
    2. Kg5 Ra6
    3. If Kh6, R4a5#
    3. If Kg6, R6a5#
    3. If any rook move, R4a5#
    """
    num_moves = 3
    best_move = "Ra4+"
    print(f"{num_moves}. Best move: {best_move}")

solve_chess_puzzle()