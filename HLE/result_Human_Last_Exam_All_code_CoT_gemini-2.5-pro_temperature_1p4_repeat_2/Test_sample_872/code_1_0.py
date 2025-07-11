import collections
from fractions import Fraction
from functools import lru_cache

# The winning combinations of board positions
WIN_COMBINATIONS = [
    [0, 1, 2], [3, 4, 5], [6, 7, 8],  # Rows
    [0, 3, 6], [1, 4, 7], [2, 5, 8],  # Columns
    [0, 4, 8], [2, 4, 6]              # Diagonals
]

def check_win(board, player):
    """Checks if the given player has won the game."""
    for combo in WIN_COMBINATIONS:
        if all(board[pos] == player for pos in combo):
            return True
    return False

@lru_cache(maxsize=None)
def get_win_prob(board, is_my_turn):
    """
    Recursively calculates the maximum win probability from a given board state.
    - board: A tuple representing the 3x3 grid.
    - is_my_turn: A boolean, True if it's the player's turn, False for computer.
    """
    # Base case: check for a terminal state (win, loss, or tie)
    if check_win(board, 'X'):  # I win
        return Fraction(1)
    if check_win(board, 'O'):  # I lose
        return Fraction(0)
    
    open_squares = [i for i, mark in enumerate(board) if mark == '.']
    if not open_squares:  # Tie game
        return Fraction(0)

    # Recursive step
    if is_my_turn:
        # I choose the move that maximizes my chance of winning.
        return max(
            get_win_prob(board[:move] + ('X',) + board[move+1:], False)
            for move in open_squares
        )
    else:  # Computer's turn
        # The computer plays a random move. Probability is the average over all possibilities.
        return sum(
            get_win_prob(board[:move] + ('O',) + board[move+1:], True)
            for move in open_squares
        ) / len(open_squares)

def solve():
    """
    Calculates the win probability for the three unique starting moves
    (corner, edge, center) and finds the maximum.
    """
    print("Calculating the maximum win probability...")
    
    # 1. My first move is a corner (position 0)
    board_corner_start = list('.........')
    board_corner_start[0] = 'X'
    prob_corner = get_win_prob(tuple(board_corner_start), is_my_turn=False)

    # 2. My first move is an edge (position 1)
    board_edge_start = list('.........')
    board_edge_start[1] = 'X'
    prob_edge = get_win_prob(tuple(board_edge_start), is_my_turn=False)
    
    # 3. My first move is the center (position 4)
    board_center_start = list('.........')
    board_center_start[4] = 'X'
    prob_center = get_win_prob(tuple(board_center_start), is_my_turn=False)

    # The optimal strategy is to choose the start position with the max win probability.
    max_prob = max(prob_corner, prob_edge, prob_center)

    print("\nWin probability if starting in a corner: {}/{}".format(prob_corner.numerator, prob_corner.denominator))
    print("Win probability if starting on an edge: {}/{}".format(prob_edge.numerator, prob_edge.denominator))
    print("Win probability if starting in the center: {}/{}".format(prob_center.numerator, prob_center.denominator))
    
    print(f"\nThe maximum chance of winning is the maximum of these probabilities.")
    print(f"max({prob_corner.numerator}/{prob_corner.denominator}, {prob_edge.numerator}/{prob_edge.denominator}, {prob_center.numerator}/{prob_center.denominator}) = {max_prob.numerator}/{max_prob.denominator}")

solve()
<<<149/210>>>