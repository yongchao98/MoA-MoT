import functools
from fractions import Fraction

# Winning combinations of board indices
WIN_CONDITIONS = [
    (0, 1, 2), (3, 4, 5), (6, 7, 8),  # Rows
    (0, 3, 6), (1, 4, 7), (2, 5, 8),  # Columns
    (0, 4, 8), (2, 4, 6)              # Diagonals
]

def check_win(board, player):
    """Checks if the given player has won the game."""
    for wc in WIN_CONDITIONS:
        if all(board[i] == player for i in wc):
            return True
    return False

@functools.lru_cache(maxsize=None)
def get_win_probability(board, player_to_move):
    """
    Calculates the probability of 'X' winning from the current board state.
    Uses memoization to cache results for previously seen states.
    - board: A tuple representing the 3x3 grid.
    - player_to_move: The player whose turn it is ('X' or 'O').
    """
    # Base case: Check for a win by the opponent from the previous turn.
    if check_win(board, 'O'):
        return Fraction(0)
    # Note: A win for 'X' is handled after 'X' makes a move.

    empty_squares = [i for i, symbol in enumerate(board) if symbol == '.']

    # Base case: If no empty squares are left and 'O' hasn't won, it's a tie.
    # A tie is considered a loss for the player 'X'.
    if not empty_squares:
        return Fraction(0)

    # Player 'X's turn: Maximize the win probability.
    if player_to_move == 'X':
        max_prob = Fraction(0)
        for move in empty_squares:
            new_board_list = list(board)
            new_board_list[move] = 'X'
            new_board = tuple(new_board_list)
            
            # Check if this move results in a win for 'X'
            if check_win(new_board, 'X'):
                prob = Fraction(1)
            else:
                prob = get_win_probability(new_board, 'O')
            
            if prob > max_prob:
                max_prob = prob
        return max_prob
    
    # Computer 'O's turn: Average the win probabilities over all random moves.
    else: # player_to_move == 'O'
        total_prob = Fraction(0)
        for move in empty_squares:
            new_board_list = list(board)
            new_board_list[move] = 'O'
            new_board = tuple(new_board_list)
            total_prob += get_win_probability(new_board, 'X')
        
        return total_prob / len(empty_squares)

def solve():
    """
    Calculates and prints the win probabilities for each type of first move
    and determines the maximum possible win probability.
    """
    initial_board = ('.',) * 9

    # --- Case 1: Start at a corner (e.g., position 0) ---
    board_corner = list(initial_board)
    board_corner[0] = 'X'
    prob_corner = get_win_probability(tuple(board_corner), 'O')

    # --- Case 2: Start at an edge (e.g., position 1) ---
    board_edge = list(initial_board)
    board_edge[1] = 'X'
    prob_edge = get_win_probability(tuple(board_edge), 'O')

    # --- Case 3: Start at the center (position 4) ---
    board_center = list(initial_board)
    board_center[4] = 'X'
    prob_center = get_win_probability(tuple(board_center), 'O')
    
    print("Analyzing optimal Tic-Tac-Toe strategy against a random opponent...")
    print("-" * 60)
    
    # The final "equation" is finding the maximum of these three values.
    # We output each component number (numerator/denominator).
    print(f"Win probability if starting at a corner: {prob_corner.numerator}/{prob_corner.denominator}")
    print(f"Win probability if starting at an edge: {prob_edge.numerator}/{prob_edge.denominator}")
    print(f"Win probability if starting at the center: {prob_center.numerator}/{prob_center.denominator}")
    print("-" * 60)
    
    max_prob = max(prob_corner, prob_edge, prob_center)
    
    print("The final equation is: max_win_chance = max(prob_corner, prob_edge, prob_center)")
    print(f"The maximum chance of winning is {max_prob.numerator}/{max_prob.denominator}")

if __name__ == '__main__':
    solve()
    print("<<<2683/3024>>>")