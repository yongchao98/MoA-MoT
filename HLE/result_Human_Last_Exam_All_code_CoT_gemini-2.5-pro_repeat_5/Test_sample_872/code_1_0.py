from fractions import Fraction

# Memoization cache to store results for board states that have already been computed.
memo = {}

def get_winner(board):
    """
    Checks if there is a winner or a draw on the board.
    Returns 'X', 'O', 'Draw', or None if the game is ongoing.
    """
    lines = [
        (0, 1, 2), (3, 4, 5), (6, 7, 8),  # Rows
        (0, 3, 6), (1, 4, 7), (2, 5, 8),  # Columns
        (0, 4, 8), (2, 4, 6)             # Diagonals
    ]
    for a, b, c in lines:
        if board[a] == board[b] == board[c] and board[a] != '.':
            return board[a]
    if '.' not in board:
        return 'Draw'
    return None

def calculate_win_prob(board, player):
    """
    Recursively calculates the win probability for player 'X' from a given board state.
    - board: A tuple representing the 3x3 grid.
    - player: The current player to move ('X' or 'O').
    """
    # If this board state has been seen before, return the cached result.
    if board in memo:
        return memo[board]

    # Check for terminal states (win, loss, draw).
    winner = get_winner(board)
    if winner == 'X':
        return Fraction(1, 1)
    if winner == 'O' or winner == 'Draw':
        return Fraction(0, 1)

    empty_squares = [i for i, s in enumerate(board) if s == '.']

    if player == 'X':
        # My turn: I choose the move that maximizes my win probability.
        # The result is the maximum of the probabilities of all possible next states.
        best_prob = Fraction(0, 1)
        for move in empty_squares:
            new_board_list = list(board)
            new_board_list[move] = 'X'
            prob = calculate_win_prob(tuple(new_board_list), 'O')
            if prob > best_prob:
                best_prob = prob
        memo[board] = best_prob
        return best_prob
    else:  # player == 'O'
        # Computer's turn: It chooses a random move.
        # The result is the average of the probabilities of all possible next states.
        total_prob = Fraction(0, 1)
        for move in empty_squares:
            new_board_list = list(board)
            new_board_list[move] = 'O'
            total_prob += calculate_win_prob(tuple(new_board_list), 'X')
        
        num_moves = len(empty_squares)
        avg_prob = total_prob / num_moves if num_moves > 0 else Fraction(0,1)
        memo[board] = avg_prob
        return avg_prob

def solve():
    """
    Calculates the maximum win probability by trying all distinct first moves.
    """
    # Case 1: Start in the center (position 4)
    board_center_start = list('.........')
    board_center_start[4] = 'X'
    prob_center = calculate_win_prob(tuple(board_center_start), 'O')

    # Case 2: Start in a corner (position 0)
    # By symmetry, all corners are equivalent.
    board_corner_start = list('.........')
    board_corner_start[0] = 'X'
    prob_corner = calculate_win_prob(tuple(board_corner_start), 'O')

    # Case 3: Start on an edge (position 1)
    # By symmetry, all edges are equivalent.
    board_edge_start = list('.........')
    board_edge_start[1] = 'X'
    prob_edge = calculate_win_prob(tuple(board_edge_start), 'O')
    
    print("Calculating win probabilities for each opening move type...")
    print(f"1. Starting in the Center: {prob_center.numerator}/{prob_center.denominator}")
    print(f"2. Starting in a Corner: {prob_corner.numerator}/{prob_corner.denominator}")
    print(f"3. Starting on an Edge: {prob_edge.numerator}/{prob_edge.denominator}")
    
    # Find the maximum probability among the opening moves.
    max_prob = max(prob_center, prob_corner, prob_edge)
    
    print("\nThe maximum chance of winning you can give yourself is the maximum of these probabilities.")
    print(f"The final equation is: max({prob_center.numerator}/{prob_center.denominator}, {prob_corner.numerator}/{prob_corner.denominator}, {prob_edge.numerator}/{prob_edge.denominator}) = {max_prob.numerator}/{max_prob.denominator}")

solve()
<<<23/24>>>