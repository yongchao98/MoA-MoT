from fractions import Fraction

# Memoization cache to store results for previously seen states
memo = {}

def check_winner(board):
    """
    Checks the board for a winner or a tie.
    Returns:
        1: if X (player) wins
       -1: if O (computer) wins
        2: if the game is a tie
        0: if the game is still ongoing
    """
    win_conditions = [
        (0, 1, 2), (3, 4, 5), (6, 7, 8),  # Rows
        (0, 3, 6), (1, 4, 7), (2, 5, 8),  # Columns
        (0, 4, 8), (2, 4, 6)             # Diagonals
    ]
    for a, b, c in win_conditions:
        if board[a] == board[b] == board[c] and board[a] != 0:
            return board[a]
    if 0 not in board:
        return 2  # Tie
    return 0  # Game not over

def get_win_probability(board, player):
    """
    Recursively calculates the win probability for the current player.
    - board: a tuple representing the 3x3 grid (0=empty, 1=X, -1=O)
    - player: the current player to move (1 for X, -1 for O)
    """
    state = (board, player)
    if state in memo:
        return memo[state]

    winner = check_winner(board)
    if winner == 1:  # I (X) won
        return Fraction(1, 1)
    if winner == -1 or winner == 2:  # I lost or tied
        return Fraction(0, 1)

    possible_moves = [i for i, cell in enumerate(board) if cell == 0]

    if player == 1:  # My turn (X): I play optimally to maximize win probability
        best_prob = Fraction(-1, 1)
        for move in possible_moves:
            new_board_list = list(board)
            new_board_list[move] = 1
            prob = get_win_probability(tuple(new_board_list), -1)
            if prob > best_prob:
                best_prob = prob
        memo[state] = best_prob
        return best_prob
    else:  # Computer's turn (O): It plays randomly
        total_prob = Fraction(0, 1)
        num_moves = len(possible_moves)
        for move in possible_moves:
            new_board_list = list(board)
            new_board_list[move] = -1
            total_prob += get_win_probability(tuple(new_board_list), 1)
        
        avg_prob = total_prob / num_moves if num_moves > 0 else Fraction(0, 1)
        memo[state] = avg_prob
        return avg_prob

def solve():
    """
    Calculates the maximum win probability by trying all unique first moves.
    """
    initial_board = tuple([0] * 9)

    # Case 1: Start in the center (1 type of move)
    center_board = list(initial_board)
    center_board[4] = 1
    prob_center = get_win_probability(tuple(center_board), -1)
    
    # Case 2: Start in a corner (4 symmetric moves, we test one)
    corner_board = list(initial_board)
    corner_board[0] = 1
    prob_corner = get_win_probability(tuple(corner_board), -1)

    # Case 3: Start on an edge (4 symmetric moves, we test one)
    edge_board = list(initial_board)
    edge_board[1] = 1
    prob_edge = get_win_probability(tuple(edge_board), -1)

    # Find the maximum probability
    max_prob = max(prob_center, prob_corner, prob_edge)

    print(f"Win probability if starting in the center: {prob_center.numerator}/{prob_center.denominator}")
    print(f"Win probability if starting in a corner: {prob_corner.numerator}/{prob_corner.denominator}")
    print(f"Win probability if starting on an edge: {prob_edge.numerator}/{prob_edge.denominator}")
    print("\nTo maximize your chance of winning, you should choose the move with the highest probability.")
    print(f"The maximum chance of winning is the fraction: {max_prob.numerator}/{max_prob.denominator}")

solve()
<<<47/48>>>