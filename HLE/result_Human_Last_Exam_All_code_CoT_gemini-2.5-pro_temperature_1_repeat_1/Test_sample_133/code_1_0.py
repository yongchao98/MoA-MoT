import copy

def get_possible_moves(board):
    """Finds the next available spot in each column."""
    moves = []
    for col in range(7):
        for row in range(5, -1, -1):
            if board[row][col] == '.':
                moves.append((row, col))
                break
    return moves

def check_win(board, player, move):
    """Checks if a player has won with the given move."""
    r, c = move
    # Check horizontal
    for i in range(4):
        if c - i >= 0 and c - i + 3 < 7:
            if all(board[r][c - i + j] == player for j in range(4)):
                return True
    # Check vertical
    if r <= 2:
        if all(board[r + j][c] == player for j in range(4)):
            return True
    # Check diagonal /
    for i in range(4):
        if r + i <= 5 and c - i >= 0 and r + i - 3 >= 0 and c - i + 3 < 7:
            if all(board[r + i - j][c - i + j] == player for j in range(4)):
                return True
    # Check diagonal \
    for i in range(4):
        if r - i >= 0 and c - i >= 0 and r - i + 3 <= 5 and c - i + 3 < 7:
            if all(board[r - i + j][c - i + j] == player for j in range(4)):
                return True
    return False

def find_optimal_moves(board, player):
    """Finds moves that guarantee a win on the next turn."""
    optimal_moves = []
    possible_moves = get_possible_moves(board)

    for move in possible_moves:
        r, c = move
        temp_board = copy.deepcopy(board)
        temp_board[r][c] = player

        # After we move, find how many ways we can win on the *next* turn
        winning_spots = set()
        next_possible_moves = get_possible_moves(temp_board)
        for next_move in next_possible_moves:
            next_r, next_c = next_move
            final_board = copy.deepcopy(temp_board)
            final_board[next_r][next_c] = player
            if check_win(final_board, player, next_move):
                winning_spots.add(next_move)
        
        # If there are 2 or more independent winning spots, this is an optimal move
        if len(winning_spots) >= 2:
            optimal_moves.append(move)

    return optimal_moves

def format_moves(moves):
    """Converts (row, col) tuples to alphanumeric format."""
    # col_map = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g'}
    # Python's chr() can do this directly: a=97
    return [f"{chr(c + 97)}{6 - r}" for r, c in moves]

# Main execution
if __name__ == "__main__":
    # Board state from the image
    # Row 1 is board[0], Row 6 is board[5]
    # Col 'a' is board[...][0], Col 'g' is board[...][6]
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]
    player = 'O'
    
    # Find and print optimal moves
    winning_moves = find_optimal_moves(board, player)
    formatted_moves = format_moves(winning_moves)
    
    # Sort for consistent output
    formatted_moves.sort()
    
    print(", ".join(formatted_moves))