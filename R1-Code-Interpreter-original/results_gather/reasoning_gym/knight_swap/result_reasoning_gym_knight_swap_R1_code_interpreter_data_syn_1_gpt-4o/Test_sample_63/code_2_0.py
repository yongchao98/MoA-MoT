def is_valid_move(start, end, board):
    # Check if the move is within bounds and the end position is empty
    if end in board and board[end] == '.':
        return True
    return False

def get_knight_moves(position):
    # Calculate all possible knight moves from a given position
    x, y = position
    moves = [
        (x + 2, y + 1), (x + 2, y - 1),
        (x - 2, y + 1), (x - 2, y - 1),
        (x + 1, y + 2), (x + 1, y - 2),
        (x - 1, y + 2), (x - 1, y - 2)
    ]
    return moves

def swap_knights():
    # Initial board setup
    board = {
        (1, 3): '.', (2, 3): 'B', (3, 3): '.', (4, 3): 'B',
        (1, 2): '.', (2, 2): 'w', (3, 2): '.', (4, 2): 'w',
        (1, 1): '.', (2, 1): '.', (3, 1): '.', (4, 1): '.'
    }
    
    # Target positions
    target = {
        (2, 3): 'w', (4, 3): 'w',
        (2, 2): 'B', (4, 2): 'B'
    }
    
    # Convert board positions to coordinates
    pos_to_coord = {'A': 1, 'B': 2, 'C': 3, 'D': 4}
    coord_to_pos = {v: k for k, v in pos_to_coord.items()}
    
    # Initial positions of knights
    white_knights = [(2, 2), (4, 2)]
    black_knights = [(2, 3), (4, 3)]
    
    # Moves list
    moves = []
    
    def dfs(w_knights, b_knights, current_moves, current_board, depth):
        # Check if the target configuration is reached
        if all(current_board[pos] == target[pos] for pos in target):
            return current_moves
        
        if depth > 10:  # Limit the depth to prevent excessive computation
            return None
        
        # White's turn
        for i, wk in enumerate(w_knights):
            possible_moves = get_knight_moves(wk)
            for move in possible_moves:
                if is_valid_move(wk, move, current_board):
                    # Make the move
                    new_board = current_board.copy()
                    new_board[move] = 'w'
                    new_board[wk] = '.'
                    new_w_knights = w_knights[:i] + [move] + w_knights[i+1:]
                    new_moves = current_moves + [f"w,{coord_to_pos[wk[0]]}{wk[1]},{coord_to_pos[move[0]]}{move[1]}"]
                    result = dfs(new_w_knights, b_knights, new_moves, new_board, depth + 1)
                    if result:
                        return result
        
        # Black's turn
        for i, bk in enumerate(b_knights):
            possible_moves = get_knight_moves(bk)
            for move in possible_moves:
                if is_valid_move(bk, move, current_board):
                    # Make the move
                    new_board = current_board.copy()
                    new_board[move] = 'B'
                    new_board[bk] = '.'
                    new_b_knights = b_knights[:i] + [move] + b_knights[i+1:]
                    new_moves = current_moves + [f"B,{coord_to_pos[bk[0]]}{bk[1]},{coord_to_pos[move[0]]}{move[1]}"]
                    result = dfs(w_knights, new_b_knights, new_moves, new_board, depth + 1)
                    if result:
                        return result
        
        return None
    
    result = dfs(white_knights, black_knights, moves, board, 0)
    return result if result else "No"

result = swap_knights()
print(result)