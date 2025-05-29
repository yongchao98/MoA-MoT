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
    
    # Simulate moves
    for _ in range(10):  # Limit to 10 moves for simplicity
        # White's turn
        for wk in white_knights:
            possible_moves = get_knight_moves(wk)
            for move in possible_moves:
                if is_valid_move(wk, move, board):
                    # Make the move
                    board[move] = 'w'
                    board[wk] = '.'
                    moves.append(f"w,{coord_to_pos[wk[0]]}{wk[1]},{coord_to_pos[move[0]]}{move[1]}")
                    white_knights.remove(wk)
                    white_knights.append(move)
                    break
        
        # Black's turn
        for bk in black_knights:
            possible_moves = get_knight_moves(bk)
            for move in possible_moves:
                if is_valid_move(bk, move, board):
                    # Make the move
                    board[move] = 'B'
                    board[bk] = '.'
                    moves.append(f"B,{coord_to_pos[bk[0]]}{bk[1]},{coord_to_pos[move[0]]}{move[1]}")
                    black_knights.remove(bk)
                    black_knights.append(move)
                    break
        
        # Check if the target configuration is reached
        if all(board[pos] == target[pos] for pos in target):
            return moves
    
    return "No"

result = swap_knights()
print(result)