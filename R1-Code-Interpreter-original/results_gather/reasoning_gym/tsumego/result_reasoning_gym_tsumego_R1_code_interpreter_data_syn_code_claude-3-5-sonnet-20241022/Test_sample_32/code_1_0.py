def create_board():
    board = [['.'] * 12 for _ in range(12)]
    # Place black stones (X)
    black_stones = [
        (0,1),  # B12
        (9,3), (9,8),  # D6, I6
        (9,9), (9,10),  # J10, K10
        (8,8), (8,11),  # I9, L9
        (7,8),  # I8
        (6,9),  # J7
        (5,3), (5,7),  # D6, H6
        (3,10),  # K4
        (2,2)   # C2
    ]
    # Place white stones (O)
    white_stones = [
        (1,0), (1,10),  # A11, K11
        (2,8),  # I10
        (3,8), (3,9), (3,10),  # I9, J9, K9
        (4,3), (4,11),  # D8, L8
        (9,0),  # A3
        (10,3)  # D2
    ]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    return board

def is_valid_pos(x, y):
    return 0 <= x < 12 and 0 <= y < 12

def get_group_liberties(board, start_x, start_y, color, visited=None):
    if visited is None:
        visited = set()
    
    if not is_valid_pos(start_x, start_y) or board[start_x][start_y] != color:
        return set()
    
    if (start_x, start_y) in visited:
        return set()
    
    visited.add((start_x, start_y))
    liberties = set()
    
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = start_x + dx, start_y + dy
        if is_valid_pos(new_x, new_y):
            if board[new_x][new_y] == '.':
                liberties.add((new_x, new_y))
            elif board[new_x][new_y] == color:
                liberties.update(get_group_liberties(board, new_x, new_y, color, visited))
    
    return liberties

def count_captures(board, move_x, move_y):
    if not is_valid_pos(move_x, move_y) or board[move_x][move_y] != '.':
        return 0
    
    # Make the move
    board[move_x][move_y] = 'X'
    captures = 0
    
    # Check all adjacent positions for white groups with no liberties
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        x, y = move_x + dx, move_y + dy
        if is_valid_pos(x, y) and board[x][y] == 'O':
            liberties = get_group_liberties(board, x, y, 'O')
            if not liberties:
                captures += 1
    
    # Undo the move
    board[move_x][move_y] = '.'
    return captures

def find_best_move():
    board = create_board()
    best_move = None
    max_captures = 0
    
    for i in range(12):
        for j in range(12):
            if board[i][j] == '.':
                captures = count_captures(board, i, j)
                if captures > max_captures:
                    max_captures = captures
                    # Convert to game coordinates (A-L, 1-12)
                    col = chr(ord('A') + j)
                    row = 12 - i
                    best_move = f"{col}{row}"
    
    return best_move, max_captures

best_move, captures = find_best_move()
print(f"Best move: {best_move} (captures {captures} stones)")