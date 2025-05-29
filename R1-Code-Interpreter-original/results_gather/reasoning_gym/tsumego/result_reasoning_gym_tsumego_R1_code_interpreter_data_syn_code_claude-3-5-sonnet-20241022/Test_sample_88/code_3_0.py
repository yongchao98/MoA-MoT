def create_board():
    # Initialize empty board
    board = [['.'] * 13 for _ in range(13)]
    
    # Place White stones (O)
    white_positions = [
        (7, 12),  # H13
        (1, 10),  # B11
        (3, 9),   # D10
        (3, 6),   # D7
        (9, 6),   # J7
        (9, 5),   # J6
        (10, 5),  # K6
        (11, 5),  # L6
        (8, 4),   # I5
        (10, 4),  # K5
        (10, 3),  # K4
        (11, 2)   # L3
    ]
    
    # Place Black stones (X)
    black_positions = [
        (11, 9),  # L10
        (4, 8),   # E9
        (4, 7),   # E8
        (10, 6),  # K7
        (11, 6),  # L7
        (1, 5),   # B6
        (9, 4),   # J5
        (10, 3),  # K4
        (5, 1)    # F2
    ]
    
    # Place stones on board (adjusting for 0-based indexing)
    for col, row in white_positions:
        board[row][col] = 'O'
    
    for col, row in black_positions:
        board[row][col] = 'X'
    
    return board

def check_liberties(board, row, col, visited=None):
    if visited is None:
        visited = set()
    
    if row < 0 or row >= 13 or col < 0 or col >= 13:
        return False
    
    if board[row][col] == '.':
        return True
    
    if (row, col) in visited:
        return False
    
    visited.add((row, col))
    stone_type = board[row][col]
    
    # Check adjacent positions
    for dr, dc in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_row, new_col = row + dr, col + dc
        if 0 <= new_row < 13 and 0 <= new_col < 13:
            if board[new_row][new_col] == '.':
                return True
            if board[new_row][new_col] == stone_type:
                if check_liberties(board, new_row, new_col, visited):
                    return True
    
    return False

def count_captures(board, row, col):
    if board[row][col] != '.':
        return 0
    
    # Temporarily place black stone
    board[row][col] = 'X'
    captures = 0
    
    # Check for captures
    for dr, dc in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_row, new_col = row + dr, col + dc
        if 0 <= new_row < 13 and 0 <= new_col < 13:
            if board[new_row][new_col] == 'O':
                if not check_liberties(board, new_row, new_col):
                    captures += 1
    
    # Remove temporary stone
    board[row][col] = '.'
    return captures

def find_best_move():
    board = create_board()
    max_captures = 0
    best_move = None
    
    for row in range(13):
        for col in range(13):
            if board[row][col] == '.':
                captures = count_captures(board, row, col)
                if captures > max_captures:
                    max_captures = captures
                    best_move = (chr(col + ord('A')), row + 1)
    
    if best_move:
        print(f"<<<{best_move[0]}{best_move[1]}>>>")
    else:
        print("<<<No capturing move found>>>")

find_best_move()