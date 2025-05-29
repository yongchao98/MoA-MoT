def create_board():
    # Initialize empty board
    board = [['.'] * 13 for _ in range(13)]
    
    # Place White stones (O)
    white_stones = """
    H13 B11 D10 D7 J7 J6 K6 L6 I5 K5 K4 L3
    """
    
    # Place Black stones (X)
    black_stones = """
    L10 E9 E8 K7 L7 B6 X6 J5 K4 F2
    """
    
    # Helper function to convert coordinates to indices
    def coord_to_index(coord):
        col = ord(coord[0]) - ord('A')
        row = int(coord[1:]) - 1
        return row, col
    
    # Place stones on board
    for stone in white_stones.split():
        row, col = coord_to_index(stone)
        board[row][col] = 'O'
    
    for stone in black_stones.split():
        row, col = coord_to_index(stone)
        board[row][col] = 'X'
    
    return board

def has_liberty(board, row, col, visited=None):
    if visited is None:
        visited = set()
    
    if row < 0 or row >= 13 or col < 0 or col >= 13:
        return False
    
    if board[row][col] == '.':
        return True
    
    if (row, col) in visited:
        return False
    
    visited.add((row, col))
    stone = board[row][col]
    
    # Check all adjacent positions
    for dr, dc in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_row, new_col = row + dr, col + dc
        if 0 <= new_row < 13 and 0 <= new_col < 13:
            if board[new_row][new_col] == '.':
                return True
            if board[new_row][new_col] == stone and has_liberty(board, new_row, new_col, visited):
                return True
    
    return False

def would_capture(board, row, col):
    if board[row][col] != '.':
        return 0
    
    # Temporarily place black stone
    board[row][col] = 'X'
    captures = 0
    
    # Check for captures in adjacent positions
    for dr, dc in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_row, new_col = row + dr, col + dc
        if 0 <= new_row < 13 and 0 <= new_col < 13:
            if board[new_row][new_col] == 'O' and not has_liberty(board, new_row, new_col):
                captures += 1
    
    # Remove temporary stone
    board[row][col] = '.'
    return captures

def find_best_capture():
    board = create_board()
    max_captures = 0
    best_move = None
    
    for row in range(13):
        for col in range(13):
            if board[row][col] == '.':
                captures = would_capture(board, row, col)
                if captures > max_captures:
                    max_captures = captures
                    best_move = (chr(col + ord('A')), row + 1)
    
    if best_move:
        print(f"<<<{best_move[0]}{best_move[1]}>>>")
    else:
        print("<<<No capturing move found>>>")

find_best_capture()