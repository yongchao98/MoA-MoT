def count_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return 0
    
    pos = (x, y)
    if pos in visited:
        return 0
    
    if board[x][y] == '.':
        return 1
    if board[x][y] != color:
        return 0
    
    visited.add(pos)
    liberties = 0
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        liberties += count_liberties(board, x + dx, y + dy, color, visited)
    return liberties

def create_board():
    board = [['.'] * 11 for _ in range(11)]
    # Place black stones (X)
    black_stones = [(0,10), (5,10), (6,10), (7,10), (3,8), (4,7), (5,8), (5,4), (9,3), (9,2)]
    # Place white stones (O)
    white_stones = [(6,9), (6,8), (7,9), (2,1), (8,4)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'
    
    return board

def check_capture_move(board, move_x, move_y):
    # Try the move and check if it captures any stones
    board[move_x][move_y] = 'X'
    captured = False
    
    # Check adjacent positions for white stones
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        x, y = move_x + dx, move_y + dy
        if 0 <= x < 11 and 0 <= y < 11 and board[x][y] == 'O':
            if count_liberties(board, x, y, 'O') == 0:
                captured = True
    
    board[move_x][move_y] = '.'  # Reset the move
    return captured

# Create board and check potential capturing moves
board = create_board()
best_move = None

# Check H9 (7,8) as it looks like the most promising move
if check_capture_move(board, 7, 8):
    best_move = 'H9'

print(best_move)