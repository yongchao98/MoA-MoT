def is_inside_board(x, y):
    return 0 <= x < 13 and 0 <= y < 13

def get_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if not is_inside_board(x, y) or board[y][x] != color:
        return 0
    
    pos = (x, y)
    if pos in visited:
        return 0
    
    visited.add(pos)
    liberties = 0
    
    # Check adjacent points
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        if is_inside_board(new_x, new_y):
            if board[new_y][new_x] == '.':
                liberties += 1
            elif board[new_y][new_x] == color:
                liberties += get_liberties(board, new_x, new_y, color, visited)
    
    return liberties

# Create the board
board = [['.'] * 13 for _ in range(13)]

# Set up the relevant stones
# Black stones (X)
black_stones = [(1, 12), (3, 3), (2, 4), (0, 4), (2, 2), (5, 3)]
# White stones (O)
white_stones = [(1, 1), (1, 2), (1, 3), (2, 3), (0, 3)]

for x, y in black_stones:
    board[y][x] = 'X'
for x, y in white_stones:
    board[y][x] = 'O'

# Check liberties of white group
white_liberties = get_liberties(board, 1, 2, 'O')
print(f"White group liberties: {white_liberties}")

# Test critical points
critical_points = [(1, 0), (2, 1), (0, 2)]
for x, y in critical_points:
    board_copy = [row[:] for row in board]
    board_copy[y][x] = 'X'
    new_liberties = get_liberties(board_copy, 1, 2, 'O')
    print(f"If Black plays at ({x}, {y}): White group liberties = {new_liberties}")