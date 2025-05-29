# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'O', 'X'],
    ['X', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', 'X', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', 'O', '.', '.', '.'],
    ['X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def check_capture(board, move):
    x, y = move
    # Place the black stone
    board[y][x] = 'X'
    
    # Check the group of white stones around the move
    def get_group(x, y, color):
        group = set()
        stack = [(x, y)]
        while stack:
            cx, cy = stack.pop()
            if (cx, cy) not in group and 0 <= cx < 9 and 0 <= cy < 9 and board[cy][cx] == color:
                group.add((cx, cy))
                stack.extend([(cx+1, cy), (cx-1, cy), (cx, cy+1), (cx, cy-1)])
        return group
    
    # Check liberties of a group
    def has_liberty(group):
        for gx, gy in group:
            for nx, ny in [(gx+1, gy), (gx-1, gy), (gx, gy+1), (gx, gy-1)]:
                if 0 <= nx < 9 and 0 <= ny < 9 and board[ny][nx] == '.':
                    return True
        return False
    
    # Check the group of white stones at G7
    group = get_group(6, 6, 'O')
    if not has_liberty(group):
        return True
    return False

# Check if the move at F7 captures stones
move = (5, 6)  # F7
capture = check_capture(board, move)
print(capture)