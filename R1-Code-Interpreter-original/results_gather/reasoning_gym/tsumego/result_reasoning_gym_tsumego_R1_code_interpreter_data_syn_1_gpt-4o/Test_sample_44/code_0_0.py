# Define the board
board = [
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.', 'O'],
    ['X', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.'],
    ['X', '.', '.', '.', '.', '.', 'X', 'X', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.']
]

# Function to check if a move at (x, y) captures any stones
def check_capture(x, y, board):
    # Check if the move is valid
    if board[x][y] != '.':
        return False
    
    # Place the stone temporarily
    board[x][y] = 'X'
    
    # Check the group of stones around the move
    def get_group(x, y, color):
        group = set()
        stack = [(x, y)]
        while stack:
            cx, cy = stack.pop()
            if (cx, cy) not in group and 0 <= cx < len(board) and 0 <= cy < len(board[0]) and board[cx][cy] == color:
                group.add((cx, cy))
                stack.extend([(cx+1, cy), (cx-1, cy), (cx, cy+1), (cx, cy-1)])
        return group
    
    # Check if the group has any liberties
    def has_liberties(group):
        for gx, gy in group:
            for nx, ny in [(gx+1, gy), (gx-1, gy), (gx, gy+1), (gx, gy-1)]:
                if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == '.':
                    return True
        return False
    
    # Check all adjacent groups of the opposite color
    captured = False
    for nx, ny in [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]:
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == 'O':
            group = get_group(nx, ny, 'O')
            if not has_liberties(group):
                captured = True
                for gx, gy in group:
                    board[gx][gy] = '.'
    
    # Remove the temporary stone
    board[x][y] = '.'
    
    return captured

# Check if playing at F7 captures any stones
capture_result = check_capture(4, 5, board)
print(capture_result)