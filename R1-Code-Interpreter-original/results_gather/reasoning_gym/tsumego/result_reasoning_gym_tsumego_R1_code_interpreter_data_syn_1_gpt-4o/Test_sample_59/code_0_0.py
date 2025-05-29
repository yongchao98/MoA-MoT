# Simulate the move at C4 and check if it captures any white stones
def check_capture(board, move):
    # Place the black stone
    board[move[0]][move[1]] = 'X'
    
    # Check if any white stones are captured
    def is_captured(x, y, visited):
        if (x, y) in visited:
            return True
        visited.add((x, y))
        
        # Check all adjacent positions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    return False
                if board[nx][ny] == 'O' and not is_captured(nx, ny, visited):
                    return False
        return True
    
    captured = False
    for x in range(len(board)):
        for y in range(len(board[0])):
            if board[x][y] == 'O' and is_captured(x, y, set()):
                captured = True
                board[x][y] = '.'
    
    return captured

# Initial board setup
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'O', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', 'O', 'X', '.'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', 'X', '.'],
    ['X', 'O', 'X', 'O', '.', '.', 'X', '.', '.']
]

# Check if playing at C4 captures any stones
move = (5, 2)  # C4 in 0-indexed coordinates
captured = check_capture(board, move)
print(captured)