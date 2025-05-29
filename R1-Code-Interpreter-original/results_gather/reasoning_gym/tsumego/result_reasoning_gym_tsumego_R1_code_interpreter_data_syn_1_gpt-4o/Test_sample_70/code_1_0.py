# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.'],
    ['.', '.', 'O', 'O', '.', '.', 'X', 'O', '.', '.'],
    ['.', '.', 'X', '.', '.', 'X', 'O', 'X', '.', '.']
]

# Function to check if a move at (x, y) captures any stones
def check_capture(x, y):
    # Place the black stone
    board[x][y] = 'X'
    
    # Check surrounding white groups for capture
    captured = False
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    def is_captured(x, y, visited):
        if (x, y) in visited:
            return True
        if x < 0 or x >= 10 or y < 0 or y >= 10:
            return False
        if board[x][y] == '.':
            return False
        if board[x][y] == 'X':
            return True
        
        visited.add((x, y))
        return all(is_captured(x + dx, y + dy, visited) for dx, dy in directions)
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and board[nx][ny] == 'O':
            if is_captured(nx, ny, set()):
                captured = True
                # Capture the stones
                def capture(x, y):
                    if board[x][y] == 'O':
                        board[x][y] = '.'
                        for dx, dy in directions:
                            capture(x + dx, y + dy)
                capture(nx, ny)
    
    return captured

# Check if placing a stone at G6 captures any stones
captured = check_capture(4, 6)  # G6 corresponds to (4, 6) in 0-indexed array
print(captured)