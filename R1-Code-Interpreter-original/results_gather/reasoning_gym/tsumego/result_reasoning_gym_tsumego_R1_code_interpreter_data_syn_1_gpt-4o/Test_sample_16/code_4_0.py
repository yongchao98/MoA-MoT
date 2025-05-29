# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', 'O', 'O', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.']
]

# Function to check if a move at (x, y) captures any stones
def check_capture(x, y):
    # Place the black stone
    board[y][x] = 'X'
    
    # Check if any white stones are captured
    def has_liberty(x, y, visited):
        if x < 0 or x >= 13 or y < 0 or y >= 13:
            return False
        if board[y][x] == '.':
            return True
        if board[y][x] == 'X' or (x, y) in visited:
            return False
        visited.add((x, y))
        return (has_liberty(x+1, y, visited) or has_liberty(x-1, y, visited) or
                has_liberty(x, y+1, visited) or has_liberty(x, y-1, visited))
    
    captured = False
    for i in range(13):
        for j in range(13):
            if board[j][i] == 'O':
                visited = set()
                if not has_liberty(i, j, visited):
                    captured = True
    
    # Remove the black stone
    board[y][x] = '.'
    return captured

# Check if placing a stone at E8 captures any stones
capture_result = check_capture(4, 5)
print(capture_result)