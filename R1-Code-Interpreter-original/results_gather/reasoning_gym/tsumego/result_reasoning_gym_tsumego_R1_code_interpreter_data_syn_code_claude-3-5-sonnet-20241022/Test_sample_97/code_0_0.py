def has_liberty(board, x, y):
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 13 and 0 <= ny < 13:  # Check board boundaries
            if board[ny][nx] == '.':
                return True
    return False

def can_capture(board, x, y):
    # Try placing a black stone at (x,y) and check if it captures any white stones
    if board[y][x] != '.':
        return False
        
    # Temporarily place black stone
    board[y][x] = 'X'
    
    # Check adjacent positions for white stones that might be captured
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 13 and 0 <= ny < 13:
            if board[ny][nx] == 'O' and not has_liberty(board, nx, ny):
                board[y][x] = '.'  # Reset the position
                return True
                
    board[y][x] = '.'  # Reset the position
    return False

# Create the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'O', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.', '.']
]

# Check all empty positions for capturing moves
for y in range(13):
    for x in range(13):
        if board[y][x] == '.' and can_capture(board, x, y):
            print(f"Capturing move found at column {chr(65+x)}, row {13-y}")