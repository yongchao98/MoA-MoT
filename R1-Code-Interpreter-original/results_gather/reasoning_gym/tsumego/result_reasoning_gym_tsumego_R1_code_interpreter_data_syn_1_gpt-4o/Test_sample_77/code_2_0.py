# Define the board state
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', 'O', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'O', '.', '.', 'X', '.', '.', '.'],
    ['.', 'O', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def is_capture(board, move, player):
    # Simulate the move
    x, y = move
    board[y][x] = player
    
    # Check for captures
    def get_neighbors(x, y):
        neighbors = []
        if x > 0: neighbors.append((x-1, y))
        if x < 10: neighbors.append((x+1, y))
        if y > 0: neighbors.append((x, y-1))
        if y < 10: neighbors.append((x, y+1))
        return neighbors
    
    def count_liberties(x, y, visited):
        if (x, y) in visited:
            return 0
        visited.add((x, y))
        liberties = 0
        for nx, ny in get_neighbors(x, y):
            if board[ny][nx] == '.':
                liberties += 1
            elif board[ny][nx] == board[y][x]:
                liberties += count_liberties(nx, ny, visited)
        return liberties
    
    # Check if the move captures any stones
    captured = False
    for nx, ny in get_neighbors(x, y):
        if board[ny][nx] != player and board[ny][nx] != '.':
            if count_liberties(nx, ny, set()) == 0:
                captured = True
                break
    
    # Revert the move
    board[y][x] = '.'
    return captured

# Check if C5 is a capturing move
move = (2, 5)  # C5
player = 'X'
capture = is_capture(board, move, player)
print(capture)