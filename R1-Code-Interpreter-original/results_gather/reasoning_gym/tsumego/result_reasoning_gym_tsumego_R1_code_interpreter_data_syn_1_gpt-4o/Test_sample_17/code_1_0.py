# Define the board
board = [
    ['O', '.', '.', '.', '.', '.', 'O', '.', 'O'],
    ['.', '.', '.', 'X', 'X', '.', '.', '.', '.'],
    ['.', 'X', 'X', 'O', 'O', 'X', '.', '.', '.'],
    ['.', 'O', 'X', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'X', '.', 'X', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', 'O'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def check_capture(board, move):
    x, y = move
    board[x][y] = 'X'  # Place the black stone
    captured = 0

    # Check surrounding groups for capture
    # Define directions for checking
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    def is_captured(x, y, visited):
        if (x, y) in visited or board[x][y] == 'X':
            return True
        if board[x][y] == '.':
            return False
        visited.add((x, y))
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 9 and 0 <= ny < 9 and not is_captured(nx, ny, visited):
                return False
        return True

    # Check each white stone if it can be captured
    for i in range(9):
        for j in range(9):
            if board[i][j] == 'O':
                visited = set()
                if is_captured(i, j, visited):
                    captured += len(visited)
                    for vx, vy in visited:
                        board[vx][vy] = '.'  # Remove captured stones

    board[x][y] = '.'  # Reset the move
    return captured

# Evaluate potential moves
# Consider more potential moves around the white stones
potential_moves = [(2, 3), (5, 6), (3, 5), (1, 5), (4, 3), (3, 4), (4, 6), (5, 4)]
best_move = None
max_captured = 0

for move in potential_moves:
    captured = check_capture(board, move)
    if captured > max_captured:
        max_captured = captured
        best_move = move

# Convert the best move to board coordinates
if best_move is not None:
    move_coords = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    best_move_coord = f"{move_coords[best_move[1]]}{9 - best_move[0]}"
    print(best_move_coord)
else:
    print("No capturing move found")