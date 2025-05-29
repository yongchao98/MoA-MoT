def is_captured(board, x, y, color):
    """Check if a group of stones is captured."""
    visited = set()
    stack = [(x, y)]
    liberties = 0
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 12 and 0 <= ny < 12:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties == 0

def simulate_move(board, move, color):
    """Simulate a move and return the number of stones captured."""
    x, y = move
    board[x][y] = color
    captured = 0
    opponent_color = 'O' if color == 'X' else 'X'
    for nx, ny in [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]:
        if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == opponent_color:
            if is_captured(board, nx, ny, opponent_color):
                captured += 1
    board[x][y] = '.'
    return captured

def find_best_move(board):
    """Find the best move for Black to capture the most White stones."""
    best_move = None
    max_captured = 0
    for x in range(12):
        for y in range(12):
            if board[x][y] == '.':
                captured = simulate_move(board, (x, y), 'X')
                if captured > max_captured:
                    max_captured = captured
                    best_move = (x, y)
    return best_move

# Initial board setup
board = [
    ['O', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'X', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', 'X', 'O', '.', '.', '.', '.', '.']
]

best_move = find_best_move(board)
# Convert to board coordinates
columns = "ABCDEFGHIJKL"
rows = "123456789101112"
move_str = f"{columns[best_move[1]]}{rows[11-best_move[0]]}"
print(move_str)