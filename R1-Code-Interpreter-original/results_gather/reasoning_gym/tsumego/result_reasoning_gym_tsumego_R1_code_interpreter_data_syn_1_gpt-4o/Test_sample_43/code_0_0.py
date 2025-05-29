# Define the board
board = [
    ['.', 'X', 'O', 'X', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', 'X', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', 'X', 'O', 'X', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', 'O', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', 'X', '.']
]

# Function to check if a group of stones is captured
def is_captured(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        # Check all adjacent positions
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 12 and 0 <= ny < 12:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Check potential moves
potential_moves = [(8, 2), (9, 5), (11, 7)]
captures = []

for move in potential_moves:
    x, y = move
    board[x][y] = 'X'  # Simulate the move
    if is_captured(board, 2, 2, 'O') or is_captured(board, 5, 5, 'O') or is_captured(board, 7, 7, 'O'):
        captures.append(move)
    board[x][y] = '.'  # Reset the move

# Output the first capturing move found
if captures:
    move = captures[0]
    print(f"<<<{chr(move[1] + 65)}{12 - move[0]}>>>")
else:
    print("No capturing move found.")