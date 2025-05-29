# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', 'O', '.', '.', '.', 'X', '.'],
    ['.', 'X', '.', '.', '.', '.', 'O', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'X', '.', 'X', '.', '.']
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
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Check potential moves
potential_moves = [(11, 5), (2, 5)]
for move in potential_moves:
    x, y = move
    board[x][y] = 'X'  # Place a black stone
    if is_captured(board, 11, 6, 'O') or is_captured(board, 2, 6, 'O'):
        print(f"Move at {chr(y + 65)}{13 - x} captures white stones.")
    board[x][y] = '.'  # Reset the board
