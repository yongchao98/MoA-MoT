# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', 'O', '.', '.'],
    ['.', 'O', 'X', 'X', '.', '.', 'O', '.', '.'],
    ['.', 'X', 'O', 'O', 'X', 'O', '.', 'O', '.'],
    ['.', 'X', 'O', '.', '.', '.', 'O', '.', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'X', '.', '.', '.', 'X'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a group of stones is captured
def has_liberties(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        # Check all adjacent positions
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[nx][ny] == '.':
                    return True
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return False

# Simulate the move at E5
board[4][4] = 'X'

# Check if the group at E6, F6, and F5 has liberties
captured = not has_liberties(board, 3, 4, 'O')

print(captured)