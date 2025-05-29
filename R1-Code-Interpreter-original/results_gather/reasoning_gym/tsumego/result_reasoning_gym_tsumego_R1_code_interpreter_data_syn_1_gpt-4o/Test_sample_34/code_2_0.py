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
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Simulate the move at E5
board[4][4] = 'X'

# Check if the group at E6, F6, and F5 is captured
captured = is_captured(board, 3, 4, 'O') and is_captured(board, 3, 5, 'O') and is_captured(board, 4, 5, 'O')

print(captured)