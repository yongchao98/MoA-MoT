# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', 'O', '.', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', 'X', 'X', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.']
]

# Function to check if a group is captured
def is_captured(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties == 0

# Simulate the move at B8
board[1][8] = 'X'  # Place black stone at B8

# Check if the group at B9, C9 is captured
captured = is_captured(board, 1, 9, 'O')

# Output the result
print("Captured:", captured)