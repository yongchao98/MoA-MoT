# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', 'O', '.', '.'],
    ['.', 'O', 'O', 'X', '.', '.', '.', '.', '.'],
    ['O', 'X', 'O', 'O', 'X', '.', '.', '.', '.'],
    ['X', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', 'O', '.'],
    ['.', '.', '.', 'O', 'X', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.']
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

        # Check all adjacent positions
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Place a black stone at D5 (index 4, 3)
board[4][3] = 'X'

# Check if the white group at C6 (index 3, 2) is captured
captured = is_captured(board, 3, 2, 'O')

print(captured)