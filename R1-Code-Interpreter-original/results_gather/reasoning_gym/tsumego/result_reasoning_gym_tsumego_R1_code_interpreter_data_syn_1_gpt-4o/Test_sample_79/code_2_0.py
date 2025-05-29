# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', 'X'],
    ['O', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.']
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
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Simulate placing a black stone at G8 (index 7, 6)
board[7][6] = 'X'

# Check if the stone at G7 is captured
captured = is_captured(board, 6, 6, 'O')

print(captured)