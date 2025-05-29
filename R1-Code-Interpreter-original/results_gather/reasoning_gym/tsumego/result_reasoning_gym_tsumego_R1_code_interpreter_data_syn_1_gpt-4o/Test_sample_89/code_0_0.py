# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', 'X'],
    ['.', '.', 'X', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', 'X', '.', '.', '.', '.'],
    ['O', 'X', 'X', 'O', 'O', 'X', '.', '.', '.'],
    ['.', '.', 'X', 'O', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.']
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

# Simulate the move E2 for Black
board[7][4] = 'X'

# Check if the White stones at D3, E3, and D2 are captured
captured = all(is_captured(board, x, y, 'O') for x, y in [(6, 3), (6, 4), (7, 3)])

print(captured)