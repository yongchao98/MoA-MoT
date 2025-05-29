# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', 'X', 'O', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', 'O', 'X', 'O', 'O', 'O', 'X', '.'],
    ['.', '.', 'X', '.', 'X', 'O', 'X', '.', '.'],
    ['.', 'O', 'X', 'X', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.']
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
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Place a black stone at D5
board[4][3] = 'X'

# Check if the white group at D6, E6, F6, E5 is captured
captured = is_captured(board, 3, 2, 'O')

print(captured)