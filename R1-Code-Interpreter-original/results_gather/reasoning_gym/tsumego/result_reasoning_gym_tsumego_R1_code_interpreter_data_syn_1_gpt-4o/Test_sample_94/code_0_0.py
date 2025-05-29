# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', 'X', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['O', 'O', '.', 'X', '.', '.', 'X', '.', 'O', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.']
]

# Function to check if a group is captured
def is_captured(board, x, y):
    # Directions for adjacent positions
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    visited = set()
    stack = [(x, y)]
    liberties = 0

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties == 0

# Check if placing a black stone at C4 captures the group
board[8][2] = 'X'  # Place black stone at C4
captured = is_captured(board, 6, 1)  # Check the group starting at B7

print(captured)