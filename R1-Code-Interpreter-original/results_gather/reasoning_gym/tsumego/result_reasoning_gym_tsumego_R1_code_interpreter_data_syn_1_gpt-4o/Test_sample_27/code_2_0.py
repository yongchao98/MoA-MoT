# Define the board as a 2D list
board = [
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', 'X']
]

# Function to check liberties of a group
def check_liberties(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        # Check all 4 directions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties

# Simulate the move D5
board[5][3] = 'X'  # D5

# Check liberties of the White group starting from D4
liberties = check_liberties(board, 6, 2, 'O')
print(liberties)