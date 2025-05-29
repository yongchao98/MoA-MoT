# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', 'O', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', 'X', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', 'X', 'O', 'O', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', 'X', 'O', '.', '.', '.']
]

# Function to check liberties of a group
def check_liberties(board, x, y):
    visited = set()
    stack = [(x, y)]
    liberties = 0
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Place a Black stone at K4 (index 9, 10)
board[9][10] = 'X'

# Check liberties of the White group around J3, K3, L3
liberties = check_liberties(board, 10, 9)  # Starting from J3 (index 10, 9)
print(liberties)