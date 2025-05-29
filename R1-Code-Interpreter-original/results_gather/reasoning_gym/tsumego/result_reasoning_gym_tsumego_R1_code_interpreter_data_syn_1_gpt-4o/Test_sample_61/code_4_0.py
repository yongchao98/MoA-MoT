# Define the board state
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', 'X', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', 'O', '.', 'O', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
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
            if 0 <= nx < 12 and 0 <= ny < 12:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Check liberties of the White group at F6 (5, 6) after Black plays at G6 (5, 7)
board[5][7] = 'X'  # Black plays at G6
liberties_after_move = check_liberties(board, 5, 6)
print(liberties_after_move)