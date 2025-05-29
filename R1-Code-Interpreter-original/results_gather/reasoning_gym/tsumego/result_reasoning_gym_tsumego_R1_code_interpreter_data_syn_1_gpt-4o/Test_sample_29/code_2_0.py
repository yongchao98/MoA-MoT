# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'X', 'X', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', 'X', 'X', '.', 'X', '.', '.', '.', '.', 'X', 'X', '.'],
    ['.', '.', 'X', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.']
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

# Place a black stone at J9
board[3][9] = 'X'

# Check liberties of the white group at H9
liberties_after_move = check_liberties(board, 3, 8)

print(liberties_after_move)