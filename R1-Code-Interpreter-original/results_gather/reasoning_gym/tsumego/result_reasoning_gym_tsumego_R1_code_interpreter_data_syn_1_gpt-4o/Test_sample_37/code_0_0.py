# Define the board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', '.', 'X', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', 'O', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', 'O', '.', 'O'],
    ['.', 'X', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'X', '.', '.']
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
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Place a black stone at C7
board[2][2] = 'X'

# Check liberties of the white group at B8, B7, B6
liberties = check_liberties(board, 1, 1, 'O')
print("Liberties of the white group at B8, B7, B6 after playing C7:", liberties)