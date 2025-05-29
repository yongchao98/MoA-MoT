# Define the board
board = [
    ['.', '.', '.', 'X', '.', '.', 'O', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', 'O', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', 'O', 'X', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.']
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
                elif board[nx][ny] == board[x][y]:
                    stack.append((nx, ny))
    return liberties

# Check liberties of the White stone at G5 after Black plays at G6
board[5][6] = 'X'  # Black plays at G6
liberties_G5 = check_liberties(board, 6, 6)

print("Liberties of White stone at G5:", liberties_G5)