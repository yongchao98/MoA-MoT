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

# Check liberties of the White group at H4, H3, G3 after Black plays at G4
board[7][6] = 'X'  # Black plays at G4
liberties_H4 = check_liberties(board, 7, 7)
liberties_H3 = check_liberties(board, 8, 7)
liberties_G3 = check_liberties(board, 8, 6)

print("Liberties of White group at H4:", liberties_H4)
print("Liberties of White group at H3:", liberties_H3)
print("Liberties of White group at G3:", liberties_G3)