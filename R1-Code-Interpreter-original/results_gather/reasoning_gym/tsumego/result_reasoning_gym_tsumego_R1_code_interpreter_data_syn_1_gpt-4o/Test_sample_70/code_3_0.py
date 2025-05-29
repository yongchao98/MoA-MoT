# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.'],
    ['.', '.', 'O', 'O', '.', '.', 'X', 'O', '.', '.'],
    ['.', '.', 'X', '.', '.', 'X', 'O', 'X', '.', '.']
]

# Function to check liberties of a group of stones
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
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Place a black stone at H2 (index 8, 7)
board[8][7] = 'X'

# Check liberties of the white stones at G3 and H3
liberties_G3 = check_liberties(board, 7, 6, 'O')
liberties_H3 = check_liberties(board, 7, 7, 'O')

print(liberties_G3, liberties_H3)