# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', 'X', '.', '.', '.', '.', 'X', '.', '.'],
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.']
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

# Place a Black stone at B9
board[0][1] = 'X'

# Check liberties of the White group starting at C9
liberties = check_liberties(board, 0, 2, 'O')
print(liberties)