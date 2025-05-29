# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', 'O', 'O', 'O', '.'],
    ['.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', '.', '.', '.', '.', '.', 'X', '.', '.']
]

# Function to count liberties of a group
def count_liberties(board, x, y, color):
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

# Make the move C6
board[4][2] = 'X'

# Count liberties of the White group starting from D7
liberties_after_move = count_liberties(board, 3, 3, 'O')
print(liberties_after_move)