# Define the board
board = [
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', 'O', '.', '.'],
    ['.', 'X', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
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
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[ny][nx] == '.':
                    liberties += 1
                elif board[ny][nx] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Check liberties for the vulnerable groups
liberties_C8 = check_liberties(board, 2, 8, 'O')
liberties_F5 = check_liberties(board, 5, 5, 'O')
liberties_H4 = check_liberties(board, 7, 4, 'O')

# Determine the best move
best_move = None
if liberties_C8 == 1:
    best_move = 'D8'  # Place at D8 to capture the group
elif liberties_F5 == 1:
    best_move = 'E5'  # Place at E5 to capture the group
elif liberties_H4 == 1:
    best_move = 'G4'  # Place at G4 to capture the group

print(best_move)