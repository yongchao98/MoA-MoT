# Define the board
board = [
    ['.', '.', '.', '.', 'O', '.', 'O', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', 'X'],
    ['.', 'O', '.', '.', '.', 'X', 'O', 'O', 'O'],
    ['.', '.', '.', 'X', '.', 'O', 'X', 'O', 'X'],
    ['.', 'O', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', 'O', '.', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.']
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

# Check the liberties of the group at F5
liberties_before = check_liberties(board, 4, 5, 'O')

# Simulate the move at F6
board[3][5] = 'X'

# Check the liberties of the group at F5 after the move
liberties_after = check_liberties(board, 4, 5, 'O')

# Determine if the move captures any stones
captured = liberties_after == 0

print("Liberties before move:", liberties_before)
print("Liberties after move:", liberties_after)
print("Captured:", captured)