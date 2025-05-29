# Define the board state
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', 'X', '.', '.', '.', 'O', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'X'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', 'O', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
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
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Check the liberties of the White group at G9, H9, G8
liberties_before = check_liberties(board, 8, 6, 'O')

# Simulate the move at G10
board[9][6] = 'X'

# Check the liberties after the move
liberties_after = check_liberties(board, 8, 6, 'O')

# Output the result
print("Liberties before move:", liberties_before)
print("Liberties after move:", liberties_after)