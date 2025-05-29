# Define the board
board = [
    ['.', '.', '.', '.', 'O', '.', '.', 'X', 'O', 'X'],
    ['X', '.', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', 'O', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'X'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check liberties of a group
def check_liberties(x, y, color):
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

# Check liberties for the stone at C4
liberties_before = check_liberties(6, 1, 'O')

# Simulate placing a black stone at C3
board[7][1] = 'X'

# Check liberties again
liberties_after = check_liberties(6, 1, 'O')

# Output the result
print("Liberties before:", liberties_before)
print("Liberties after:", liberties_after)