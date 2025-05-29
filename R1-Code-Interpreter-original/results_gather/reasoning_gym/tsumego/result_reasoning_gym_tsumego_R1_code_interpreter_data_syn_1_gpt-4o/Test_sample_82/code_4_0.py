# Define the board
board = [
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', 'X', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.']
]

# Function to check if a move captures stones
def captures_stones(board, move):
    x, y = move
    # Check the group at (C8, D8, C9)
    group = [(2, 8), (3, 8), (2, 9)]
    liberties = [(3, 7)]
    if (x, y) in liberties:
        return all(board[gy][gx] == 'O' for gx, gy in group)
    return False

# Check if playing at D7 captures the group
move = (3, 7)  # D7
result = captures_stones(board, move)
print(result)