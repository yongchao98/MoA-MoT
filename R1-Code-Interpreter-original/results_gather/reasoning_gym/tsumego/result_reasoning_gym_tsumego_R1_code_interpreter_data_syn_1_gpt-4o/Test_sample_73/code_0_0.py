# Define the board
board = [
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', 'X', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', 'O', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', 'X', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'X', 'O', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Define potential moves and their impact
potential_moves = {
    'F3': [(3, 4), (3, 5), (3, 6), (3, 7), (3, 8)],  # Group at D4-I4
    'F12': [(2, 5)],  # Single stone at F11
    'C9': [(5, 2)]  # Single stone at C8
}

# Function to count captured stones
def count_captured_stones(move, group):
    captured = 0
    for r, c in group:
        if board[r][c] == 'O':
            captured += 1
    return captured

# Evaluate each move
best_move = None
max_captured = 0

for move, group in potential_moves.items():
    captured = count_captured_stones(move, group)
    if captured > max_captured:
        max_captured = captured
        best_move = move

print(best_move)