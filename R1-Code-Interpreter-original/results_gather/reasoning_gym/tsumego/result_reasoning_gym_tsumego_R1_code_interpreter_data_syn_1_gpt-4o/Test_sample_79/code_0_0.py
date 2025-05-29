# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', 'X'],
    ['O', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.']
]

# Define potential moves and the stones they capture
potential_moves = {
    'F7': 2,  # Captures G7, H7
    'I4': 2,  # Captures H4, J4
    'G2': 3   # Captures H2, I2, J2
}

# Find the move that captures the most stones
best_move = max(potential_moves, key=potential_moves.get)
print(best_move)