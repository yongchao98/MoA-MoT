# Initial matrix
matrix = [
    [0, 0, 0, 8, 0],
    [0, 4, 0, 0, 0],
    [0, 0, 2, 0, 8],
    [0, 0, 0, 0, 0],
    [0, 16, 0, 0, 0]
]

# Initial position of the number 16
position = (4, 1)

# Move sequence
moves = 'lrrul'

# Move deltas
move_deltas = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Simulate the moves
for move in moves:
    delta = move_deltas[move]
    new_position = (position[0] + delta[0], position[1] + delta[1])
    
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        # Move the number
        position = new_position

# Update the matrix
matrix[4][1] = 0  # Clear the initial position
matrix[position[0]][position[1]] = 16  # Set the new position

print(matrix)