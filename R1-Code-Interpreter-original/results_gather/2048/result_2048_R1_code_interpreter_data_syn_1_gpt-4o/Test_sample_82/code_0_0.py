# Initial matrix
matrix = [
    [0, 8, 16, 2, 2],
    [2, 8, 0, 0, 0],
    [0, 4, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [8, 0, 8, 0, 0]
]

# Initial position of the number to move
position = (1, 3)

# Sequence of moves
moves = 'uldlur'

# Define the movement directions
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Simulate the movement
for move in moves:
    # Calculate new position
    new_position = (position[0] + directions[move][0], position[1] + directions[move][1])
    
    # Check if the new position is within bounds
    if 0 <= new_position[0] < 5 and 0 <= new_position[1] < 5:
        # Move the number
        matrix[new_position[0]][new_position[1]] = matrix[position[0]][position[1]]
        matrix[position[0]][position[1]] = 0
        position = new_position

# Print the final matrix
print(matrix)