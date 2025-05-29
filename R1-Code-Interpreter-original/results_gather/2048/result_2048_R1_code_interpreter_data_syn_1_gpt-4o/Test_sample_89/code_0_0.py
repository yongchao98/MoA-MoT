# Initial matrix
matrix = [
    [0, 2, 0, 4, 0],
    [0, 16, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [32, 0, 8, 0, 0],
    [0, 0, 16, 0, 0]
]

# Initial position
position = (1, 3)

# Movement sequence
sequence = 'rurlrr'

# Define the movement directions
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Function to move the number
def move_number(matrix, position, sequence):
    x, y = position
    for move in sequence:
        dx, dy = directions[move]
        new_x, new_y = x + dx, y + dy
        # Check if the new position is within bounds
        if 0 <= new_x < 5 and 0 <= new_y < 5:
            # Check if the numbers can be combined
            if matrix[new_x][new_y] == matrix[x][y] and matrix[x][y] != 0:
                matrix[new_x][new_y] *= 2
                matrix[x][y] = 0
            else:
                # Move the number
                matrix[new_x][new_y], matrix[x][y] = matrix[x][y], matrix[new_x][new_y]
            x, y = new_x, new_y
    return matrix

# Apply the movement
result_matrix = move_number(matrix, position, sequence)

# Print the resulting matrix
print(result_matrix)