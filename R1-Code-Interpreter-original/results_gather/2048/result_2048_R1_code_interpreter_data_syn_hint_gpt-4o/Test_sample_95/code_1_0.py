# Initial matrix
matrix = [
    [4, 32, 0, 16, 0],
    [16, 0, 0, 16, 0],
    [16, 4, 16, 0, 0],
    [0, 0, 32, 16, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number to move
start_pos = (2, 2)  # Adjusted to 0-based index
number = matrix[start_pos[0]][start_pos[1]]

# Move sequence
moves = 'rrdluu'

# Directions mapping
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Current position
current_pos = start_pos

# Clear the initial position
matrix[current_pos[0]][current_pos[1]] = 0

# Process each move
for move in moves:
    dx, dy = directions[move]
    new_pos = (current_pos[0] + dx, current_pos[1] + dy)
    
    # Check if the new position is within bounds
    if 0 <= new_pos[0] < 5 and 0 <= new_pos[1] < 5:
        # Move the number
        if matrix[new_pos[0]][new_pos[1]] == 0:
            current_pos = new_pos
        elif matrix[new_pos[0]][new_pos[1]] == number:
            # Combine numbers
            number *= 2
            current_pos = new_pos
        else:
            # If different number, stop moving
            break

# Place the number in the final position
matrix[current_pos[0]][current_pos[1]] = number

# Print the final matrix
print(matrix)