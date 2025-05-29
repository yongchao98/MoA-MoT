# Initial matrix
matrix = [
    [0, 0, 8, 32, 0],
    [8, 0, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [4, 4, 0, 32, 4],
    [2, 32, 0, 0, 0]
]

# Initial position of the number to move
position = (4, 1)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'lrrur'

# Function to perform a move
def move(position, direction):
    row, col = position
    if direction == 'l':
        col -= 1
    elif direction == 'r':
        col += 1
    elif direction == 'u':
        row -= 1
    elif direction == 'd':
        row += 1
    return row, col

# Simulate the moves
for move_direction in moves:
    new_position = move(position, move_direction)
    new_row, new_col = new_position
    current_value = matrix[position[0]][position[1]]
    target_value = matrix[new_row][new_col]
    
    # Check if we can combine
    if current_value == target_value:
        matrix[new_row][new_col] *= 2
        matrix[position[0]][position[1]] = 0
    else:
        matrix[new_row][new_col] = current_value
        matrix[position[0]][position[1]] = 0
    
    # Update the position
    position = new_position

# Print the final matrix
print(matrix)