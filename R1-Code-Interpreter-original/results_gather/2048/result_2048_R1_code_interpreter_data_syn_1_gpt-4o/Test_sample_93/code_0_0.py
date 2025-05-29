# Initial matrix
matrix = [
    [0, 0, 0, 0, 2],
    [0, 32, 2, 0, 0],
    [0, 4, 16, 0, 2],
    [32, 0, 32, 2, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number to move
current_position = (2, 2)

# Move sequence
moves = 'duldrl'

# Direction vectors for moves
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Function to perform the move
def move_number(matrix, current_position, moves):
    x, y = current_position
    for move in moves:
        dx, dy = directions[move]
        new_x, new_y = x + dx, y + dy
        # Check if the new position is within bounds
        if 0 <= new_x < 5 and 0 <= new_y < 5:
            # Check if the numbers can combine
            if matrix[new_x][new_y] == matrix[x][y]:
                matrix[new_x][new_y] *= 2
                matrix[x][y] = 0
            else:
                # Move the number to the new position
                matrix[new_x][new_y], matrix[x][y] = matrix[x][y], matrix[new_x][new_y]
            # Update current position
            x, y = new_x, new_y
    return matrix

# Perform the moves
result_matrix = move_number(matrix, current_position, moves)

# Print the final matrix
print(result_matrix)