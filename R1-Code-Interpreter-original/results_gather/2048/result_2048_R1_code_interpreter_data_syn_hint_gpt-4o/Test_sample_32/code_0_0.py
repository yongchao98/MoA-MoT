# Initial matrix
matrix = [
    [0, 0, 0],
    [0, 8, 0],
    [32, 32, 0]
]

# Initial position of the number to move
position = (2, 1)  # (3, 2) in 0-indexed

# Move sequence
moves = 'uur'

# Function to move the number
def move_number(matrix, position, moves):
    for move in moves:
        x, y = position
        if move == 'u':
            new_position = (x - 1, y)
        elif move == 'd':
            new_position = (x + 1, y)
        elif move == 'l':
            new_position = (x, y - 1)
        elif move == 'r':
            new_position = (x, y + 1)
        
        # Check if the new position is within bounds
        if 0 <= new_position[0] < 3 and 0 <= new_position[1] < 3:
            # Move the number
            if matrix[new_position[0]][new_position[1]] == matrix[x][y]:
                # Combine numbers
                matrix[new_position[0]][new_position[1]] *= 2
                matrix[x][y] = 0
            elif matrix[new_position[0]][new_position[1]] == 0:
                # Move to empty space
                matrix[new_position[0]][new_position[1]] = matrix[x][y]
                matrix[x][y] = 0
            position = new_position
    
    return matrix

# Apply the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)