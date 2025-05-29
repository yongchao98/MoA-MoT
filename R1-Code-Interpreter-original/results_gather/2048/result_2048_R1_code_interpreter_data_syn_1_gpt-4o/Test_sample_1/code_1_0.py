# Initial matrix
matrix = [
    [8, 32],
    [0, 16]
]

# Move sequence
moves = 'dr'

# Initial position of the number to move
position = (0, 0)  # Adjusted to 0-based index

# Function to perform the moves
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'd':
            x += 1
        elif move == 'r':
            y += 1
        elif move == 'u':
            x -= 1
        elif move == 'l':
            y -= 1
        
        # Check if the move is within bounds
        if 0 <= x < 2 and 0 <= y < 2:
            # Move the number
            matrix[x][y], matrix[position[0]][position[1]] = matrix[position[0]][position[1]], 0
            position = (x, y)
    
    return matrix

# Perform the moves
result_matrix = move_number(matrix, position, moves)

# Print the result
print(result_matrix)