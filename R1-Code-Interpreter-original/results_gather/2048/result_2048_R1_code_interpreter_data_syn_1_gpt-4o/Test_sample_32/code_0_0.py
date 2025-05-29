# Initial matrix
matrix = [
    [0, 0, 0],
    [0, 8, 0],
    [32, 32, 0]
]

# Move sequence
moves = 'uur'

# Initial position of the number to move
position = (2, 1)  # (row, column) in 0-indexed

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    for move in moves:
        row, col = position
        if move == 'u':
            new_position = (row - 1, col)
        elif move == 'd':
            new_position = (row + 1, col)
        elif move == 'l':
            new_position = (row, col - 1)
        elif move == 'r':
            new_position = (row, col + 1)
        
        # Move the number
        matrix[new_position[0]][new_position[1]] = matrix[row][col]
        matrix[row][col] = 0
        position = new_position
    
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)