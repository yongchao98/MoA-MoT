# Initial matrix
matrix = [
    [0, 16, 32],
    [0, 0, 32],
    [16, 4, 2]
]

# Initial position of the number to move
initial_position = (2, 1)  # (row, column) in 0-indexed

# Move sequence
moves = 'dlr'

# Function to perform the move
def move_number(matrix, position, moves):
    row, col = position
    for move in moves:
        if move == 'd':  # move down
            if row < 2:
                row += 1
        elif move == 'u':  # move up
            if row > 0:
                row -= 1
        elif move == 'l':  # move left
            if col > 0:
                col -= 1
        elif move == 'r':  # move right
            if col < 2:
                col += 1
        
        # Check for combination
        if matrix[row][col] == matrix[position[0]][position[1]]:
            matrix[row][col] *= 2
            matrix[position[0]][position[1]] = 0
        else:
            # Move the number
            matrix[row][col] = matrix[position[0]][position[1]]
            matrix[position[0]][position[1]] = 0
        
        # Update the position
        position = (row, col)
    
    return matrix

# Perform the moves
final_matrix = move_number(matrix, initial_position, moves)

# Print the final matrix
print(final_matrix)