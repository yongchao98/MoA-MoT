# Initial matrix
matrix = [
    [0, 16, 32, 0],
    [0, 32, 8, 0],
    [0, 8, 0, 0],
    [0, 0, 0, 8]
]

# Initial position of the number to move
position = (1, 1)  # (row, column) in 0-indexed

# Move sequence
moves = 'urdr'

# Function to perform the move
def move_number(matrix, position, moves):
    for move in moves:
        row, col = position
        if move == 'u':
            new_row, new_col = row - 1, col
        elif move == 'd':
            new_row, new_col = row + 1, col
        elif move == 'l':
            new_row, new_col = row, col - 1
        elif move == 'r':
            new_row, new_col = row, col + 1
        else:
            continue
        
        # Check if the new position is within bounds
        if 0 <= new_row < len(matrix) and 0 <= new_col < len(matrix[0]):
            # Check if the numbers can combine
            if matrix[new_row][new_col] == matrix[row][col]:
                matrix[new_row][new_col] *= 2
                matrix[row][col] = 0
            elif matrix[new_row][new_col] == 0:
                # Move the number to the new position
                matrix[new_row][new_col] = matrix[row][col]
                matrix[row][col] = 0
            # Update the position
            position = (new_row, new_col)
    
    return matrix

# Perform the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)