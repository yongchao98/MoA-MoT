# Initial matrix
matrix = [
    [32, 8, 32],
    [0, 32, 16],
    [0, 0, 2]
]

# Move sequence
sequence = 'uul'

# Initial position of the number to move
pos = (1, 1)  # (row, column) in 0-indexed

# Function to move the number in the matrix
def move_number(matrix, pos, sequence):
    for move in sequence:
        r, c = pos
        if move == 'u':
            new_pos = (r - 1, c)
        elif move == 'd':
            new_pos = (r + 1, c)
        elif move == 'l':
            new_pos = (r, c - 1)
        elif move == 'r':
            new_pos = (r, c + 1)
        
        # Check if the new position is within bounds
        if 0 <= new_pos[0] < len(matrix) and 0 <= new_pos[1] < len(matrix[0]):
            # Check if the numbers can combine
            if matrix[new_pos[0]][new_pos[1]] == matrix[r][c]:
                matrix[new_pos[0]][new_pos[1]] *= 2
                matrix[r][c] = 0
            elif matrix[new_pos[0]][new_pos[1]] == 0:
                matrix[new_pos[0]][new_pos[1]] = matrix[r][c]
                matrix[r][c] = 0
            pos = new_pos
    return matrix

# Apply the move sequence
result_matrix = move_number(matrix, pos, sequence)

# Print the resulting matrix
print(result_matrix)