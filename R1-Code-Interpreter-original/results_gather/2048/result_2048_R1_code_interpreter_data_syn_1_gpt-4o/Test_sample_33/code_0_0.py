# Initial matrix
matrix = [
    [8, 0, 4],
    [32, 16, 16],
    [0, 0, 2]
]

# Move sequence
moves = 'lll'

# Function to move left
def move_left(matrix, row):
    # Find the first non-zero element from the right
    for col in range(2, -1, -1):
        if matrix[row][col] != 0:
            # Move it to the leftmost available position
            for target_col in range(3):
                if matrix[row][target_col] == 0:
                    matrix[row][target_col] = matrix[row][col]
                    matrix[row][col] = 0
                    break
                elif matrix[row][target_col] == matrix[row][col]:
                    # Combine if they are the same
                    matrix[row][target_col] *= 2
                    matrix[row][col] = 0
                    break
            break

# Apply the moves
for move in moves:
    if move == 'l':
        move_left(matrix, 2)

# Print the final matrix
print(matrix)